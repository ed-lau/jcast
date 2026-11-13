# -*- coding: utf-8 -*-

"""Packaging check for PR 10 (reports/2026-11-13-packaging-check.md): JCAST builds, installs and passes its tests on
each Python version, at the newest dependency versions and at the lowest that its bounds allow.

    build       the sdist and the wheel, with `uv build`, from the files git tracks: the working tree's, or a
                commit's (--ref). The wheel must hold only the jcast package; the sdist must hold every tracked
                file under tests/, so that the tests run from it.
    install     for each interpreter (--python), a fresh virtual environment with the wheel and its `test` extra:
                - newest: the newest versions within the bounds;
                - lowest: the lowest version of each dependency within the bounds (`--resolution lowest-direct`),
                  from wheels only for the packages that need compiling (--only_binary), so that each is the
                  oldest release with a wheel for that Python; the other packages at their newest.
    run         `jcast --help`, `python -m jcast --help` and the version, then the tests, from the sdist with its
                jcast/ removed, so that they import the installed wheel.

Writes <out>/check_<label>.txt: the files in each package, then per environment the installed versions and the
test result.

Usage (from the repository root):
    python reports/packaging/check_install.py --uv <uv> --python <python> [--python <python> ...] \\
        [--resolution newest|lowest] [--ref <commit>] [--jobs N] [--pytest_args "-m 'not network'"] \\
        [--work <dir>] [--out <dir>] [--label <name>]
"""

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# the dependencies, and the packages that need compiling: only their wheels are used for the lowest versions
SHOWN = ('biopython', 'gtfparse', 'pandas', 'numpy', 'scipy', 'requests', 'platformdirs', 'tqdm', 'matplotlib',
         'scikit-learn', 'pytest')
COMPILED = ('biopython', 'pandas', 'numpy', 'scipy', 'matplotlib', 'scikit-learn', 'contourpy', 'kiwisolver',
            'pillow')


def sh(cmd, **kwargs):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kwargs)


def build(uv, python, ref, work):
    """ the sdist and wheel from the tracked files, and the list of tracked files under tests/ """
    src = work / 'src'
    shutil.rmtree(src, ignore_errors=True)
    src.mkdir(parents=True)
    if ref:
        archive = subprocess.run(['git', 'archive', ref], cwd=REPO, capture_output=True, check=True).stdout
        subprocess.run(['tar', '-x', '-C', src], input=archive, check=True)
        tracked = sh(['git', 'ls-tree', '-r', '--name-only', ref, 'tests'], cwd=REPO, check=True).stdout.split()
    else:
        files = sh(['git', 'ls-files', '-z'], cwd=REPO, check=True).stdout.split('\0')
        for f in filter(None, files):
            (src / f).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO / f, src / f)
        tracked = [f for f in files if f.startswith('tests/')]
    dist = work / 'dist'
    shutil.rmtree(dist, ignore_errors=True)
    done = sh([uv, 'build', '--python', python, '--out-dir', dist, src])
    if done.returncode:
        sys.exit(done.stderr)
    return next(dist.glob('*.tar.gz')), next(dist.glob('*.whl')), sorted(tracked)


def check_packages(sdist, wheel, tracked):
    lines = []
    names = zipfile.ZipFile(wheel).namelist()
    top = sorted({n.split('/')[0] for n in names})
    lines.append(f'{wheel.name}: {wheel.stat().st_size / 1e3:.0f} kB, {len(names)} files; top level: {", ".join(top)}')
    with tarfile.open(sdist) as t:
        members = [m.name.split('/', 1)[1] for m in t.getmembers() if m.isfile()]
    missing = [f for f in tracked if f not in members]
    tests = [m for m in members if m.startswith('tests/')]
    lines.append(f'{sdist.name}: {sdist.stat().st_size / 1e6:.1f} MB, {len(members)} files, of which {len(tests)} '
                 f'under tests/; tracked files under tests/ missing from it: {len(missing)}')
    lines += [f'  missing: {f}' for f in missing]
    extra = sorted(set(tests) - set(tracked))
    lines += [f'  untracked file in the sdist: {f}' for f in extra]
    return lines


def environment(uv, python, resolution, sdist, wheel, work, pytest_args):
    """ a fresh virtual environment with the wheel; its versions, the command line and the tests """
    version = sh([python, '-c', 'import sys, platform; print(sys.version.split()[0], platform.machine())']).stdout
    name = f'{resolution}-{version.split()[0]}'
    venv = work / f'venv-{name}'
    shutil.rmtree(venv, ignore_errors=True)
    sh([uv, 'venv', '--python', python, venv], check=True)
    py = venv / 'bin' / 'python'
    if resolution == 'lowest':
        binary = [a for p in COMPILED for a in ('--only-binary', p)]
        installed = sh([uv, 'pip', 'install', '--python', py, '-r', REPO / 'pyproject.toml', '--extra', 'test',
                        '--resolution', 'lowest-direct'] + binary)
        if not installed.returncode:
            installed = sh([uv, 'pip', 'install', '--python', py, '--no-deps', wheel])
    else:
        installed = sh([uv, 'pip', 'install', '--python', py, f'jcast[test] @ {wheel.as_uri()}'])
    lines = [f'== {resolution}, Python {version.strip()}']
    if installed.returncode:
        return name, lines + ['  install failed:'] + ['    ' + l for l in installed.stderr.strip().splitlines()[-15:]]

    listed = {p['name'].lower(): p['version']
              for p in json.loads(sh([uv, 'pip', 'list', '--python', py, '--format', 'json']).stdout)}
    lines.append('  ' + ', '.join(f'{p} {listed.get(p, "-")}' for p in SHOWN))

    commands = {
        'jcast --help': [venv / 'bin' / 'jcast', '--help'],
        'python -m jcast --help': [py, '-m', 'jcast', '--help'],
        'version': [py, '-c', 'import jcast, importlib.metadata as m; '
                              'print(jcast.__version__, m.version("jcast"), jcast.__file__)'],
    }
    for label, cmd in commands.items():
        done = sh(cmd, cwd=work)
        shown = done.stdout.strip().replace(f'{venv}/', '<venv>/') if label == 'version' else f'exit {done.returncode}'
        lines.append(f'  {label}: {shown}' + (f' ({done.stderr.strip().splitlines()[-1]})' if done.returncode else ''))

    tests = work / f'tests-{name}'
    shutil.rmtree(tests, ignore_errors=True)
    with tarfile.open(sdist) as t:
        t.extractall(tests, filter='data')
    root = next(tests.iterdir())
    shutil.rmtree(root / 'jcast')
    started = time.monotonic()
    done = sh([py, '-m', 'pytest', '-q', '-p', 'no:cacheprovider'] + shlex.split(pytest_args), cwd=root,
              env={'PATH': f'{venv / "bin"}:/usr/bin:/bin', 'MPLBACKEND': 'Agg', 'HOME': str(work / 'home')})
    summary = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr.strip()
    lines.append(f'  tests: {re.sub(r" in [0-9.]+s.*", "", summary)} ({time.monotonic() - started:.0f} s)')
    if done.returncode:
        lines += ['    ' + l for l in done.stdout.splitlines() if l.startswith(('FAILED', 'ERROR'))][:30]
    return name, lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--uv', required=True)
    parser.add_argument('--python', action='append', required=True, help='an interpreter; repeat for each version')
    parser.add_argument('--resolution', choices=('newest', 'lowest'), default='newest')
    parser.add_argument('--ref', help='build from this commit instead of the working tree')
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--pytest_args', default="-m 'not network'")
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_check_install')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'check_install')
    parser.add_argument('--label', help='the evidence file is check_<label>.txt [default: the resolution]')
    args = parser.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / 'home').mkdir(exist_ok=True)
    sdist, wheel, tracked = build(args.uv, args.python[0], args.ref, args.work)
    head = sh(['git', 'rev-parse', '--short', args.ref or 'HEAD'], cwd=REPO).stdout.strip()
    changed = sh(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=REPO).stdout.strip()
    source = f'commit {head}' if args.ref or not changed else f'the working tree ({head} and uncommitted changes)'
    lines = [f'Built from {source} with '
             f'{sh([args.uv, "--version"]).stdout.strip()}; resolution: {args.resolution}; '
             f'pytest {args.pytest_args}', ''] + check_packages(sdist, wheel, tracked) + ['']

    with ThreadPoolExecutor(args.jobs) as pool:
        results = pool.map(lambda p: environment(args.uv, p, args.resolution, sdist, wheel, args.work,
                                                 args.pytest_args), args.python)
        for _, env_lines in results:
            lines += env_lines + ['']

    text = '\n'.join(lines).rstrip() + '\n'
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f'check_{args.label or args.resolution}.txt').write_text(text)


if __name__ == '__main__':
    main()
