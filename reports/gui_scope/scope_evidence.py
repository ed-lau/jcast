# -*- coding: utf-8 -*-

"""Evidence for the GUI and command-line scope of v0.5.0 (reports/2026-11-13-gui-cli-scope.md).

    typer       For each Typer version listed, its dependencies and whether it vendors Click (from the wheel's
                metadata, downloaded from PyPI). In a virtual environment per version in --parse, how typer_cli.py
                parses the command lines JCAST users write today, and whether a hand-written [default: ...] survives
                in the help.
    preview     The -m model fitted from an rMATS folder alone, as a GUI preview would before a run: the time to read
                the tables, compute the read counts and fit, and the cutoff. On the bundled chr15 folder, and on each
                --rmats folder.

Needs network access (PyPI). Writes only to --work and --out.

Usage (from the repository root, in an environment where jcast is importable):
    python reports/gui_scope/scope_evidence.py --work <dir> [--python <interpreter for the venvs>] [--rmats <dir> ...]
        [--out <dir>]
"""

import argparse
import logging
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

METADATA_VERSIONS = ['0.24.0', '0.25.1', '0.26.0', '0.28.0', '0.29.0', '0.29.1']
PARSE_VERSIONS = ['0.26.0', '0.29.1']
COMMAND_LINES = [['run', 'a', 'b', 'c', '-s', 'SE'],
                 ['run', 'a', 'b', 'c', '-s', 'SE', 'MXE'],
                 ['run', 'a', 'b', 'c', '-s', 'SE', '-s', 'MXE'],
                 ['run', 'a', 'b', 'c', '-q', '0', '0.05'],
                 ['run', 'a', 'b', 'c', '--g_or_ln', '1'],
                 ['a', 'b', 'c'],
                 ['gui']]
HELP_PHRASE = 'which keeps every junction'


def wheel_metadata(python, version, work) -> str:
    """ the wheel's Requires-Dist without extras, and whether it holds typer/_click """
    target = work / 'wheels' / version
    subprocess.run([python, '-m', 'pip', 'download', '-q', '--no-deps', '--only-binary', ':all:', '-d', str(target),
                    f'typer=={version}'], check=True, capture_output=True)
    with zipfile.ZipFile(next(target.glob('*.whl'))) as whl:
        meta = whl.read(next(n for n in whl.namelist() if n.endswith('.dist-info/METADATA'))).decode()
        vendored = any(n.startswith('typer/_click/') for n in whl.namelist())
    requires = [line.split(': ', 1)[1] for line in meta.splitlines()
                if line.startswith('Requires-Dist: ') and 'extra ==' not in line]
    return f'typer {version}: vendors Click: {"yes" if vendored else "no"}; requires: {", ".join(requires)}'


def venv_with(python, version, work) -> Path:
    """ a virtual environment holding typer==version and its dependencies """
    venv = work / f'venv-typer-{version}'
    if not (venv / 'bin' / 'python').exists():
        subprocess.run([python, '-m', 'venv', str(venv)], check=True)
        subprocess.run([str(venv / 'bin' / 'python'), '-m', 'pip', 'install', '-q', '--disable-pip-version-check',
                        f'typer=={version}'], check=True)
    return venv / 'bin' / 'python'


def parse_cases(venv_python) -> list:
    lines = []
    freeze = subprocess.run([str(venv_python), '-m', 'pip', 'list', '--format=freeze', '--disable-pip-version-check'],
                            capture_output=True, text=True, check=True).stdout.split()
    lines.append('  installed: ' + ', '.join(p for p in freeze if not p.startswith('pip==')))
    env = dict(os.environ, COLUMNS='100', NO_COLOR='1')
    for argv in COMMAND_LINES:
        done = subprocess.run([str(venv_python), str(HERE / 'typer_cli.py')] + argv, capture_output=True, text=True,
                              env=env)
        out = [line.strip(' │╭╮╰╯─') for line in (done.stdout + done.stderr).splitlines()]
        message = next((line for line in out if line.startswith('parsed:')), None) \
            or next((line for line in out if line.strip() not in ('', 'Error') and not line.startswith(('Usage:', 'Try '))),
                    '')
        lines.append(f'  {" ".join(argv):<28} exit {done.returncode}  {message.strip()}')
    help_text = subprocess.run([str(venv_python), str(HERE / 'typer_cli.py'), 'run', '-h'], capture_output=True,
                               text=True, env=env).stdout
    read_line = next(line.strip(' │') for line in help_text.splitlines() if '--read' in line)
    lines.append(f'  run -h, the -r line: {" ".join(read_line.split())}')
    lines.append(f'  hand-written "[default: 1, {HELP_PHRASE}]" in the help: '
                 f'{"kept" if HELP_PHRASE in help_text else "dropped"}')
    return lines


def preview(folder) -> str:
    """ the -m model from the rMATS folder alone, timed """
    from jcast import model
    from jcast.junctions import RmatsResults

    started = time.monotonic()
    results = RmatsResults(logging.getLogger('scope_evidence'), str(folder),
                           splice_types={'MXE', 'SE', 'RI', 'A5SS', 'A3SS'}, others_if_present=True)
    read = time.monotonic()
    results.get_junction_count_array()
    counted = time.monotonic()
    _, params, cutoff = model.general_mixture_model(min_jc_array=results.min_jc_array)
    fitted = time.monotonic()
    name = folder.relative_to(REPO) if folder.is_relative_to(REPO) else folder
    return (f'  {name}: {len(results.min_jc_array):,} junctions, cutoff {cutoff}; '
            f'read {read - started:.1f} s, counts {counted - read:.1f} s, fit {fitted - counted:.1f} s, '
            f'{fitted - started:.1f} s in all')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--work', required=True, type=Path, help='scratch directory for wheels and venvs')
    parser.add_argument('--python', default=sys.executable, help='interpreter to make the venvs with')
    parser.add_argument('--rmats', nargs='*', default=[], type=Path, help='more rMATS folders to time the preview on')
    parser.add_argument('--out', type=Path, default=HERE / 'scope_evidence')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    lines = ['Typer releases: dependencies, and whether Click is vendored (wheel metadata from PyPI)']
    lines += ['  ' + wheel_metadata(args.python, v, args.work) for v in METADATA_VERSIONS]
    for version in PARSE_VERSIONS:
        lines += ['', f'typer_cli.py on typer {version}']
        lines += parse_cases(venv_with(args.python, version, args.work))

    lines += ['', 'The -m model from the rMATS folder alone (no GTF or genome), as a GUI preview']
    for folder in [REPO / 'tests' / 'data' / 'rmats'] + [f.resolve() for f in args.rmats]:
        lines.append(preview(folder))

    text = '\n'.join(lines) + '\n'
    (args.out / 'scope_evidence.txt').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
