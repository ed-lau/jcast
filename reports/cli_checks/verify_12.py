# -*- coding: utf-8 -*-

"""Verification of PR 12 (reports/2026-11-13-command-line-checks.md): only the rMATS files of the -s splice types
must exist; with -m, the others are read too, when present, so the model's cutoff is unchanged.

Runs `python -m jcast <rmats> <chr15 GTF> <chr15 genome> -c --uniprot <chr15 Swiss-Prot> -o <dir> <args>` with the
code before and after (PYTHONPATH set to each checkout, so the code under test can be frozen in a git worktree):

    same        on the bundled rMATS folder (all five files), with and without -m and -s: every output file but the log
                and model.pdf is byte-identical, and so are the cutoff and the log's summary. The log lines that differ
                are listed (timestamps removed).
    se_only     on a folder holding only SE.MATS.JC.txt, with -s SE: the code before fails; after, the output equals
                that of the bundled folder with -s SE. With -m, the model fits on the SE junctions alone.

Writes <out>/verify_12.txt.

Usage (from the repository root):
    python reports/cli_checks/verify_12.py --before <checkout> --after <checkout> [--work <dir>] [--out <dir>]
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
INPUTS = [DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
          DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
CUTOFF = re.compile(r'Translating junctions with a read count of at least (\d+)')
FUNC = re.compile(r', func=<function run_jcast at 0x[0-9a-f]+>')
SET = re.compile(r"splice_type=\{([^}]*)\}")

SAME = {'default': [], '-m': ['-m'], '-s SE -m': ['-s', 'SE', '-m'], '-s SE RI': ['-s', 'SE', 'RI']}


def run(repo, rmats, extra, out, cache_dir):
    """ the run's exit status, its last line of stderr, and its output directory """
    cmd = [sys.executable, '-m', 'jcast', rmats] + INPUTS + ['-o', out, '--cache_dir', cache_dir] + extra
    done = subprocess.run([str(c) for c in cmd], cwd=out.parent, capture_output=True, text=True,
                          env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg'))
    error = done.stderr.strip().splitlines()[-1] if done.returncode and done.stderr.strip() else ''
    return done.returncode, error, out


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
            if p.name not in NOT_COMPARED}


def log_lines(out):
    """ the log's lines, without their timestamps, with the output directory written <out> """
    text = (out / 'logfile.log').read_text().replace(str(out), '<out>')
    return [STAMP.sub('', line) for line in text.splitlines()]


def arguments(line):
    """ the Arguments line without the function's address, with the -s set in a fixed order """
    line = FUNC.sub('', line)
    return SET.sub(lambda m: 'splice_type={' + ', '.join(sorted(m.group(1).split(', '))) + '}', line)


def summary(out):
    lines = log_lines(out)
    start = next(i for i, l in enumerate(lines) if l.endswith('Summary.'))
    end = next(i for i, l in enumerate(lines) if l.startswith('INFO~Time: '))
    return lines[start:end]


def cutoff(out):
    found = CUTOFF.search((out / 'logfile.log').read_text())
    return int(found.group(1)) if found else None


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + changes to jcast/' if d.stdout.strip() else '')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_12')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_12')
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    cache_dir = args.work / 'cache'
    se_only = args.work / 'se_only'
    se_only.mkdir(parents=True)
    shutil.copy(DATA / 'rmats' / 'SE.MATS.JC.txt', se_only)

    # warm the GTF cache once, so that the parallel runs do not all build it
    run(args.after, DATA / 'rmats', ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)

    jobs = {}
    for label, extra in SAME.items():
        for side in ('before', 'after'):
            jobs[(label, side)] = (getattr(args, side), DATA / 'rmats', extra)
    for side in ('before', 'after'):
        jobs[('SE only, -s SE', side)] = (getattr(args, side), se_only, ['-s', 'SE'])
    jobs[('SE only, -s SE -m', 'after')] = (args.after, se_only, ['-s', 'SE', '-m'])
    jobs[('-s SE', 'after')] = (args.after, DATA / 'rmats', ['-s', 'SE'])

    def one(key):
        repo, rmats, extra = jobs[key]
        out = args.work / f'{key[0]}-{key[1]}'.replace(' ', '_').replace('-', '')
        return key, run(repo, rmats, extra, out, cache_dir)

    with ThreadPoolExecutor(4) as pool:
        results = dict(pool.map(one, jobs))

    lines = [f'before: {commit(args.before)}; after: {commit(args.after)}',
             'runs: <rmats> <chr15 GTF> <chr15 genome> -c --uniprot <chr15 Swiss-Prot> <args>', '',
             '== The bundled rMATS folder (five files): before and after']
    for label in SAME:
        (rb, eb, b), (ra, ea, a) = results[(label, 'before')], results[(label, 'after')]
        if rb or ra:
            lines.append(f'{label}: exit {rb} / {ra} {eb} {ea}')
            continue
        fb, fa = files(b), files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lines.append(f'{label}: {len(fa)} output files, identical: {"yes" if not differ else "no, " + ", ".join(differ)}; '
                     f'cutoff {cutoff(b)} / {cutoff(a)}; summary identical: {"yes" if summary(b) == summary(a) else "no"}')
        lb, la = log_lines(b), log_lines(a)
        args_b = [arguments(l) for l in lb if l.startswith('INFO~Arguments: ')]
        args_a = [arguments(l) for l in la if l.startswith('INFO~Arguments: ')]
        lines.append(f'  log, Arguments: identical but for g_or_ln=0 before: '
                     f'{"yes" if [l.replace(", g_or_ln=0", "") for l in args_b] == args_a else "no"}')
        skip = ('INFO~Arguments: ', 'INFO~Time: ')
        lines += [f'  log, before only: {l}' for l in lb if l not in la and not l.startswith(skip)]
        lines += [f'  log, after only:  {l}' for l in la if l not in lb and not l.startswith(skip)]

    lines += ['', '== A folder holding only SE.MATS.JC.txt']
    rb, eb, _ = results[('SE only, -s SE', 'before')]
    lines.append(f'before, -s SE: exit {rb}: {eb}')
    ra, ea, a = results[('SE only, -s SE', 'after')]
    _, _, full = results[('-s SE', 'after')]
    same = files(a) == files(full) and summary(a) == summary(full)
    lines.append(f'after, -s SE: exit {ra}; output and summary identical to the bundled folder with -s SE: '
                 f'{"yes" if same else "no"}')
    ra, ea, a = results[('SE only, -s SE -m', 'after')]
    _, _, five = results[('-s SE -m', 'after')]
    lines.append(f'after, -s SE -m: exit {ra}; cutoff {cutoff(a)} (the model on the SE junctions alone; '
                 f'{cutoff(five)} with the five files present)')

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'verify_12.txt').write_text(text)


if __name__ == '__main__':
    main()
