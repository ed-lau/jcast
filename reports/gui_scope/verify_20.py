# -*- coding: utf-8 -*-

"""Verification of PR 20 (reports/2026-11-13-gui-cli-scope.md §11): when the -m model has no fit, for any reason, the
run warns and falls back to -r.

Runs `python -m jcast <rmats> <chr15 GTF> <chr15 genome> -c --uniprot <chr15 Swiss-Prot> -o <dir> <args>` with the
code before and after, and optionally with an earlier checkout (the code before PR 19), each from a scratch directory
with PYTHONPATH set to its checkout:

    chr15       the bundled rMATS folder, default and -m: output files byte-identical, logs identical but for
                timestamps, the output directory and the time taken.
    degenerate  encode_human_stomach (no fit) with -m -s RI -r 100000: the same output; the log names -r as the
                cutoff's option.
    empty       encode_mouse_fetalheart (five header-only rMATS files), with and without -m: the exit status, the
                last error or warning, and the cutoff logged.

Writes <out>/verify_20.txt.

Usage (from the repository root):
    python reports/gui_scope/verify_20.py --before <checkout> --after <checkout> [--earlier <checkout>]
        [--work <dir>] [--out <dir>]
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
INPUTS = [DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
          DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']
CASES = {'chr15, default': (DATA / 'rmats', []),
         'chr15, -m': (DATA / 'rmats', ['-m']),
         'encode_human_stomach, -m -s RI -r 100000': (DATA / 'model' / 'encode_human_stomach',
                                                      ['-m', '-s', 'RI', '-r', '100000']),
         'header-only files, -m': (DATA / 'model' / 'encode_mouse_fetalheart', ['-m']),
         'header-only files, -m -r 3': (DATA / 'model' / 'encode_mouse_fetalheart', ['-m', '-r', '3']),
         'header-only files, without -m': (DATA / 'model' / 'encode_mouse_fetalheart', [])}
COMPARED = ['chr15, default', 'chr15, -m', 'encode_human_stomach, -m -s RI -r 100000']
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
SET = re.compile(r"splice_type=\{([^}]*)\}")
CUTOFF = re.compile(r'Translating junctions with a read count of at least \d+ \(-[mr]\)')


def env(repo):
    return dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg', COLUMNS='1000')


def run(repo, rmats, extra, out, cache_dir):
    """ the exit status, the last error or warning, and the output directory. Run from out: python -m puts the
    current directory on sys.path before PYTHONPATH """
    out.mkdir(parents=True)
    cmd = [sys.executable, '-m', 'jcast', rmats] + INPUTS + ['-o', out, '--cache_dir', cache_dir] + extra
    done = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=out, env=env(repo))
    lines = [l for l in done.stderr.splitlines() if re.search(r'Error|WARNING', l)]
    return done.returncode, lines[-1].strip() if lines else '', out


def imported(repo, work):
    return subprocess.run([sys.executable, '-c', 'import jcast; print(jcast.__file__)'], capture_output=True,
                          text=True, cwd=work, env=env(repo)).stdout.strip()


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
            if p.is_file() and p.name not in NOT_COMPARED}


def log_lines(out):
    text = (out / 'logfile.log').read_text().replace(str(out), '<out>')
    text = SET.sub(lambda m: 'splice_type={' + ', '.join(sorted(m.group(1).split(', '))) + '}', text)
    return [STAMP.sub('', line) for line in text.splitlines() if not STAMP.sub('', line).startswith('INFO~Time: ')]


def cutoff(out):
    found = CUTOFF.search((out / 'logfile.log').read_text()) if (out / 'logfile.log').exists() else None
    return found.group(0).replace('Translating junctions with a read count of ', '') if found else '-'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--earlier', help='a checkout of the code before PR 19, for the header-only files')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_20')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_20')
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    cache_dir = args.work / 'cache'
    sides = ['before', 'after'] + (['earlier'] if args.earlier else [])

    # warm the GTF cache, so that no run's log says it wrote it
    run(args.after, DATA / 'rmats', ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)

    lines = [f'{side}: {imported(getattr(args, side), args.work).replace(str(args.work), "<work>")}' for side in sides]
    results = {}
    for label, (rmats, extra) in CASES.items():
        for side in sides:
            if side == 'earlier' and 'header-only' not in label:
                continue
            out = args.work / re.sub(r'[^A-Za-z0-9]+', '_', f'{side} {label}')
            results[(label, side)] = run(getattr(args, side), rmats, extra, out, cache_dir)

    lines += ['', '== Runs that worked before: before and after']
    for label in COMPARED:
        (rb, _, b), (ra, _, a) = results[(label, 'before')], results[(label, 'after')]
        fb, fa = files(b), files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lb, la = log_lines(b), log_lines(a)
        lines.append(f'{label}: exit {rb} / {ra}; {len(fa)} output files, identical: '
                     f'{"yes" if not differ else "no, " + ", ".join(differ)}; cutoff {cutoff(b)} / {cutoff(a)}')
        lines += [f'  log, before only: {l}' for l in lb if l not in la]
        lines += [f'  log, after only:  {l}' for l in la if l not in lb]

    lines += ['', '== Header-only rMATS files (no junctions): exit status, last error or warning, cutoff']
    for label in [l for l in CASES if 'header-only' in l]:
        lines.append(label)
        for side in sides:
            r, message, out = results[(label, side)]
            tsv = out / 'junctions.tsv'
            rows = len(tsv.read_text().splitlines()) - 1 if tsv.exists() else '-'
            lines.append(f'  {side}: exit {r}; {message or "(nothing)"}; cutoff {cutoff(out)}; junctions.tsv rows '
                         f'{rows}')

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'verify_20.txt').write_text(text)


if __name__ == '__main__':
    main()
