# -*- coding: utf-8 -*-

"""Verification of PR 24 (reports/2026-11-15-read-count-model-audit.md §17): JCAST reads rMATS folders made without
the statistics (--statoff) or from one sample.

With the code before and after, each from a scratch directory with PYTHONPATH set to its checkout:

    chr15       the bundled rMATS folder, default and -m, run with
                `python -m jcast <rmats> <chr15 GTF> <chr15 genome> -c --uniprot <chr15 Swiss-Prot> -o <dir> <args>`:
                output files byte-identical, logs identical but for timestamps, the output directory and the time taken.
    bundled     every bundled rMATS folder (tests/data/rmats and tests/data/model), read as the GUI's Model page and
                -m read them (model.preview): the events, the junctions -r 1 -q 0 1 keeps, and the -m cutoff.
    folders     (--folder, whole-genome rMATS folders) the same, and the events left out (after only).

And (--union) one run of the code after on a --statoff folder, against two earlier runs at the same settings: one on
the folder made with statistics, one on the events the statistics left out (dropped_events.py --dropped). Its
junctions.tsv should be the union of theirs, but for the FDR (NA) and new_records (which junction writes a sequence
first depends on the order), and each output FASTA's sequences the union of theirs.

Writes <out>/verify_24.txt.

Usage (from the repository root):
    python reports/model_audit/verify_24.py --before <checkout> --after <checkout>
        [--folder <label>=<rMATS folder> ...] [--union <statoff run>,<statistics run>,<dropped run>]
        [--work <dir>] [--out <dir>]
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
INPUTS = [DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
          DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']
RUNS = {'chr15, default': [], 'chr15, -m': ['-m']}
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
SET = re.compile(r"splice_type=\{([^}]*)\}")

# model.preview and model.filter_counts as both checkouts have them (PR 19)
PREVIEW = """
import json, sys
from jcast import model
counts, _, _, cutoff, no_fit = model.preview(sys.argv[1])
kept = model.filter_counts(counts, 1, (0, 1))
from jcast.junctions import RmatsResults
import logging
left = getattr(RmatsResults(logging.getLogger('verify'), sys.argv[1]), 'left_out', None)
print(json.dumps({'events': len(counts), 'no FDR': int(counts.fdr.isna().sum()),
                  'read count 0': int((counts.min_jc == 0).sum()), 'kept by -r 1 -q 0 1': int(kept.kept.sum()),
                  '-m cutoff': None if cutoff is None else int(cutoff), 'no fit': no_fit,
                  'left out': None if left is None else sum(left.values())}))
"""


def env(repo):
    return dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg', COLUMNS='1000')


def run(repo, rmats, extra, out, cache_dir):
    """ the exit status and the output directory. Run from out: python -m puts the current directory on sys.path
    before PYTHONPATH """
    out.mkdir(parents=True)
    cmd = [sys.executable, '-m', 'jcast', rmats] + INPUTS + ['-o', out, '--cache_dir', cache_dir] + extra
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=out, env=env(repo)).returncode, out


def preview(repo, rmats, work) -> dict:
    done = subprocess.run([sys.executable, '-c', PREVIEW, str(rmats)], capture_output=True, text=True, cwd=work,
                          env=env(repo))
    if done.returncode:
        return {'error': done.stderr.strip().splitlines()[-1]}
    return json.loads(done.stdout.strip().splitlines()[-1])


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


def sequences(path) -> set:
    with open(path) as fh:
        return {str(r.seq) for r in SeqIO.parse(fh, 'fasta')}


def union(statoff, stat, dropped) -> list:
    """ the --statoff run against the union of the other two """
    lines = [f'statoff run: {statoff}', f'statistics run: {stat}', f'dropped events run: {dropped}']
    key = ['junction_type', 'junction_id']
    read = {name: pd.read_table(Path(d) / 'junctions.tsv', dtype=str, keep_default_na=False)
            for name, d in [('statoff', statoff), ('stat', stat), ('dropped', dropped)]}
    both = pd.concat([read['stat'], read['dropped']], ignore_index=True)
    lines.append(f"junctions: statoff run {len(read['statoff'])}; statistics run {len(read['stat'])} + dropped "
                 f"events run {len(read['dropped'])} = {len(both)}; the same junctions: "
                 f"{set(map(tuple, read['statoff'][key].values)) == set(map(tuple, both[key].values))}")
    lines.append(f"FDR NA in the statoff run: {(read['statoff'].fdr == 'NA').all()}")
    compared = [c for c in read['statoff'].columns if c not in ('fdr', 'new_records')]
    a = read['statoff'][compared].sort_values(key).reset_index(drop=True)
    b = both[compared].sort_values(key).reset_index(drop=True)
    differ = [c for c in compared if not a[c].equals(b[c])] if len(a) == len(b) else ['(lengths differ)']
    lines.append(f"every other column but new_records identical ({', '.join(compared[2:])}): "
                 f"{'yes' if not differ else 'no: ' + ', '.join(differ)}")
    lines.append('fates: ' + ', '.join(f'{f} {n}' for f, n in read['statoff'].fate.value_counts().items()))
    for name in sorted(p.name for p in Path(statoff).glob('*.fasta')):
        s = sequences(Path(statoff) / name)
        u = set().union(*(sequences(Path(d) / name) for d in (stat, dropped) if (Path(d) / name).exists()))
        lines.append(f'{name}: {len(s)} sequences; the union of the other two: {s == u}')
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--folder', action='append', default=[], help='<label>=<rMATS folder>')
    parser.add_argument('--union', help='<statoff run>,<statistics run>,<dropped events run>')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_24')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_24')
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    cache_dir = args.work / 'cache'

    # warm the GTF cache, so that no run's log says it wrote it
    run(args.after, DATA / 'rmats', ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)

    lines = [f'{side}: {imported(getattr(args, side), args.work).replace(str(args.work), "<work>")}'
             for side in ('before', 'after')]

    lines += ['', '== chr15 runs: before and after']
    for label, extra in RUNS.items():
        (rb, b), (ra, a) = (run(getattr(args, side), DATA / 'rmats', extra,
                                args.work / re.sub(r'[^A-Za-z0-9]+', '_', f'{side} {label}'), cache_dir)
                            for side in ('before', 'after'))
        fb, fa = files(b), files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lb, la = log_lines(b), log_lines(a)
        lines.append(f'{label}: exit {rb} / {ra}; {len(fa)} output files, identical: '
                     f'{"yes" if not differ else "no, " + ", ".join(differ)}')
        lines += [f'  log, before only: {l}' for l in lb if l not in la]
        lines += [f'  log, after only:  {l}' for l in la if l not in lb]

    folders = [('bundled ' + p.name, p) for p in [DATA / 'rmats'] + sorted(q for q in (DATA / 'model').iterdir()
                                                                          if q.is_dir())]
    folders += [(f.split('=', 1)[0], Path(f.split('=', 1)[1])) for f in args.folder]
    lines += ['', '== rMATS folders read (model.preview): before / after']
    for label, folder in folders:
        b, a = preview(args.before, folder, args.work), preview(args.after, folder, args.work)
        left = a.pop('left out', None)
        b.pop('left out', None)
        if a == b:
            found = (f"identical: events {a['events']}, kept {a['kept by -r 1 -q 0 1']}, -m cutoff {a['-m cutoff']}"
                     if 'error' not in a else f"both fail: {a['error']}")
        else:
            found = '; '.join(f'{k} {b.get(k)} / {v}' for k, v in a.items() if b.get(k) != v)
        lines.append(f'{label}: {found}; left out (after) {left}')

    if args.union:
        lines += ['', '== A run on a --statoff folder, against the runs on its two parts']
        lines += union(*args.union.split(','))

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'verify_24.txt').write_text(text)


if __name__ == '__main__':
    main()
