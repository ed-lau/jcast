# -*- coding: utf-8 -*-

"""Verification of PR 19 (reports/2026-11-13-gui-cli-scope.md §10): what the GUI needs from the pipeline, with no
output change.

With the code before and after (PYTHONPATH set to each checkout, so the code under test can be frozen in a git
worktree):

    runs        the chr15 test set, run as `python -m jcast <rmats> <GTF> <genome> -c --uniprot <Swiss-Prot> <args>`:
                every output file but the log and model.pdf is byte-identical, and so is the log but for timestamps,
                the output directory and the time taken.
    plot        model.pdf's figure, for the -m model of each bundled rMATS folder, saved as PNG (Figure.savefig
                patched): byte-identical.
    counts      the read count array of the -m model (RmatsResults.get_junction_count_array), on each bundled rMATS
                folder and each --rmats folder: identical values, dtype and shape; the time it takes, before and after
                back to back; and the cutoff, before, and from model.preview after.

Writes <out>/verify_19.txt.

Usage (from the repository root):
    python reports/gui_scope/verify_19.py --before <checkout> --after <checkout> [--rmats <dir> ...] [--work <dir>]
        [--out <dir>]
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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
RMATS = DATA / 'rmats'
INPUTS = [RMATS, DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
          DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']
RUNS = {'default': [], '-m': ['-m'], '-r 30 -q 0 0.05': ['-r', '30', '-q', '0', '0.05'], '-p 2': ['-p', '2']}
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
SET = re.compile(r"splice_type=\{([^}]*)\}")
FOLDERS = [RMATS] + sorted(p for p in (DATA / 'model').iterdir() if p.is_dir())

# run in each checkout: the plots as PNG, and the read count arrays, as JSON
PROBE = r'''
import hashlib, json, logging, sys, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.figure
from jcast import model
from jcast.junctions import RmatsResults

mode, png_dir, folders = sys.argv[1], sys.argv[2], sys.argv[3:]
saved = []
def savefig(self, fname, **kwargs):
    path = f'{png_dir}/{len(saved)}.png'
    self.canvas.print_figure(path, dpi=100, format='png', metadata={'Software': None})
    saved.append(path)
matplotlib.figure.Figure.savefig = savefig

import jcast, inspect
found = [{'jcast': jcast.__file__, 'counts_from': 'iloc per row' if 'iloc' in inspect.getsource(
    RmatsResults.get_junction_count_array) else 'columns'}]
for folder in folders:
    started = time.monotonic()
    results = RmatsResults(logging.getLogger('verify'), folder, splice_types=(), others_if_present=True)
    read = time.monotonic()
    results.get_junction_count_array()
    counted = time.monotonic()
    a = results.min_jc_array
    try:
        ln_sjc, params, cutoff = model.general_mixture_model(min_jc_array=a)
    except Exception as e:      # no junctions: the model fails, before and after
        params, cutoff = None, f'{type(e).__name__}'
    row = dict(folder=folder, junctions=len(a), dtype=str(a.dtype), shape=list(a.shape),
               sha256=hashlib.sha256(a.tobytes()).hexdigest(), counts_s=round(counted - read, 2), cutoff=cutoff)
    if mode == 'after':
        try:
            row['preview_cutoff'] = model.preview(folder)[3]
        except Exception as e:
            row['preview_cutoff'] = f'{type(e).__name__}'
    if params is not None and mode != 'counts_only':
        model.plot_general_mixture_model(ln_sjc, params, cutoff, write_dir=png_dir, filename='model')
        row['png'] = hashlib.sha256(open(saved[-1], 'rb').read()).hexdigest()
    found.append(row)
print(json.dumps(found))
'''


# Each process runs in a directory without a jcast package: python -m and -c put the current directory on sys.path
# before PYTHONPATH, so from the repository root both sides would import the working tree's jcast

def run(repo, extra, out, cache_dir):
    cmd = [sys.executable, '-m', 'jcast'] + INPUTS + ['-o', out, '--cache_dir', cache_dir] + extra
    out.mkdir(parents=True, exist_ok=True)
    done = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=out,
                          env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg'))
    return done.returncode, out


def probe(repo, mode, png_dir, folders):
    png_dir.mkdir(parents=True, exist_ok=True)
    done = subprocess.run([sys.executable, '-c', PROBE, mode, str(png_dir)] + [str(f) for f in folders],
                          capture_output=True, text=True, check=True, cwd=png_dir,
                          env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve())))
    return json.loads(done.stdout.splitlines()[-1])


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
            if p.is_file() and p.name not in NOT_COMPARED}


def log_lines(out):
    """ the log's lines, without timestamps or the time taken, the output directory written <out>, and -s's default
    set, whose order changes with each process, sorted """
    text = (out / 'logfile.log').read_text().replace(str(out), '<out>')
    text = SET.sub(lambda m: 'splice_type={' + ', '.join(sorted(m.group(1).split(', '))) + '}', text)
    return [STAMP.sub('', line) for line in text.splitlines() if not STAMP.sub('', line).startswith('INFO~Time: ')]


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    return c.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--rmats', nargs='*', default=[], type=Path, help='more rMATS folders, e.g. a whole genome')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_19')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_19')
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    cache_dir = args.work / 'cache'
    run(args.after, ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)

    jobs = [(label, side) for label in RUNS for side in ('before', 'after')]
    with ThreadPoolExecutor(4) as pool:
        results = dict(zip(jobs, pool.map(
            lambda job: run(getattr(args, job[1]), RUNS[job[0]], args.work / f'{job[1]}_{job[0]}'.replace(' ', '_'),
                            cache_dir), jobs)))

    lines = [f'before: {commit(args.before)}; after: {commit(args.after)}', '',
             '== chr15 runs (<inputs> -c --uniprot <chr15 Swiss-Prot> <args>): before and after']
    for label in RUNS:
        (rb, b), (ra, a) = results[(label, 'before')], results[(label, 'after')]
        fb, fa = files(b), files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lines.append(f'{label}: exit {rb} / {ra}; {len(fa)} output files, identical: '
                     f'{"yes" if not differ else "no, " + ", ".join(differ)}; log identical but for timestamps, the '
                     f'output directory, the time taken and the order of the -s set: {"yes" if log_lines(b) == log_lines(a) else "no"}')

    # the bundled folders: plots and counts; then each --rmats folder, counts only, before and after back to back
    bundled = {side: probe(getattr(args, side), side, args.work / f'png_{side}', FOLDERS) for side in ('before', 'after')}
    lines += ['', '== The -m model of each bundled rMATS folder: model.pdf drawn as PNG, and the read count array']
    lines += [f'{side}: jcast from {bundled[side][0]["jcast"]}, read counts from {bundled[side][0]["counts_from"]}'
              for side in ('before', 'after')]
    for b, a in zip(bundled['before'][1:], bundled['after'][1:]):
        name = Path(b['folder']).relative_to(REPO)
        same = all(b[k] == a[k] for k in ('junctions', 'dtype', 'shape', 'sha256', 'cutoff'))
        png = 'no fit' if 'png' not in b else ('identical' if b['png'] == a.get('png') else 'DIFFERENT')
        lines.append(f'{name}: {a["junctions"]:,} junctions; counts identical: {"yes" if same else "NO"}; cutoff '
                     f'{b["cutoff"]} / {a["cutoff"]}, preview {a["preview_cutoff"]}; plot: {png}')

    if args.rmats:
        lines += ['', '== More rMATS folders: the read count array, and the time to compute it (before, then after)']
    for folder in args.rmats:
        b = probe(args.before, 'counts_only', args.work / 'png_counts', [folder.resolve()])[1]
        a = probe(args.after, 'after', args.work / 'png_counts', [folder.resolve()])[1]
        same = all(b[k] == a[k] for k in ('junctions', 'dtype', 'shape', 'sha256', 'cutoff'))
        lines.append(f'{folder}: {a["junctions"]:,} junctions; counts identical: {"yes" if same else "NO"} '
                     f'({a["dtype"]}, shape {tuple(a["shape"])}); cutoff {b["cutoff"]} / {a["cutoff"]}, preview '
                     f'{a["preview_cutoff"]}; counts in {b["counts_s"]} s / {a["counts_s"]} s')

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'verify_19.txt').write_text(text)


if __name__ == '__main__':
    main()
