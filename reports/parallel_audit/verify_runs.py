# -*- coding: utf-8 -*-

"""Verification of PRs 9 and 15 (reports/2026-11-12-parallel-audit.md §8-9): whole-genome runs through the real CLI.

Writes a seeded random subset of an rMATS run (the same rows of each of the five files), then runs
`python -m jcast <subset> <gtf> <genome> -c -r 1 -q 0 1 --uniprot <fasta> [-p N]` once per run spec, each from its
own repository checkout, and compares every run with the first:

    output files     every file in the output directory except the log, byte for byte
    log              every log message in order, without its timestamp; messages of the GTF and genome loads differ
                     between versions only in wording and are listed separately
    time             wall time of the whole run, and of the junction loop (first to last junction in the log)

Writes <out_name>/verify.txt.

Usage (from the repository root):
    python reports/parallel_audit/verify_runs.py --rmats <dir> --gtf <gtf> --genome <fa> --uniprot <fasta> \\
        --sample 2000 --run before=<checkout of the earlier commit> --run after=. [--run after_p16=.:16] --out_name <name>
"""

import argparse
import datetime
import hashlib
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']


def subset(rmats, n, seed, out):
    """ the same seeded random sample of rows from each rMATS file, header kept, file order kept """
    out.mkdir()
    lines = {t: open(Path(rmats) / f'{t}.MATS.JC.txt').read().splitlines(keepends=True) for t in TYPES}
    rows = [(t, i) for t in TYPES for i in range(1, len(lines[t]))]
    keep = set(random.Random(seed).sample(rows, min(n, len(rows))))
    for t in TYPES:
        (out / f'{t}.MATS.JC.txt').write_text(lines[t][0] + ''.join(
            line for i, line in enumerate(lines[t]) if i and (t, i) in keep))
    return len(keep)


def log_messages(path):
    """ (timestamp, message) per log line; continuation lines belong to the line before """
    out = []
    for line in open(path):
        parts = line.rstrip('\n').split('~')
        try:
            stamp = datetime.datetime.strptime(parts[0], '%Y-%m-%d %H:%M:%S,%f')
            out.append([stamp, '~'.join(parts[1:])])
        except ValueError:
            if out:
                out[-1][1] += '\n' + line.rstrip('\n')
    return out


def main_(a):
    for k in ['rmats', 'gtf', 'genome', 'uniprot']:
        setattr(a, k, str(Path(getattr(a, k)).resolve()))
    res = HERE / a.out_name
    res.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='verify_runs_', dir=a.work))
    n = subset(a.rmats, a.sample, a.seed, work / 'rmats')
    out = [f'rMATS subset: {n} junctions (seed {a.seed}) of {a.rmats}\ngtf {a.gtf}\ngenome {a.genome}\n'
           f'uniprot {a.uniprot}\nmachine: {os.cpu_count()} CPUs; load average at start {os.getloadavg()[0]:.1f}\n']

    runs = []
    for spec in a.run:
        label, repo = spec.split('=', 1)
        repo, _, processes = repo.partition(':')
        repo = Path(repo).resolve()
        commit = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'],
                                capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'],
                               capture_output=True, text=True).stdout.strip()
        out_dir = work / label
        out_dir.mkdir()     # jcast does not create it (REVIEW.md MED-4)
        cmd = [sys.executable, '-m', 'jcast', str(work / 'rmats'), a.gtf, a.genome, '-c', '-r', '1', '-q', '0', '1',
               '--uniprot', a.uniprot, '-o', str(out_dir)] + (['-p', processes] if processes else [])
        t0 = time.perf_counter()
        subprocess.run(cmd, cwd=work, env=dict(os.environ, PYTHONPATH=str(repo)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        wall = time.perf_counter() - t0
        logs = log_messages(out_dir / 'logfile.log')
        doing = [s for s, m in logs if '>>>>>> Doing' in m]
        loop = (logs[-1][0] - doing[0]).total_seconds() if doing else float('nan')
        files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                 for p in sorted(out_dir.iterdir()) if p.name != 'logfile.log'}
        runs.append(dict(label=label, commit=commit + (' + working tree' if dirty else ''), processes=processes or '-',
                         wall=wall, loop=loop, files=files, messages=[m for _, m in logs]))

    first = runs[0]
    start = next(i for i, m in enumerate(first['messages']) if '>>>>>> Doing' in m)
    out.append('run\tcommit\tprocesses\twall_s\tjunction_loop_s\tjunctions_per_s\tfiles_identical\t'
               'messages_identical_from_first_junction')
    for r in runs:
        s = next(i for i, m in enumerate(r['messages']) if '>>>>>> Doing' in m)
        out.append(f"{r['label']}\t{r['commit']}\t{r['processes']}\t{r['wall']:.0f}\t{r['loop']:.0f}\t"
                   f"{n / r['loop']:.1f}\t{r['files'] == first['files']}\t"
                   f"{r['messages'][s:] == first['messages'][start:]}")
    out.append('')
    for r in runs[1:]:
        if r['files'] != first['files']:
            out.append(f"{r['label']}: files differing from {first['label']}: "
                       f"{sorted(k for k in set(r['files']) | set(first['files']) if r['files'].get(k) != first['files'].get(k))}")
        head_a, head_b = first['messages'][:start], r['messages'][:next(
            i for i, m in enumerate(r['messages']) if '>>>>>> Doing' in m)]
        if head_a != head_b:
            out.append(f"{r['label']}: start-up messages (before the first junction) differ:\n  - " +
                       '\n  - '.join(m.split('~module')[0][:150] for m in head_a if m not in head_b) + '\n  + ' +
                       '\n  + '.join(m.split('~module')[0][:150] for m in head_b if m not in head_a))
    out.append(f"output files: {', '.join(f'{k} {v}' for k, v in first['files'].items())}")
    (res / 'verify.txt').write_text('\n'.join(out) + '\n')
    print((res / 'verify.txt').read_text())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--rmats', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--uniprot', required=True)
    p.add_argument('--sample', type=int, default=2000)
    p.add_argument('--seed', type=int, default=3)
    p.add_argument('--run', action='append', required=True, help='label=repo[:processes]; the first is the reference')
    p.add_argument('--work', default=None)
    p.add_argument('--out_name', required=True)
    main_(p.parse_args())
