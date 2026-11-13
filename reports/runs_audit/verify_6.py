# -*- coding: utf-8 -*-

"""Verification of PR 6 (reports/2026-11-12-thresholds-runs-canonical-audit.md §9): runs pairs of jcast CLI runs
and compares them.

PR 6 keeps a junction if its read count is at least the cutoff, where it used to need one read more, and returns
the model's cutoff in reads (one less than before). So the code after, with -r N, must give exactly the output of
the code before with -r N-1; and with -m (cutoff c in reads), that of the code before with -r c-1. PR 6 also applies
the filters first, which changes no output file.

Each run is `python -m jcast <rmats> <gtf> <genome> -c -q 0 1 --uniprot <fasta> -o <dir> --cache_dir <dir> <args>`
with PYTHONPATH set to a checkout (so the code under test can be frozen in a git worktree). For each pair, compares
every output file byte for byte (model.pdf and the log excepted), counts fates from each log, and times the junction
loop. Writes <out>/verify_6.txt.

Usage (from the repository root):
    python reports/runs_audit/verify_6.py --rmats <dir> --gtf <gtf> --genome <fa> --uniprot <fasta> \\
        --pair 'label|<before checkout>|<before args>|<after checkout>|<after args>' [...] --out <dir>
An after-args value of '-m' takes the before run's cutoff from the after run's log: write '-r {cutoff-1}' as the
before args.
"""

import argparse
import collections
import datetime
import hashlib
import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

FATE = re.compile(r'~INFO~(SUCCESS TIER 1 ALT TERMINUS|SUCCESS TIER [12]|PARTIAL TIER 4|FAILURE \(exception\)|FAILURE|'
                  r'SKIPPED\. [^.~]*)')
CUTOFF = re.compile(r'Translating junctions with a read count of at least (\d+)')


def run(a, repo, extra, out, cache_dir):
    cmd = [sys.executable, '-m', 'jcast', a.rmats, a.gtf, a.genome, '-c', '-q', '0', '1', '--uniprot', a.uniprot,
           '-o', str(out), '--cache_dir', str(cache_dir)] + shlex.split(extra)
    subprocess.run(cmd, cwd=out.parent, env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg'),
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    text = (out / 'logfile.log').read_text()
    stamps = [datetime.datetime.strptime(line[:23], '%Y-%m-%d %H:%M:%S,%f') for line in text.splitlines()
              if '>>>>>> Doing' in line]
    end = datetime.datetime.strptime([l for l in text.splitlines() if l[:4].isdigit()][-1][:23], '%Y-%m-%d %H:%M:%S,%f')
    cutoff = CUTOFF.search(text)
    return dict(fates=collections.Counter(FATE.findall(text)),
                files={p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in sorted(out.iterdir())
                       if p.name not in ('logfile.log', 'model.pdf')},
                loop=(end - stamps[0]).total_seconds() if stamps else float('nan'),
                cutoff=int(cutoff.group(1)) if cutoff else None)


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + working tree' if d.stdout.strip() else '')


def main_(a):
    for k in ['rmats', 'gtf', 'genome', 'uniprot']:
        setattr(a, k, str(Path(getattr(a, k)).resolve()))
    work = Path(tempfile.mkdtemp(prefix='verify_6_', dir=a.work))
    cache_dir = work / 'cache'
    lines = [f'rmats {a.rmats}\ngtf {a.gtf}\ngenome {a.genome}\nuniprot {a.uniprot}', '']
    for spec in a.pair:
        label, before_repo, before_args, after_repo, after_args = spec.split('|')
        after = run(a, after_repo, after_args, work / f'{label}_after', cache_dir)
        before_args = before_args.replace('{cutoff-1}', str((after['cutoff'] or 0) - 1))
        before = run(a, before_repo, before_args, work / f'{label}_before', cache_dir)
        differ = sorted(k for k in set(before['files']) | set(after['files'])
                        if before['files'].get(k) != after['files'].get(k))
        lines.append(f'{label}: before {commit(before_repo)} `{before_args}`; after {commit(after_repo)} `{after_args}`'
                     f' (cutoff logged: {after["cutoff"]})')
        lines.append(f'  output files: {len(after["files"])}; identical: {not differ}' +
                     (f'; differing: {differ}' if differ else ''))
        for name in sorted(set(before['fates']) | set(after['fates'])):
            lines.append(f'  {name}: before {before["fates"][name]}, after {after["fates"][name]}')
        lines.append(f'  junction loop: before {before["loop"]:.0f} s, after {after["loop"]:.0f} s')
        lines.append('')
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'verify_6.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--rmats', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--uniprot', required=True)
    p.add_argument('--pair', action='append', required=True)
    p.add_argument('--work', default=None)
    p.add_argument('--out', required=True)
    main_(p.parse_args())
