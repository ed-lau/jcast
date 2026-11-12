# -*- coding: utf-8 -*-

"""Verification of PR 16 (reports/2026-11-12-log-audit.md §7): junctions.tsv, the small log, --verbose.

For each case, runs `python -m jcast <rmats> <gtf> <genome> -c --uniprot <fasta> -o <dir> <args>` three times, with
PYTHONPATH set to a checkout (so the code under test can be frozen in a git worktree): the code before; the code
after; and the code after with --verbose. Then checks:

    output      every output file but the log, model.pdf and junctions.tsv is byte-identical before and after (with
                and without --verbose); junctions.tsv is identical with and without --verbose
    fates       junctions.tsv has one row per junction, in the order of the log before, with the same fate
    records     its new_records add up to the records in the psq_*.fasta files; its canonical sources are those of
                uniprot_canonicals.tsv
    verbose     the --verbose log has the messages of the log before, in the same order, once the line suffix, the
                header and the summary are set aside, and the DEBUG message 'Chosen coding region …', which the
                code before never logged (its logger was at INFO level)
    sizes       the log before, the log after, the --verbose log and junctions.tsv
    time        the junction loop: before, from the log's time stamps; after, from the summary

Cases: chr15 (the test set) with -r 1, -m, and -r 1 -p 3; with --rmats/--gtf/--genome/--uniprot, also a whole genome
with -m -p 16. Writes <out>/verify_16.txt.

Usage (from the repository root):
    python reports/log_audit/verify_16.py --before <checkout> --after <checkout> \\
        [--rmats <dir> --gtf <gtf> --genome <fa> --uniprot <fasta>] [--work <dir>] [--out <dir>]
"""

import argparse
import collections
import csv
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
CHR15 = dict(rmats=DATA / 'rmats', gtf=DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
             genome=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz',
             uniprot=DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz')

NOT_COMPARED = ('logfile.log', 'model.pdf', 'junctions.tsv')
STAMP = re.compile(r'^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3})~(\w+)~(.*)$')
SUFFIX = re.compile(r'~module:\w+~function:\w+$')
DOING = re.compile(r'^>>>>>> Doing (\w+) junction (\S+) for gene')


def run(paths, repo, extra, out, cache_dir):
    cmd = [sys.executable, '-m', 'jcast', str(paths['rmats']), str(paths['gtf']), str(paths['genome']), '-c',
           '--uniprot', str(paths['uniprot']), '-o', str(out)] + (['--cache_dir', str(cache_dir)] if cache_dir else [])
    subprocess.run(cmd + extra, cwd=out.parent, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg'))
    return out


def records(path):
    """ log records: (time stamp, level, message), a message's continuation lines joined to it """
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            m = STAMP.match(line)
            if m:
                out.append([m.group(1), m.group(2), m.group(3)])
            elif out:
                out[-1][2] += '\n' + line
    return [(t, level, SUFFIX.sub('', msg)) for t, level, msg in out]


def fates_before(log, fates):
    """ (type, id, fate name) per junction, from a log before PR 16 """
    heads = [m.split('{')[0].rstrip(' \n') for m, _ in fates.NAMES]
    rows, current = [], None
    for _, _, msg in log:
        m = DOING.match(msg)
        if m:
            current = m.groups()
        elif current and any(msg.startswith(h) for h in heads):
            # the fate messages end with ' \n\n'; the suffix removed, the message keeps ' \n'
            rows.append((*current, fates.name(msg.rstrip('\n') + '\n\n')[0]))
            current = None
    return rows


# logged at DEBUG level before PR 16 too, so never written: --verbose writes it
NEWLY_LOGGED = 'Chosen coding region is '


def junction_messages(log):
    """ the messages from the first junction to the last, without trailing blanks, levels aside """
    msgs = [msg.rstrip() for _, _, msg in log]
    first = next(i for i, m in enumerate(msgs) if m.startswith('>>>>>> Doing'))
    last = max(i for i, m in enumerate(msgs) if not m.startswith(('Canonicals by source', 'Summary.', 'Time: ')))
    return [m for m in msgs[first:last + 1] if m and not m.startswith(NEWLY_LOGGED)]


def loop_seconds(log):
    """ before: from the first junction's first record to its last junction's last record """
    t = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S,%f') for s, _, m in log]
    first = next(i for i, (_, _, m) in enumerate(log) if m.startswith('>>>>>> Doing'))
    end = next(i for i, (_, _, m) in enumerate(log) if m.startswith('Canonicals by source'))
    return (t[end - 1] - t[first]).total_seconds()


def compare(label, before, after, verbose, fates, lines):
    lines.append(f'{label}:')
    files = sorted(p.name for p in before.iterdir() if p.name not in NOT_COMPARED)
    for other in (after, verbose):
        assert sorted(p.name for p in other.iterdir() if p.name not in NOT_COMPARED) == files, other
    differ = [f for f in files if not ((before / f).read_bytes() == (after / f).read_bytes()
                                       == (verbose / f).read_bytes())]
    lines.append(f'  output files: {len(files)}; identical before, after and after --verbose: '
                 f'{"yes" if not differ else "no, " + str(differ)}')
    same_table = (after / 'junctions.tsv').read_bytes() == (verbose / 'junctions.tsv').read_bytes()
    lines.append(f'  junctions.tsv identical with and without --verbose: {"yes" if same_table else "NO"}')

    log_before = records(before / 'logfile.log')
    with open(after / 'junctions.tsv') as fh:
        table = list(csv.DictReader(fh, delimiter='\t'))
    old = fates_before(log_before, fates)
    new = [(r['junction_type'], r['junction_id'], r['fate']) for r in table]
    lines.append(f'  junctions: {len(new):,} rows; the log before names {len(old):,}; same junctions, order and fates: '
                 f'{"yes" if old == new else "NO"}')
    if old != new:
        diff = [(o, n) for o, n in zip(old, new) if o != n][:5]
        lines.append(f'    first differences (before, after): {diff}')
    counts = collections.Counter(r['fate'] for r in table)
    lines.append('  fates: ' + ', '.join(f'{name} {counts[name]:,}' for _, name in fates.NAMES if counts[name]))

    n_records = sum(sum(1 for line in open(p) if line.startswith('>')) for p in after.glob('psq_*.fasta'))
    n_new = sum(int(r['new_records']) for r in table)
    with open(after / 'uniprot_canonicals.tsv') as fh:
        chosen = {r['gene_id']: (r['source'], r['accession']) for r in csv.DictReader(fh, delimiter='\t')}
    agree = all(chosen[r['gene_id']] == (r['canonical_source'], r['canonical_accession']) for r in table
                if r['canonical_source'])
    lines.append(f'  new_records add up to the records in psq_*.fasta: {n_new:,} vs {n_records:,} '
                 f'({"yes" if n_new == n_records else "NO"}); canonical sources as in uniprot_canonicals.tsv: '
                 f'{"yes" if agree else "NO"}')

    log_verbose = records(verbose / 'logfile.log')
    same = junction_messages(log_before) == junction_messages(log_verbose)
    newly = sum(1 for _, _, m in log_verbose if m.startswith(NEWLY_LOGGED))
    lines.append(f'  --verbose log: the messages of the log before, from the first junction to the last: '
                 f'{"yes" if same else "NO"}; and {newly:,} "{NEWLY_LOGGED}…" messages')
    if not same:
        a, b = junction_messages(log_before), junction_messages(log_verbose)
        i = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        lines.append(f'    first difference at {i}: {a[i:i + 1]!r} vs {b[i:i + 1]!r}')
    levels = collections.Counter(level for _, level, _ in records(after / 'logfile.log'))
    lines.append(f'  log after, records by level: {dict(sorted(levels.items()))}')

    sizes = {name: (d / f).stat().st_size for name, d, f in [('log before', before, 'logfile.log'),
                                                             ('log after', after, 'logfile.log'),
                                                             ('log after --verbose', verbose, 'logfile.log'),
                                                             ('junctions.tsv', after, 'junctions.tsv')]}
    lines.append('  sizes: ' + '; '.join(f'{k} {v:,} bytes ({v / max(len(new), 1):,.0f} per junction)'
                                         for k, v in sizes.items()))
    time_line = [m for _, _, m in records(after / 'logfile.log') if m.startswith('Time: ')]
    lines.append(f'  junction loop: before {loop_seconds(log_before):.0f} s (time stamps); after: {time_line[0]}')
    lines.append('')


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + working tree' if d.stdout.strip() else '')


def main_(a):
    sys.path.insert(0, str(Path(a.after).resolve()))
    from jcast import fates

    work = Path(tempfile.mkdtemp(prefix='verify_16_', dir=a.work))
    lines = [f'before: {commit(a.before)}; after: {commit(a.after)}', '']
    cases = [('chr15 -r 1', CHR15, ['-r', '1'], work / 'cache'),
             ('chr15 -m', CHR15, ['-m'], work / 'cache'),
             ('chr15 -r 1 -p 3', CHR15, ['-r', '1', '-p', '3'], work / 'cache')]
    if a.rmats:
        genome = {k: Path(getattr(a, k)).resolve() for k in ('rmats', 'gtf', 'genome', 'uniprot')}
        cases.append(('genome -m -p 16', genome, ['-m', '-p', '16'], a.cache_dir))
    for label, paths, extra, cache_dir in cases:
        name = label.replace(' ', '_')
        before = run(paths, a.before, extra, work / f'{name}_before', cache_dir)
        after = run(paths, a.after, extra, work / f'{name}_after', cache_dir)
        verbose = run(paths, a.after, extra + ['--verbose'], work / f'{name}_verbose', cache_dir)
        compare(label, before, after, verbose, fates, lines)
        print('\n'.join(lines), flush=True)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'verify_16.txt').write_text('\n'.join(lines) + '\n')
    # examples: the new log and table of the chr15 -m run
    log = (work / 'chr15_-m_after' / 'logfile.log').read_text()
    for path, name in [(work, '<work>'), (Path(a.after).resolve(), '<after>'), (REPO, '<repo>')]:
        log = log.replace(str(path), name)
    (out / 'logfile_chr15_-m.txt').write_text(log)
    shutil.copy(work / 'chr15_-m_after' / 'junctions.tsv', out / 'junctions_chr15_-m.tsv')
    if a.keep:
        print(f'runs kept in {work}')
    else:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--before', required=True)
    p.add_argument('--after', required=True)
    p.add_argument('--rmats')
    p.add_argument('--gtf')
    p.add_argument('--genome')
    p.add_argument('--uniprot')
    p.add_argument('--cache_dir', default=None, help='for the whole genome [default: the user cache directory]')
    p.add_argument('--work', default=None)
    p.add_argument('--keep', action='store_true')
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'verify_16'))
    main_(p.parse_args())
