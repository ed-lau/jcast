# -*- coding: utf-8 -*-

"""Verification of PR 17 (reports/2026-11-12-ensembl-canonical-audit.md §10): the GTF canonical from the
Ensembl_canonical tag, without readthrough transcripts.

    chr15       runs `python -m jcast <chr15 test set> -c --uniprot <fasta> -o <dir> <args>` with the code before and
                after (PYTHONPATH set to each checkout, so the code under test can be frozen in a git worktree), with
                -r 1 and with -m. Every output file but the log and model.pdf is byte-identical: the bundled Ensembl
                89 GTF has no Ensembl_canonical tag.
    genome      from three whole-genome runs with -m -c --uniprot -p 16 (not bundled), given as output directories:
                the run before (--genome_before), the audit's prototype of the decided rule (--genome_prototype,
                `JCAST_AUDIT_EC_FALLBACK=no_readthrough` with prototype.diff) and the run after (--genome_after):
                - after and prototype: every output file but the log and model.pdf is byte-identical, and so is
                  the log's summary;
                - before and after: the files that differ and the genes whose canonical changed (as in
                  `audit_canonical.py compare`); every junctions.tsv column but the canonical and the records
                  written is unchanged; jcast_uniprot_merged.fasta and psq_isoforms.tsv are unchanged;
                - the time each run took translating junctions, from its log.

Writes <out>/verify_17.txt and, with the whole genome, <out>/compare_pr17.txt.

Usage (from the repository root):
    python reports/canonical_audit/verify_17.py --before <checkout> --after <checkout> \\
        [--genome_before <dir> --genome_prototype <dir> --genome_after <dir>] [--work <dir>] [--out <dir>]
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
CHR15 = [DATA / 'rmats', DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
         DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
         DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']

NOT_COMPARED = ('logfile.log', 'model.pdf')
UNCHANGED = ('jcast_uniprot_merged.fasta', 'psq_isoforms.tsv')
# the junctions.tsv columns a change of canonical can change
CANONICAL_COLUMNS = ('canonical_source', 'canonical_accession', 'files', 'new_records')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~\w+~')


def run(repo, extra, out, cache_dir):
    cmd = [sys.executable, '-m', 'jcast'] + [str(c) for c in CHR15] + ['-o', str(out), '--cache_dir', str(cache_dir)]
    subprocess.run(cmd + extra, cwd=out.parent, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg'))
    return out


def outputs(d):
    return sorted(p.name for p in d.iterdir() if p.name not in NOT_COMPARED)


def differing(a, b):
    assert outputs(a) == outputs(b), (a, b)
    return [f for f in outputs(a) if (a / f).read_bytes() != (b / f).read_bytes()]


def log_messages(d):
    """ the log's messages, a message's continuation lines joined to it """
    msgs = []
    for line in (d / 'logfile.log').read_text().split('\n'):
        if STAMP.match(line):
            msgs.append(STAMP.sub('', line))
        elif msgs:
            msgs[-1] += '\n' + line
    return msgs


def summary(d):
    return next(m for m in log_messages(d) if m.startswith('Summary.'))


def time_line(d):
    return next(m for m in log_messages(d) if m.startswith('Time: ')).rstrip()


def table(d):
    with open(d / 'junctions.tsv') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def genome(a, lines):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import audit_canonical

    before, prototype, after = (Path(x).resolve() for x in (a.genome_before, a.genome_prototype, a.genome_after))
    lines.append('genome -m -c --uniprot -p 16:')
    differ = differing(prototype, after)
    lines.append(f'  after and the prototype (no_readthrough): {len(outputs(after))} output files; identical: '
                 f'{"yes" if not differ else "NO, " + str(differ)}; log summary identical: '
                 f'{"yes" if summary(prototype) == summary(after) else "NO"}')

    rows_before, rows_after = table(before), table(after)
    others = [c for c in rows_before[0] if c not in CANONICAL_COLUMNS]
    same = len(rows_before) == len(rows_after) and all(
        [r[c] for c in others] == [s[c] for c in others] for r, s in zip(rows_before, rows_after))
    changed = sum(1 for r, s in zip(rows_before, rows_after) if r != s)
    lines.append(f'  before and after: junctions.tsv {len(rows_after):,} rows, {changed:,} differ; every column but '
                 f'{", ".join(CANONICAL_COLUMNS)} unchanged (fates, phases, coding regions): '
                 f'{"yes" if same else "NO"}')
    for f in UNCHANGED:
        lines.append(f'  before and after: {f} identical: '
                     f'{"yes" if (before / f).read_bytes() == (after / f).read_bytes() else "NO"}')
    for label, d in (('before', before), ('prototype', prototype), ('after', after)):
        lines.append(f'  {label}: {time_line(d)}')

    audit_canonical.compare(argparse.Namespace(before=before, after=after, label='pr17', out=a.out))
    compared = Path(a.out) / 'compare_pr17.txt'
    text = compared.read_text().replace(f'{before} -> {after}', '<PR 16 run> -> <PR 17 run>')
    compared.write_text(text)
    lines.append(f'  before and after, by gene: compare_pr17.txt')
    lines.append('')


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + working tree' if d.stdout.strip() else '')


def main_(a):
    work = Path(tempfile.mkdtemp(prefix='verify_17_', dir=a.work))
    lines = [f'before: {commit(a.before)}; after: {commit(a.after)}', '']
    for label, extra in (('chr15 -r 1', ['-r', '1']), ('chr15 -m', ['-m'])):
        name = label.replace(' ', '_')
        before = run(a.before, extra, work / f'{name}_before', work / 'cache')
        after = run(a.after, extra, work / f'{name}_after', work / 'cache')
        differ = differing(before, after)
        lines.append(f'{label}: {len(outputs(after))} output files; identical before and after: '
                     f'{"yes" if not differ else "NO, " + str(differ)}; log summary identical: '
                     f'{"yes" if summary(before) == summary(after) else "NO"}')
        print('\n'.join(lines), flush=True)
    lines.append('')

    Path(a.out).mkdir(parents=True, exist_ok=True)
    if a.genome_after:
        genome(a, lines)
    (Path(a.out) / 'verify_17.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--before', required=True)
    p.add_argument('--after', required=True)
    p.add_argument('--genome_before')
    p.add_argument('--genome_prototype')
    p.add_argument('--genome_after')
    p.add_argument('--work', default=None)
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'verify_17'))
    main_(p.parse_args())
