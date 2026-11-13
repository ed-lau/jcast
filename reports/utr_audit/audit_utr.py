# -*- coding: utf-8 -*-

"""Audit of PR 4a: trimming every exon to the coding region and skipping events that lie outside it.

Runs every junction in tests/data/rmats through main._translate_one twice, each in its own process:

    baseline    a git worktree of BASELINE_COMMIT (PR 4: corrected tx0/tx1, per-type _trim)
    candidate   this repository's working tree (PR 4a: strand-agnostic _trim, SKIPPED_OUTSIDE_CDS)

and writes:

    junctions.tsv     one row per junction: bounds, trimmed coordinates, phase and fate under both versions,
                      and whether exons the candidate removes overlap any annotated CDS
    fasta_diff.tsv    every output sequence present under one version but not the other
    summary.txt       the aggregate numbers quoted in the report

Canonical sequences come from the GTF only and UniProt is stubbed out, so the run is deterministic and needs
no network.

Usage (from the repository root):
    python reports/utr_audit/audit_utr.py [--baseline COMMIT]
"""

import argparse
import collections
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
BASELINE_COMMIT = 'da72a0a'
EXONS = ['anc', 'alt1', 'alt2', 'down']
COORDS = [f'{e}_{x}' for e in EXONS for x in ('es', 'ee')]
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']


# ---- worker: runs inside one code version ------------------------------------------------------------------

def worker(out_dir):
    """Run the pipeline with whichever jcast is first on sys.path; dump per-junction results as JSON."""
    import logging
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    # record whether trim_cds falls back to the second start site
    original = Junction._get_translated_region

    def recording(self, gtf, startsite_index=0):
        if startsite_index == 1:
            self.retry_fired = True
        return original(self, gtf, startsite_index=startsite_index)
    Junction._get_translated_region = recording

    out_dir = Path(out_dir)
    work = Path(tempfile.mkdtemp(prefix='utr_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask='NnXx', qvalue=[0, 1])

    rows = []
    for jt, rma in zip(TYPES, rmats.rmats_list):
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            before = {c: int(getattr(jx, c)) for c in COORDS}
            fate = main._translate_one(jx, gtf, genome, args, pred_bound=args.read)
            rows.append({
                'type': jt, 'rmats_id': jx.name, 'gene_id': jx.gene_id, 'gene_symbol': jx.gene_symbol,
                'strand': jx.strand, 'sum_sjc': int(jx.sum_sjc), 'fate': short_fate(fate),
                'tx0': jx.tx0, 'tx1': jx.tx1, 'phase': jx.phase, 'retry_fired': getattr(jx, 'retry_fired', False),
                'before': before, 'after': {c: int(getattr(jx, c)) for c in COORDS},
            })
    (out_dir / 'junctions.json').write_text(json.dumps(rows, default=int))


def short_fate(fate):
    head = fate.strip().split('\n')[0]
    if head.startswith('SKIPPED'):
        return head.split('. ')[0] + ': ' + head.split('. ')[1].rstrip('.')
    return head.split('.')[0]


def run_variant(repo, out_dir):
    out_dir.mkdir(parents=True)
    env = dict(os.environ, PYTHONPATH=str(repo), MPLBACKEND='Agg')
    subprocess.run([sys.executable, __file__, '--worker', str(out_dir)], env=env, check=True, cwd=str(out_dir),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return {(r['type'], r['rmats_id']): r for r in json.loads((out_dir / 'junctions.json').read_text())}


# ---- comparison --------------------------------------------------------------------------------------------

def any_cds_intervals(gtf_path):
    """gene_id -> list of (start, end, transcript_name, biotype, tsl) for every CDS row, any transcript."""
    import pandas as pd
    from gtfparse import read_gtf
    gtf = read_gtf(str(gtf_path))
    gtf = pd.DataFrame(gtf) if not isinstance(gtf, pd.DataFrame) else gtf
    cds = gtf[gtf['feature'] == 'CDS']
    out = collections.defaultdict(list)
    for r in cds.itertuples():
        out[r.gene_id].append((int(r.start), int(r.end), r.transcript_name, r.transcript_biotype,
                               str(getattr(r, 'transcript_support_level', ''))))
    return out


def exon_changes(before, after):
    dropped = [e for e in EXONS if after[f'{e}_es'] == -1 and before[f'{e}_es'] > 0]
    trimmed = [e for e in EXONS if e not in dropped and
               (after[f'{e}_es'], after[f'{e}_ee']) != (before[f'{e}_es'], before[f'{e}_ee'])]
    return dropped, trimmed


def utr_exon_kept(row):
    lo, hi = row['tx0'], row['tx1']
    return any(row['after'][f'{e}_es'] > 0 and row['after'][f'{e}_ee'] > 0 and
               ((lo is not None and row['after'][f'{e}_ee'] < lo) or (hi is not None and row['after'][f'{e}_es'] > hi))
               for e in EXONS)


def read_fastas(out_dir):
    from Bio import SeqIO
    return {p.name: {str(r.seq): r.id for r in SeqIO.parse(p, 'fasta')} for p in sorted(out_dir.glob('psq_*.fasta'))}


def table(title, counter):
    lines = [title] + [f'  {v:5d}  {" | ".join(map(str, k)) if isinstance(k, tuple) else k}'
                       for k, v in sorted(counter.items(), key=lambda kv: str(kv[0]))]
    return '\n'.join(lines) + '\n'


def main_(baseline_commit):
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / 'reports' / 'tx1_audit'))
    work = Path(tempfile.mkdtemp(prefix='utr_audit_'))
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), baseline_commit],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        base = run_variant(baseline_repo, work / 'out_baseline')
        cand = run_variant(REPO, work / 'out_candidate')
    finally:
        subprocess.run(['git', '-C', str(REPO), 'worktree', 'remove', '--force', str(baseline_repo)], check=True)

    import logging
    import audit_tx1
    from jcast.annots import ReadAnnotations, ReadGenome
    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    proteins = audit_tx1.annotated_proteins(gtf, genome)
    cds = collections.defaultdict(list)
    for r in gtf.annot.query('feature == "CDS"').itertuples():
        cds[r.gene_id].append((int(r.start), int(r.end), r.transcript_name, r.transcript_biotype,
                               r.transcript_support_level))

    rows, out = [], []
    for key, b in base.items():
        c = cand[key]
        dropped, trimmed = exon_changes(b['after'], c['after'])
        coding_hits = sorted({f'{t}:{bt}:TSL{tsl}' for e in dropped
                              for s, t_end, t, bt, tsl in cds.get(c['gene_id'], [])
                              if s <= b['after'][f'{e}_ee'] and t_end >= b['after'][f'{e}_es']})
        rows.append({
            'type': key[0], 'rmats_id': key[1], 'gene_symbol': c['gene_symbol'], 'gene_id': c['gene_id'],
            'strand': c['strand'], 'sum_sjc': c['sum_sjc'], 'tx0': c['tx0'], 'tx1': c['tx1'],
            'baseline_fate': b['fate'], 'candidate_fate': c['fate'],
            'exons_dropped': '+'.join(dropped), 'exons_trimmed': '+'.join(trimmed),
            'dropped_overlaps_cds': ';'.join(coding_hits),
            'baseline_utr_exon_kept': utr_exon_kept(b), 'candidate_utr_exon_kept': utr_exon_kept(c),
            'baseline_phase': b['phase'], 'candidate_phase': c['phase'],
            'baseline_retry': b['retry_fired'], 'candidate_retry': c['retry_fired'],
            **{f'baseline_{k}': v for k, v in b['after'].items()},
            **{f'candidate_{k}': v for k, v in c['after'].items()},
        })

    with open(OUT / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    fb, fc = read_fastas(work / 'out_baseline'), read_fastas(work / 'out_candidate')
    diff = []
    for name in sorted(set(fb) | set(fc)):
        b_seqs, c_seqs = fb.get(name, {}), fc.get(name, {})
        for seq in sorted(set(b_seqs) - set(c_seqs)):
            diff.append({'file': name, 'change': 'lost', 'header': b_seqs[seq], 'length': len(seq),
                         'matches_annotated_protein': seq in proteins, 'sequence': seq})
        for seq in sorted(set(c_seqs) - set(b_seqs)):
            diff.append({'file': name, 'change': 'gained', 'header': c_seqs[seq], 'length': len(seq),
                         'matches_annotated_protein': seq in proteins, 'sequence': seq})
    with open(OUT / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'change', 'header', 'length', 'matches_annotated_protein',
                                           'sequence'], delimiter='\t')
        w.writeheader()
        w.writerows(diff)

    changed = [r for r in rows if r['baseline_fate'] != r['candidate_fate']]
    out.append(f'baseline commit: {baseline_commit}; junctions: {len(rows)}; fate changed: {len(changed)}\n')
    out.append(table('fate counts: baseline', collections.Counter(r['baseline_fate'] for r in rows)))
    out.append(table('fate counts: candidate', collections.Counter(r['candidate_fate'] for r in rows)))
    out.append(table('fate transitions (baseline | candidate)', collections.Counter(
        (r['baseline_fate'], r['candidate_fate']) for r in changed)))
    out.append(table('coordinate changes (type | strand | dropped | trimmed)', collections.Counter(
        (r['type'], r['strand'], r['exons_dropped'] or '-', r['exons_trimmed'] or '-')
        for r in rows if r['exons_dropped'] or r['exons_trimmed'])))
    out.append(table('junctions keeping an exon entirely outside [tx0, tx1] (baseline, candidate)', collections.Counter(
        (r['baseline_utr_exon_kept'], r['candidate_utr_exon_kept']) for r in rows)))
    hits = [r for r in rows if r['dropped_overlaps_cds']]
    out.append(f'junctions with a dropped exon overlapping any annotated CDS: {len(hits)}\n')
    out.append(table('  by overlapping transcript (name:biotype:TSL)', collections.Counter(
        h for r in hits for h in r['dropped_overlaps_cds'].split(';'))))
    out.append(table('fate transitions of junctions with a dropped exon overlapping a CDS', collections.Counter(
        (r['gene_symbol'], r['type'], r['rmats_id'], r['baseline_fate'], r['candidate_fate']) for r in hits)))
    out.append(table('second-start-site retry fired (baseline, candidate)', collections.Counter(
        (r['baseline_retry'], r['candidate_retry']) for r in rows)))
    out.append(table('GTF phase available (baseline, candidate)', collections.Counter(
        (r['baseline_phase'] is not None, r['candidate_phase'] is not None) for r in rows)))
    out.append(table('FASTA records (file | baseline | candidate)', collections.Counter(
        {(n, len(fb.get(n, {})), len(fc.get(n, {}))): 1 for n in sorted(set(fb) | set(fc))})))
    out.append(table('FASTA diff (file | change | matches annotated protein)', collections.Counter(
        (d['file'], d['change'], d['matches_annotated_protein']) for d in diff)))
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--worker':
        worker(sys.argv[2])
    else:
        parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
        parser.add_argument('--baseline', default=BASELINE_COMMIT, help='commit to compare against')
        main_(parser.parse_args().baseline)
