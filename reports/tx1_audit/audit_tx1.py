# -*- coding: utf-8 -*-

"""Audit of the translation-bound (tx0/tx1) logic in Junction._get_translated_region.

Runs every junction in tests/data/rmats through the tiering logic twice — once with the shipped
_get_translated_region ("current") and once with stop-codon-based bounds ("corrected") — and writes:

    junctions.tsv     one row per junction: coordinates, bounds, trimming and fate under both variants
    fasta_diff.tsv    every output sequence present under one variant but not the other
    summary.txt       the aggregate numbers quoted in the report

Canonical sequences come from the GTF only and UniProt is stubbed out, so the run is deterministic and needs
no network. Tier assignment does not depend on the canonical source; only stitched-vs-orphan routing does.

Usage (from the repository root):
    python reports/tx1_audit/audit_tx1.py
"""

import argparse
import collections
import copy
import csv
import logging
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from jcast import main, params  # noqa: E402
from jcast.annots import ReadAnnotations, ReadGenome  # noqa: E402
from jcast.junctions import Junction, RmatsResults  # noqa: E402
from jcast.sequences import Sequence  # noqa: E402

OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
COORDS = ['anc_es', 'anc_ee', 'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'down_es', 'down_ee']

# Which coordinates each junction type trims at the lower (tx0) and upper (tx1) genomic bound; see Junction._trim.
LOWER_TRIMMED = {'MXE': ['anc_es', 'anc_ee', 'alt1_es', 'alt2_es'], 'SE': ['anc_es', 'anc_ee', 'alt1_es'],
                 'RI': ['anc_es'], 'A5SS': ['alt1_es', 'alt2_es'], 'A3SS': ['anc_es']}
UPPER_TRIMMED = {'MXE': ['down_es', 'down_ee', 'alt1_ee', 'alt2_ee'], 'SE': ['down_es', 'down_ee', 'alt1_ee'],
                 'RI': ['down_ee'], 'A5SS': ['anc_ee'], 'A3SS': ['alt1_ee', 'alt2_ee']}


def corrected_get_translated_region(self, gtf, startsite_index: int = 0):
    """Stop-codon-based bounds. tx0/tx1 are the lower/upper genomic bounds of the CDS, as _trim expects:

        + strand: tx0 = start_codon.start (most upstream first), tx1 = max(stop_codon.end)
        - strand: tx0 = min(stop_codon.start),                   tx1 = start_codon.end (most upstream first)

    Each bound falls back to None (no trimming on that side) independently when it cannot be determined.
    """
    gtf0 = gtf.annot.query('gene_id == @self.gene_id')

    def codons(feature):
        tsl = params.tsl_threshold  # must be local to this frame for the @tsl reference in query()
        coding = gtf0.query(f'feature == "{feature}" & transcript_biotype == "protein_coding"')
        if 'transcript_support_level' not in gtf0.columns:
            return coding
        return coding.query('transcript_support_level <= @tsl')

    starts, stops = codons('start_codon'), codons('stop_codon')
    self.num_start_codons = len(starts['start'].drop_duplicates())

    if self.strand == '+':
        five_prime = sorted(starts['start'].drop_duplicates())
        self.tx0 = five_prime[startsite_index] if len(five_prime) > startsite_index else None
        self.tx1 = int(stops['end'].max()) if len(stops) else None
    else:
        five_prime = sorted(starts['end'].drop_duplicates(), reverse=True)
        self.tx1 = five_prime[startsite_index] if len(five_prime) > startsite_index else None
        self.tx0 = int(stops['start'].min()) if len(stops) else None
    return True


def cds_facts(gtf, gene_id):
    """Annotated CDS span (incl. stop codon) and stop codon positions for the protein-coding transcripts used."""
    tsl = params.tsl_threshold
    g = gtf.annot.query('gene_id == @gene_id & transcript_biotype == "protein_coding" & '
                        'transcript_support_level <= @tsl')
    cds = g.query('feature == "CDS"')
    stops = g.query('feature == "stop_codon"')
    starts = g.query('feature == "start_codon"')
    return {
        'cds_min': int(min(cds['start'].min(), stops['start'].min() if len(stops) else cds['start'].min()))
        if len(cds) else '',
        'cds_max': int(max(cds['end'].max(), stops['end'].max() if len(stops) else cds['end'].max()))
        if len(cds) else '',
        'n_start_codons': len(starts['start'].drop_duplicates()),
        'n_stop_codons': len(stops['start'].drop_duplicates()),
        'cds_intervals': list(zip(cds['start'], cds['end'])),
    }


def inside_cds(pos, intervals):
    return pos is not None and any(s < pos < e for s, e in intervals)


def annotated_proteins(gtf, genome):
    """Every protein encoded by a CDS-bearing transcript in the GTF (any biotype, any TSL), translated directly."""
    proteins = set()
    for _, t in gtf.annot.query('feature == "CDS"').groupby('transcript_id'):
        t = t.sort_values('start')
        chrom = genome.genome[str(t['seqname'].iloc[0])].seq
        nt = ''.join(str(chrom[s - 1:e]) for s, e in zip(t['start'], t['end']))
        if t['strand'].iloc[0] == '-':
            nt, frame = str(Seq(nt).reverse_complement()), int(t['frame'].iloc[-1])
        else:
            frame = int(t['frame'].iloc[0])
        nt = nt[frame:]
        proteins.add(str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip('*'))
    return proteins


def exon_dropped(before, after):
    """True if trimming removed an exon entirely (coordinates set to -1)."""
    return any(after[c] == -1 and before[c] > 0 for c in COORDS)


def utr_exon_kept(row):
    """True if a retained exon lies entirely outside [tx0, tx1], i.e. is pure UTR that _trim left in place."""
    lo, hi = row['tx0'], row['tx1']
    if lo is None or hi is None:
        return False
    after = row['after']
    return any(after[f'{e}_es'] > 0 and (after[f'{e}_ee'] < lo or after[f'{e}_es'] > hi)
               for e in ['anc', 'alt1', 'alt2', 'down'])


def run(variant, gtf, genome, rmats, out_dir):
    if variant == 'corrected':
        Junction._get_translated_region = corrected_get_translated_region
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask='NnXx', qvalue=[0, 1])
    rows = {}
    for jt, rma in zip(['MXE', 'SE', 'RI', 'A5SS', 'A3SS'], rmats.rmats_list):
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            before = {c: getattr(jx, c) for c in COORDS}
            fate = main._translate_one(jx, gtf, genome, args, pred_bound=args.read).split('.')[0]
            after = {c: getattr(jx, c) for c in COORDS}
            rows[(jt, jx.name)] = {
                'jx': jx, 'before': before, 'after': after, 'fate': fate,
                'tx0': jx.tx0, 'tx1': jx.tx1, 'phase': jx.phase,
                'lower_trim': any(before[c] != after[c] for c in LOWER_TRIMMED[jt]),
                'upper_trim': any(before[c] != after[c] for c in UPPER_TRIMMED[jt]),
            }
    return rows


def read_fastas(out_dir):
    return {p.name: {str(r.seq): r.id for r in SeqIO.parse(p, 'fasta')} for p in sorted(out_dir.glob('psq_*.fasta'))}


def main_():
    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False

    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    work = Path(tempfile.mkdtemp(prefix='tx1_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')

    original = Junction._get_translated_region
    results, fastas = {}, {}
    for variant in ['current', 'corrected']:
        out_dir = work / variant
        out_dir.mkdir()
        results[variant] = run(variant, gtf, genome, rmats, out_dir)
        fastas[variant] = read_fastas(out_dir)
    Junction._get_translated_region = original

    cur, cor = results['current'], results['corrected']
    proteins = annotated_proteins(gtf, genome)
    facts_cache = {}

    with open(OUT / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['type', 'rmats_id', 'gene_id', 'gene_symbol', 'strand', 'sum_sjc', 'fdr',
                    *[f'input_{c}' for c in COORDS],
                    'cds_min', 'cds_max', 'n_start_codons', 'n_stop_codons',
                    'current_tx0', 'current_tx1', 'current_tx_lower_inside_cds', 'current_tx_upper_inside_cds',
                    'current_lower_trim', 'current_upper_trim', 'current_phase', 'current_fate',
                    *[f'current_{c}' for c in COORDS],
                    'corrected_tx0', 'corrected_tx1', 'corrected_lower_trim', 'corrected_upper_trim',
                    'corrected_phase', 'corrected_fate', *[f'corrected_{c}' for c in COORDS],
                    'current_exon_dropped', 'corrected_exon_dropped', 'corrected_utr_exon_kept',
                    'fate_changed', 'coords_changed'])
        for key in cur:
            a, b = cur[key], cor[key]
            jx = a['jx']
            if jx.gene_id not in facts_cache:
                facts_cache[jx.gene_id] = cds_facts(gtf, jx.gene_id)
            f = facts_cache[jx.gene_id]
            w.writerow([key[0], key[1], jx.gene_id, jx.gene_symbol, jx.strand, jx.sum_sjc, jx.fdr,
                        *[a['before'][c] for c in COORDS],
                        f['cds_min'], f['cds_max'], f['n_start_codons'], f['n_stop_codons'],
                        a['tx0'], a['tx1'], inside_cds(a['tx0'], f['cds_intervals']),
                        inside_cds(a['tx1'], f['cds_intervals']),
                        a['lower_trim'], a['upper_trim'], a['phase'], a['fate'], *[a['after'][c] for c in COORDS],
                        b['tx0'], b['tx1'], b['lower_trim'], b['upper_trim'], b['phase'], b['fate'],
                        *[b['after'][c] for c in COORDS],
                        exon_dropped(a['before'], a['after']), exon_dropped(b['before'], b['after']), utr_exon_kept(b),
                        a['fate'] != b['fate'], a['after'] != b['after']])

    with open(OUT / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['file', 'present_in', 'matches_annotated_protein', 'header', 'sequence'])
        for name in sorted(set(fastas['current']) | set(fastas['corrected'])):
            A, B = fastas['current'].get(name, {}), fastas['corrected'].get(name, {})
            for seq in sorted(set(A) - set(B)):
                w.writerow([name, 'current_only', seq in proteins, A[seq], seq])
            for seq in sorted(set(B) - set(A)):
                w.writerow([name, 'corrected_only', seq in proteins, B[seq], seq])

    lines = []
    p = lines.append
    p(f'junctions: {len(cur)}')
    for variant, rows in results.items():
        p(f'\n[{variant}] fates: {dict(sorted(collections.Counter(r["fate"] for r in rows.values()).items()))}')
        by_type = collections.defaultdict(lambda: [0, 0, 0])
        for (jt, _), r in rows.items():
            by_type[jt][0] += 1
            by_type[jt][1] += r['lower_trim']
            by_type[jt][2] += r['upper_trim']
        p(f'[{variant}] per type [n, lower-bound trims, upper-bound trims]: {dict(by_type)}')
    p(f'\ncoordinates differ after trimming: {sum(cur[k]["after"] != cor[k]["after"] for k in cur)}')
    trans = collections.Counter((cur[k]['fate'], cor[k]['fate']) for k in cur if cur[k]['fate'] != cor[k]['fate'])
    p(f'fate changed: {sum(trans.values())}')
    for (x, y), n in trans.most_common():
        p(f'  {n:3d}  {x} -> {y}')
    success = ('SUCCESS TIER 1', 'SUCCESS TIER 2', 'SUCCESS TIER 3')
    for variant, rows in results.items():
        ok = [r for r in rows.values() if r['fate'] in success]
        p(f'[{variant}] T1-T3 successes whose trimming removed a whole exon (artefact or legitimate UTR removal): '
          f'{sum(exon_dropped(r["before"], r["after"]) for r in ok)} / {len(ok)}')

    p(f'[corrected] junctions retaining a pure-UTR exon outside [tx0, tx1]: '
      f'{sum(utr_exon_kept(r) for r in cor.values())}  (not meaningful for current: its bounds are degenerate)')

    p(f'\nannotated proteins in GTF (truth set): {len(proteins)}')
    p('records exactly matching an annotated protein [current, corrected]:')
    for name in ['psq_T1.fasta', 'psq_T2.fasta', 'psq_T4.fasta', 'psq_canonical.fasta']:
        counts = []
        for variant in ['current', 'corrected']:
            seqs = fastas[variant].get(name, {})
            counts.append(f'{sum(q in proteins for q in seqs)}/{len(seqs)}')
        p(f'  {name:24s} {counts[0]:>9s} {counts[1]:>9s}')
    new_t1 = set(fastas['corrected']['psq_T1.fasta']) - set(fastas['current']['psq_T1.fasta'])
    p(f'  T1 sequences gained by the correction that match an annotated protein: '
      f'{sum(q in proteins for q in new_t1)}/{len(new_t1)}')

    p('\nFASTA files [current records, corrected records, current-only seqs, corrected-only seqs]:')
    for name in sorted(set(fastas['current']) | set(fastas['corrected'])):
        A, B = fastas['current'].get(name, {}), fastas['corrected'].get(name, {})
        p(f'  {name:24s} {len(A):4d} {len(B):4d} {len(set(A) - set(B)):4d} {len(set(B) - set(A)):4d}')
    (OUT / 'summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    shutil.rmtree(work)


if __name__ == '__main__':
    main_()
