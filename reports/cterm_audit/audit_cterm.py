# -*- coding: utf-8 -*-

"""Audit of how alternative C-termini are tiered, and evaluation of an annotated-stop rule for them.

Every junction in tests/data/rmats is run through the shipped tiering logic with the corrected tx0/tx1 bounds
(reports/tx1_audit, PR 4). Each slice is then classified by where its first in-frame stop codon falls, and the
proposed rule is simulated:

    A slice whose first in-frame stop codon is an annotated stop_codon of a protein_coding transcript of the
    same gene ends there naturally: it is not a PTC, whether or not the stop is the slice's last codon. Such a
    "C-anchored" slice is stitched to the canonical at its N-terminus only. If its C-terminus is the canonical's
    it goes to psq_T1.fasta as today; otherwise to psq_T1_alt_cterm.fasta (or psq_T1_alt_cterm_orphan.fasta).
    Junctions qualify if both slices are either C-anchored or PTC-free; the T2 frameshift test is not applied
    when a slice is C-anchored, because slices ending at different stops differ in length arbitrarily.

Every stitched product is compared with the full set of proteins translated from the GTF's CDS records.
Outputs:

    slices.tsv      one row per slice: stop position and class, products and their validation
    junctions.tsv   one row per junction: actual fate (corrected bounds), slice states, proposed fate
    summary.txt     the aggregate numbers quoted in the report

Canonical sequences come from the GTF only and UniProt is stubbed out, so the run is deterministic and needs
no network.

Usage (from the repository root):
    python reports/cterm_audit/audit_cterm.py
"""

import argparse
import collections
import csv
import logging
import shutil
import sys
import tempfile
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'reports' / 'tx1_audit'))

import audit_tx1  # noqa: E402  (provides the corrected tx0/tx1 bounds)
from jcast import helpers as h, main, params  # noqa: E402
from jcast.annots import ReadAnnotations, ReadGenome  # noqa: E402
from jcast.junctions import Junction, RmatsResults  # noqa: E402
from jcast.sequences import Sequence  # noqa: E402

OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
STOPS = {'TAA', 'TAG', 'TGA'}
L = params.aa_stitch_length
READ, MASK = 1, 'NnXx'  # same settings as tests/test_pipeline.py and the tx1 audit
ACCEPTED = {'open', 'terminal', 'A'}

STATE_LABEL = {
    'open': 'no stop codon in slice',
    'terminal': 'stop is the last codon of the slice',
    'A': 'A: annotated protein_coding stop codon',
    'B': 'B: stop codon annotated only on non-protein_coding transcripts',
    'C': "C: unannotated stop in the slice's 3'-most exon",
    'D': "D: unannotated stop upstream of a junction in the slice",
}


# ---- annotation --------------------------------------------------------------------------------------------

def annotated_proteins(gtf, genome):
    """protein sequence -> sorted list of 'transcript_name:biotype' for every CDS-bearing transcript."""
    proteins = collections.defaultdict(set)
    for _, t in gtf.annot.query('feature == "CDS"').groupby('transcript_id'):
        t = t.sort_values('start')
        chrom = genome.genome[str(t['seqname'].iloc[0])].seq
        nt = ''.join(str(chrom[s - 1:e]) for s, e in zip(t['start'], t['end']))
        if t['strand'].iloc[0] == '-':
            nt, frame = str(Seq(nt).reverse_complement()), int(t['frame'].iloc[-1])
        else:
            frame = int(t['frame'].iloc[0])
        nt = nt[frame:]
        aa = str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip('*')
        proteins[aa].add(f"{t['transcript_name'].iloc[0]}:{t['transcript_biotype'].iloc[0]}")
    return {k: sorted(v) for k, v in proteins.items()}


def annotated_stop_sites(gtf):
    """gene_id -> {5'-most base of the stop codon in transcript orientation: set of transcript biotypes}."""
    sites = collections.defaultdict(lambda: collections.defaultdict(set))
    for _, r in gtf.annot.query('feature == "stop_codon"').iterrows():
        sites[r['gene_id']][r['start'] if r['strand'] == '+' else r['end']].add(r['transcript_biotype'])
    return sites


def product_class(product, proteins):
    if product is None:
        return 'unstitchable', ''
    matches = proteins.get(product)
    if not matches:
        return 'novel', ''
    if any(m.endswith(':protein_coding') for m in matches):
        return 'annotated protein_coding', ';'.join(matches)
    return 'annotated other biotype', ';'.join(matches)


# ---- slices ------------------------------------------------------------------------------------------------

def segments(jx, i):
    """Genomic (start, end) of the exons concatenated into slice i (0 or 1), in genomic order; see make_slice_localgenome."""
    t, s, alt = jx.junction_type, jx.strand, f'alt{i + 1}'
    if t in ('MXE', 'SE', 'RI'):
        names = ['anc', alt, 'down']
    elif t == 'A5SS':
        names = [alt, 'anc'] if s == '+' else ['anc', alt]
    else:
        names = ['anc', alt] if s == '+' else [alt, 'anc']
    return [(getattr(jx, f'{n}_es'), getattr(jx, f'{n}_ee')) for n in names
            if getattr(jx, f'{n}_es') > 0 and getattr(jx, f'{n}_ee') > 0]


def stitch(canonical, aa, n_only):
    """Sequence.stitch_to_canonical_aa on plain strings. n_only mirrors slice_has_ptc=True minus the both-ends try."""
    if not canonical or len(aa) < L:
        return None
    m0, m1 = canonical.find(aa[:L]), canonical.find(aa[-L:])
    if not n_only and m0 != -1 and m1 != -1:
        return canonical[:m0] + aa + canonical[m1 + L:]
    if n_only and m0 != -1:
        return canonical[:m0] + aa
    return None


def classify_slice(jx, nt_record, segs, sites):
    s = str(nt_record.seq).upper()
    if jx.strand == '-':
        s = h.get_complementary(s)
    positions = [c for es, ee in segs for c in range(es, ee + 1)]
    final = segs[-1] if jx.strand == '+' else segs[0]
    stop = next((k for k in range(jx.phase, len(s) - 2, 3) if s[k:k + 3] in STOPS), None)

    row = {'nt_len': len(s), 'stop_offset_nt': stop, 'stop_genomic': None, 'stop_biotypes': '',
           'stop_in_3prime_exon': None}
    if stop is None:
        state = 'open'
    else:
        genomic = positions[stop] if jx.strand == '+' else positions[len(positions) - 1 - stop]
        biotypes = sites.get(jx.gene_id, {}).get(genomic, set())
        in_final = stop >= len(s) - (final[1] - final[0] + 1)
        row.update(stop_genomic=genomic, stop_biotypes=';'.join(sorted(biotypes)), stop_in_3prime_exon=in_final)
        if stop > len(s) - 5:  # make_pep(terminate=True) accepts a stop in the last codon
            state = 'terminal'
        elif 'protein_coding' in biotypes:
            state = 'A'
        elif biotypes:
            state = 'B'
        elif in_final:
            state = 'C'
        else:
            state = 'D'
    row['state'] = state
    row['aa_terminate'] = h.make_pep(nt_record.seq, jx.strand, jx.phase, terminate=True)
    row['aa_to_stop'] = h.make_pep(nt_record.seq, jx.strand, jx.phase, terminate=False)
    return row


# ---- run ---------------------------------------------------------------------------------------------------

def main_():
    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))
    Junction._get_translated_region = audit_tx1.corrected_get_translated_region

    work = Path(tempfile.mkdtemp(prefix='cterm_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    proteins = annotated_proteins(gtf, genome)
    sites = annotated_stop_sites(gtf)
    args = argparse.Namespace(out=str(work), read=READ, model=False, canonical=True, mask=MASK, qvalue=[0, 1])

    slice_rows, junction_rows = [], []
    for jt, rma in zip(['MXE', 'SE', 'RI', 'A5SS', 'A3SS'], rmats.rmats_list):
        for i in range(len(rma)):
            record = rma.iloc[i].to_dict()
            fate = main._translate_one(Junction(**record), gtf, genome, args, pred_bound=READ).split('.')[0]

            jx = Junction(**record)
            jx.trim_cds(gtf)
            jx.get_translated_phase(gtf)
            base = {'type': jt, 'rmats_id': jx.name, 'gene_symbol': jx.gene_symbol, 'gene_id': jx.gene_id,
                    'strand': jx.strand, 'sum_sjc': jx.sum_sjc, 'phase': jx.phase, 'actual_fate': fate}
            if fate.startswith('SKIPPED') or jx.phase not in (0, 1, 2):
                junction_rows.append({**base, 'analysed': False, 'proposed_fate': fate})
                continue

            sq = Sequence(jx)
            sq.make_slice_localgenome(genome.genome)
            sq.get_canonical_aa(gtf=gtf, genome_index=genome.genome)
            canonical = str(sq.canonical_aa.seq)
            frameshift = sq.frameshift

            slices = []
            for k, nt in enumerate([sq.slice1_nt, sq.slice2_nt]):
                segs = segments(jx, k)
                r = classify_slice(jx, nt, segs, sites) if segs else {'state': 'empty', 'aa_terminate': '',
                                                                       'aa_to_stop': ''}
                r.update(base, slice=k + 1, frameshift=frameshift)

                # what the shipped stitching does with this slice: both ends if PTC-free, N-only if PTC (T2/T4)
                if r['aa_terminate']:
                    current = stitch(canonical, r['aa_terminate'], n_only=False)
                    m1 = canonical.find(r['aa_terminate'][-L:]) if current else -1
                    r['current_tail_after_stop'] = (r['state'] == 'terminal' and current is not None
                                                    and m1 + L < len(canonical))
                else:
                    current = stitch(canonical, r['aa_to_stop'], n_only=True)
                    r['current_tail_after_stop'] = False
                r['current_product_class'], r['current_product_matches'] = product_class(current, proteins)

                # the same slice under an N-only stitch, for every slice that has a stop codon
                has_stop = r['state'] in ('terminal', 'A', 'B', 'C', 'D')
                n_only = stitch(canonical, r['aa_to_stop'], n_only=True) if has_stop else None
                r['n_stitch_class'], r['n_stitch_matches'] = product_class(n_only, proteins) if has_stop else ('', '')
                r['c_anchored'] = r['state'] == 'A' or (r['state'] == 'terminal'
                                                        and 'protein_coding' in r['stop_biotypes'].split(';'))
                r['ends_at_canonical_cterm'] = (r['c_anchored'] and len(r['aa_to_stop']) >= L
                                                and canonical.endswith(r['aa_to_stop'][-L:]))
                # consistency: for a C-anchored slice that stitches today, the N-only product must be identical
                r['n_only_equals_current'] = (r['c_anchored'] and r['state'] == 'terminal' and current is not None
                                              and n_only == current)
                slices.append(r)

            states = [r['state'] for r in slices]
            eligible = fate == 'SUCCESS TIER 1' or (all(s in ACCEPTED for s in states) and 'A' in states)
            for r in slices:
                if not eligible:
                    r['proposed_file'] = ''
                elif r['c_anchored']:
                    if r['n_stitch_class'] == 'unstitchable':
                        r['proposed_file'] = 'psq_T1_alt_cterm_orphan'
                    elif r['ends_at_canonical_cterm']:
                        r['proposed_file'] = 'psq_T1'
                    else:
                        r['proposed_file'] = 'psq_T1_alt_cterm'
                else:
                    r['proposed_file'] = 'psq_T1' if r['current_product_class'] != 'unstitchable' else 'psq_T1_orphan'
                if fate == 'SUCCESS TIER 1' and r['state'] in ('open', 'terminal'):
                    r['current_file'] = 'psq_T1' if r['current_product_class'] != 'unstitchable' else 'psq_T1_orphan'
                else:
                    r['current_file'] = ''
            files = [r['proposed_file'] for r in slices]
            if not eligible:
                proposed = fate
            elif any(f.startswith('psq_T1_alt_cterm') for f in files):
                proposed = 'SUCCESS TIER 1 ALT CTERM'
            else:
                proposed = 'SUCCESS TIER 1'
            for r in slices:
                r['proposed_fate'] = proposed
            slice_rows.extend(slices)
            junction_rows.append({**base, 'analysed': True, 'slice1_state': states[0], 'slice2_state': states[1],
                                  'frameshift': frameshift, 'proposed_fate': proposed})

    write_tsv(OUT / 'slices.tsv', slice_rows, [
        'type', 'rmats_id', 'gene_symbol', 'gene_id', 'strand', 'slice', 'sum_sjc', 'phase', 'frameshift',
        'actual_fate', 'nt_len', 'state', 'stop_offset_nt', 'stop_genomic', 'stop_biotypes', 'stop_in_3prime_exon',
        'current_product_class', 'current_product_matches', 'current_tail_after_stop', 'current_file',
        'n_stitch_class', 'n_stitch_matches', 'c_anchored', 'ends_at_canonical_cterm', 'n_only_equals_current',
        'proposed_fate', 'proposed_file',
        'aa_to_stop'])
    write_tsv(OUT / 'junctions.tsv', junction_rows, [
        'type', 'rmats_id', 'gene_symbol', 'gene_id', 'strand', 'sum_sjc', 'phase', 'analysed', 'actual_fate',
        'slice1_state', 'slice2_state', 'frameshift', 'proposed_fate'])
    (OUT / 'summary.txt').write_text(summarise(slice_rows, junction_rows))
    print((OUT / 'summary.txt').read_text())


def write_tsv(path, rows, cols):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def table(title, counter):
    lines = [title] + [f'  {v:5d}  {" | ".join(map(str, k)) if isinstance(k, tuple) else k}'
                       for k, v in sorted(counter.items(), key=lambda kv: str(kv[0]))]
    return '\n'.join(lines) + '\n'


def summarise(slices, junctions):
    analysed = [j for j in junctions if j['analysed']]
    out = [f'junctions: {len(junctions)}; analysed (passed filters, GTF phase available): {len(analysed)}\n']
    out.append(table('actual fate of junctions not analysed', collections.Counter(
        j['actual_fate'] for j in junctions if not j['analysed'])))

    stopped = [s for s in slices if s['state'] in 'ABCD']
    out.append(table('slice stop class (analysed junctions)', collections.Counter(STATE_LABEL.get(s['state'], s['state'])
                                                                                  for s in slices)))
    out.append(table('stop class x product of an N-only stitch', collections.Counter(
        (STATE_LABEL[s['state']], s['n_stitch_class']) for s in stopped)))
    out.append(table('stop class x actual fate of the junction', collections.Counter(
        (STATE_LABEL[s['state']], s['actual_fate']) for s in stopped)))
    out.append(table('junction fate transitions (actual -> proposed)', collections.Counter(
        (j['actual_fate'], j['proposed_fate']) for j in analysed)))

    def validation(s):
        return s['n_stitch_class'] if s['c_anchored'] else s['current_product_class']

    moved = [s for s in slices if s['actual_fate'] != 'SUCCESS TIER 1' and s['proposed_fate'].startswith('SUCCESS TIER 1')]
    out.append(table('rescued junctions (not T1 today): slice proposed file x product validation', collections.Counter(
        (s['proposed_file'], validation(s)) for s in moved)))
    kept = [s for s in slices if s['actual_fate'] == 'SUCCESS TIER 1']
    out.append(table('T1 today: current file -> proposed file x product validation', collections.Counter(
        (s['current_file'], s['proposed_file'], validation(s)) for s in kept)))
    out.append(table('terminal-stop slices by stop biotype', collections.Counter(
        s['stop_biotypes'] or '(unannotated)' for s in slices if s['state'] == 'terminal')))
    anchored_stitched = [s for s in slices if s['c_anchored'] and s['state'] == 'terminal'
                         and s['current_product_class'] != 'unstitchable']
    out.append(f'consistency: C-anchored terminal slices that stitch today: {len(anchored_stitched)}; '
               f'N-only product identical: {sum(s["n_only_equals_current"] for s in anchored_stitched)}\n')
    out.append(table('proposed output file, all eligible slices', collections.Counter(
        s['proposed_file'] for s in slices if s['proposed_file'])))
    out.append(table('current output file, T1 slices today', collections.Counter(
        s['current_file'] for s in slices if s['current_file'])))
    t1 = [s for s in slices if s['actual_fate'] == 'SUCCESS TIER 1']
    out.append(table('baseline: slices of junctions that are SUCCESS TIER 1 today, by product validation',
                     collections.Counter(s['current_product_class'] for s in t1)))
    out.append(table('existing defect: PTC-free slice ending at a stop codon, stitched at both ends, with canonical '
                     'residues appended after the stop', collections.Counter(
                         (s['actual_fate'], s['gene_symbol'], s['type'], s['rmats_id'], s['slice'])
                         for s in slices if s.get('current_tail_after_stop'))))
    return '\n'.join(out)


if __name__ == '__main__':
    main_()
