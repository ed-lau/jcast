# -*- coding: utf-8 -*-

"""Audit for PR 4c: alternative C- and N-termini, evaluated on the current pipeline (after PR 4b).

Every junction in tests/data/rmats is run through the shipped main._translate_one, and its actual fate and FASTA
output are recorded. Each slice is then classified at both ends, and the proposed routing is simulated:

    C-anchored  the slice's first in-frame stop codon is an annotated stop_codon of a protein_coding transcript of
                the gene (internal or in the last codon). The slice ends there; it is not a PTC.
    N-anchored  the slice lacks the gene's start codon (the most upstream one, as used for the coding bounds) while
                the other slice contains it, and it contains another annotated protein-coding start codon. The slice
                is translated from that start (in the frame the start defines).

    product     C-anchored only: joined to the canonical at its N-terminus only
                N-anchored only: joined at its C-terminus only
                both:            the slice itself, start codon to stop codon (no joining)
                neither:         as today (both ends)
    file        by which terminus differs from the canonical: psq_T1 (neither), psq_T1_alt_cterm, psq_T1_alt_nterm,
                psq_T1_alt_both; unjoinable anchored slices go to the matching *_orphan file

A junction is routed this way if at least one slice is anchored, both slices translate to an acceptable end
(no stop codon, a last-codon stop, or an annotated protein-coding stop), and it is not a frameshift. The frameshift
test is skipped only when a slice ends at an internal annotated stop (the slices then end at different stops, so their
lengths differ arbitrarily); otherwise it compares the slices' lengths from their translation starts (frame or
alternative start codon) to the slice end. Otherwise its fate is unchanged. Annotated
sites are taken at two TSL policies: 'tsl2' (TSL <= params.tsl_threshold, as for the coding bounds) and 'any'.

Every product is compared with all proteins translated from the GTF's CDS records. Outputs:

    slices.tsv      one row per slice and policy: frames, start/stop classes, current and proposed product and file
    junctions.tsv   one row per junction and policy: actual fate, slice states, proposed fate
    summary.txt     the aggregate numbers quoted in the report

Usage (from the repository root):
    python reports/termini_audit/audit_termini.py
"""

import argparse
import collections
import csv
import logging
import re
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'reports' / 'tx1_audit'))
sys.path.insert(0, str(REPO / 'reports' / 'cterm_audit'))

from jcast import helpers as h, main, params  # noqa: E402
from jcast.annots import ReadAnnotations, ReadGenome  # noqa: E402
from jcast.junctions import Junction, RmatsResults  # noqa: E402
from jcast.sequences import Sequence  # noqa: E402

OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']
STOPS = {'TAA', 'TAG', 'TGA'}
L = params.aa_stitch_length
POLICIES = ['tsl2', 'any']
ACCEPTABLE_END = {'open', 'terminal', 'A'}
HEADER = re.compile(r'\|(MXE|SE|RI|A5SS|A3SS)([12])\|(\d+)\|')


# ---- annotation ------------------------------------------------------------------------------------------------

def codon_sites(gtf, feature):
    """policy -> gene_id -> {5' base of the codon in transcript orientation: set of transcript biotypes}."""
    rows = gtf.annot.query(f'feature == "{feature}"')
    out = {p: collections.defaultdict(lambda: collections.defaultdict(set)) for p in POLICIES}
    for r in rows.itertuples():
        site = r.start if r.strand == '+' else r.end
        for p in POLICIES:
            if p == 'tsl2' and not r.transcript_support_level <= params.tsl_threshold:
                continue
            out[p][r.gene_id][site].add(r.transcript_biotype)
    return out


def classify(product, proteins):
    if product is None:
        return 'unjoined'
    m = proteins.get(product, [])
    if not m:
        return 'new'
    return 'annotated protein_coding' if any(x.endswith(':protein_coding') for x in m) else 'annotated other'


# ---- stitching (plain strings; cf. Sequence.stitch_to_canonical_aa) -----------------------------------------------

def stitch(canonical, aa, mode):
    if not canonical or len(aa) < L:
        return None
    m0, m1 = canonical.find(aa[:L]), canonical.find(aa[-L:])
    if mode == 'both':
        return canonical[:m0] + aa + canonical[m1 + L:] if m0 != -1 and m1 != -1 else None
    if mode == 'n':
        return canonical[:m0] + aa if m0 != -1 else None
    if mode == 'c':
        return aa + canonical[m1 + L:] if m1 != -1 else None
    return aa  # 'none': self-contained


def target_file(product, canonical, anchored):
    if product is None:
        return 'psq_T1_alt_orphan' if anchored else 'psq_T1_orphan'
    if product == canonical:
        return 'psq_canonical'
    n_same, c_same = product[:L] == canonical[:L], product[-L:] == canonical[-L:]
    return {(True, True): 'psq_T1', (True, False): 'psq_T1_alt_cterm',
            (False, True): 'psq_T1_alt_nterm', (False, False): 'psq_T1_alt_both'}[(n_same, c_same)]


# ---- per slice ---------------------------------------------------------------------------------------------------

def slice_facts(jx, k, nt_record, exons, frame, stops, starts, gene_start):
    """Everything about one slice that does not depend on the other slice."""
    s = str(nt_record.seq).upper()
    if jx.strand == '-':
        s = h.get_complementary(s)
    positions = [c for es, ee in exons for c in range(es, ee + 1)]
    assert len(positions) == len(s), (jx.junction_type, jx.name, k)

    def genomic(offset):
        return positions[offset] if jx.strand == '+' else positions[len(positions) - 1 - offset]

    def offset(pos):
        i = positions.index(pos)
        return i if jx.strand == '+' else len(positions) - 1 - i

    contains_gene_start = gene_start is not None and gene_start in set(positions)
    # annotated protein-coding start codons inside the slice, as offsets; only real ATGs at that offset
    inside = [offset(p) for p, bt in starts.items() if 'protein_coding' in bt and p in set(positions)]
    alt_starts = sorted(o for o in inside if s[o:o + 3] == 'ATG')
    return {'s': s, 'genomic': genomic, 'frame': frame, 'contains_gene_start': contains_gene_start,
            'alt_starts': alt_starts, 'non_atg_starts': len(inside) - len(alt_starts), 'stops': stops}


def translate_from(f, start, nt_record, strand):
    """Translation from offset `start` (a frame or an ATG offset): first stop, stop class, products."""
    s = f['s']
    stop = next((i for i in range(start, len(s) - 2, 3) if s[i:i + 3] in STOPS), None)
    if stop is None:
        state, biotypes = 'open', set()
    else:
        biotypes = f['stops'].get(f['genomic'](stop), set())
        if stop > len(s) - 5:
            state = 'terminal'
        elif 'protein_coding' in biotypes:
            state = 'A'
        elif biotypes:
            state = 'B'
        else:
            state = 'D'  # unannotated stop
    c_anchored = state == 'A' or (state == 'terminal' and 'protein_coding' in biotypes)
    aa = h.make_pep(nt_record.seq, strand, start, terminate=False)  # make_pep accepts any start offset
    return {'stop_state': state, 'c_anchored': c_anchored, 'aa': aa}


# ---- run -----------------------------------------------------------------------------------------------------------

def main_():
    import audit_cterm  # annotated_proteins with transcript names and biotypes

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    work = Path(tempfile.mkdtemp(prefix='termini_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    proteins = audit_cterm.annotated_proteins(gtf, genome)
    stop_sites, start_sites = codon_sites(gtf, 'stop_codon'), codon_sites(gtf, 'start_codon')

    out_dir = work / 'current'
    out_dir.mkdir()
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask='NnXx', qvalue=[0, 1])

    slice_rows, junction_rows = [], []
    for jt, rma in zip(TYPES, rmats.rmats_list):
        for i in range(len(rma)):
            record = rma.iloc[i].to_dict()
            fate = main._translate_one(Junction(**record), gtf, genome, args, pred_bound=args.read).split('.')[0]

            jx = Junction(**record)
            jx.trim_cds(gtf)
            jx.get_translated_phase(gtf)
            base = {'type': jt, 'rmats_id': jx.name, 'gene_symbol': jx.gene_symbol, 'strand': jx.strand,
                    'sum_sjc': jx.sum_sjc, 'actual_fate': fate}
            if fate == 'SKIPPED':
                junction_rows.append({**base, 'policy': '-', 'analysed': False, 'proposed_fate': fate})
                continue

            sq = Sequence(jx)
            sq.make_slice_localgenome(genome.genome)
            sq.get_canonical_aa(gtf=gtf, genome_index=genome.genome)
            canonical = str(sq.canonical_aa.seq)
            gene_start = jx.tx0 if jx.strand == '+' else jx.tx1
            exons = jx.slice_exons()

            for policy in POLICIES:
                facts = [slice_facts(jx, k, nt, exons[k], jx.slice_phases[k],
                                     stop_sites[policy].get(jx.gene_id, {}), start_sites[policy].get(jx.gene_id, {}),
                                     gene_start)
                         for k, nt in enumerate([sq.slice1_nt, sq.slice2_nt])]
                start_less = [not f['contains_gene_start'] and other['contains_gene_start']
                              for f, other in zip(facts, facts[::-1])]
                rows = []
                for k, (f, nt) in enumerate(zip(facts, [sq.slice1_nt, sq.slice2_nt])):
                    n_anchored = start_less[k] and bool(f['alt_starts'])
                    if n_anchored:
                        t = translate_from(f, f['alt_starts'][0], nt, jx.strand)
                    elif f['frame'] in (0, 1, 2):
                        t = translate_from(f, f['frame'], nt, jx.strand)
                    else:
                        t = {'stop_state': 'no frame', 'c_anchored': False, 'aa': ''}
                    mode = {(True, True): 'none', (True, False): 'c', (False, True): 'n',
                            (False, False): 'both'}[(n_anchored, t['c_anchored'])]
                    product = stitch(canonical, t['aa'], mode) if t['aa'] else None
                    anchored = n_anchored or t['c_anchored']

                    # today's rule on the same slice: its frame; both ends if it reads through, else N-only (T4)
                    today = translate_from(f, f['frame'], nt, jx.strand) if f['frame'] in (0, 1, 2) else None
                    if today is None or not today['aa']:
                        today_product = None
                    elif today['stop_state'] in ('open', 'terminal'):
                        today_product = stitch(canonical, today['aa'], 'both')
                    else:
                        today_product = stitch(canonical, today['aa'], 'n')
                    n_only = stitch(canonical, t['aa'], 'n') if t['aa'] else None
                    tx_start = f['alt_starts'][0] if n_anchored else f['frame']
                    rows.append({**base, 'policy': policy, 'slice': k + 1, 'frame': f['frame'],
                                 'coding_nt': len(f['s']) - tx_start if tx_start is not None else None,
                                 'contains_gene_start': f['contains_gene_start'], 'start_less': start_less[k],
                                 'alt_starts_in_slice': len(f['alt_starts']), 'n_anchored': n_anchored,
                                 'stop_state': t['stop_state'], 'c_anchored': t['c_anchored'], 'join_mode': mode,
                                 'aa_len': len(t['aa']), 'product_class': classify(product, proteins),
                                 'non_atg_starts': f['non_atg_starts'],
                                 'today_stop_state': today['stop_state'] if today else 'no frame',
                                 'today_product_class': classify(today_product, proteins) if today else 'no frame',
                                 'today_is_canonical': today_product == canonical if today_product else False,
                                 'n_only_class': classify(n_only, proteins),
                                 'product_matches': ';'.join(proteins.get(product, [])) if product else '',
                                 'proposed_file_if_routed': target_file(product, canonical, anchored),
                                 'aa': t['aa']})

                anchored_any = any(r['n_anchored'] or r['c_anchored'] for r in rows)
                acceptable = all(r['stop_state'] in ACCEPTABLE_END for r in rows)
                internal_stop = any(r['stop_state'] == 'A' for r in rows)
                frameshift = acceptable and not internal_stop and \
                    (rows[0]['coding_nt'] - rows[1]['coding_nt']) % 3 != 0
                routed = anchored_any and acceptable and not frameshift
                for r in rows:
                    r['frameshift_blocks'] = anchored_any and acceptable and frameshift
                if routed:
                    alt = any(r['proposed_file_if_routed'].startswith('psq_T1_alt') for r in rows)
                    proposed = 'SUCCESS TIER 1 ALT TERMINUS' if alt else 'SUCCESS TIER 1'
                else:
                    proposed = fate
                for r in rows:
                    r['routed'] = routed
                    r['proposed_fate'] = proposed
                    r['proposed_file'] = r['proposed_file_if_routed'] if routed else ''
                slice_rows.extend(rows)
                junction_rows.append({**base, 'policy': policy, 'analysed': True,
                                      'slice1': f"{rows[0]['stop_state']}/{'N' if rows[0]['n_anchored'] else '-'}",
                                      'slice2': f"{rows[1]['stop_state']}/{'N' if rows[1]['n_anchored'] else '-'}",
                                      'routed': routed, 'proposed_fate': proposed})

    # current output, by junction: which files each junction's records went to, and annotation
    current = collections.defaultdict(list)
    for p in sorted(out_dir.glob('psq_*.fasta')):
        for r in SeqIO.parse(p, 'fasta'):
            m = HEADER.search(r.id)
            if m:
                current[(m.group(1), m.group(3), int(m.group(2)))].append((p.stem, classify(str(r.seq), proteins)))
    for r in slice_rows:
        recs = current.get((r['type'], r['rmats_id'], r['slice']), [])
        r['current_files'] = ';'.join(sorted({f for f, _ in recs}))
        r['current_class'] = ';'.join(sorted({c for _, c in recs}))

    cols = ['policy', 'type', 'rmats_id', 'gene_symbol', 'strand', 'slice', 'sum_sjc', 'actual_fate', 'frame',
            'contains_gene_start', 'start_less', 'alt_starts_in_slice', 'n_anchored', 'stop_state', 'c_anchored',
            'join_mode', 'aa_len', 'coding_nt', 'frameshift_blocks', 'routed', 'proposed_fate', 'proposed_file',
            'product_class', 'product_matches',
            'today_stop_state', 'today_product_class', 'today_is_canonical', 'n_only_class', 'non_atg_starts',
            'current_files', 'current_class', 'aa']
    with open(OUT / 'slices.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(slice_rows)
    with open(OUT / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['policy', 'type', 'rmats_id', 'gene_symbol', 'strand', 'sum_sjc',
                                           'analysed', 'actual_fate', 'slice1', 'slice2', 'routed', 'proposed_fate'],
                           delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(junction_rows)
    (OUT / 'summary.txt').write_text(summarise(slice_rows, junction_rows, current))
    print((OUT / 'summary.txt').read_text())


def table(title, counter):
    lines = [title] + [f'  {v:5d}  {" | ".join(map(str, k)) if isinstance(k, tuple) else k}'
                       for k, v in sorted(counter.items(), key=lambda kv: str(kv[0]))]
    return '\n'.join(lines) + '\n'


def summarise(slices, junctions, current):
    out = [f'junctions: {len({(j["type"], j["rmats_id"]) for j in junctions})}; analysed (not skipped): '
           f'{len({(j["type"], j["rmats_id"]) for j in junctions if j["analysed"]})}\n']
    out.append(table('current output records by file and annotation (all junctions)', collections.Counter(
        rec for recs in current.values() for rec in recs)))
    for policy in POLICIES:
        S = [s for s in slices if s['policy'] == policy]
        J = [j for j in junctions if j['policy'] == policy]
        out.append(f'=================== policy: {policy} ===================\n')
        out.append(table('slice: stop class (from its frame, or from the alternative start if N-anchored)',
                         collections.Counter(s['stop_state'] for s in S)))
        out.append(table('slice: start codon | N-anchored', collections.Counter(
            ('contains gene start' if s['contains_gene_start'] else
             ('start-less' if s['start_less'] else 'gene start upstream / absent'),
             'alt start in slice' if s['alt_starts_in_slice'] else '-', s['n_anchored']) for s in S)))
        out.append(table('anchored slices: join mode x product class', collections.Counter(
            (s['join_mode'], s['product_class']) for s in S if s['n_anchored'] or s['c_anchored'])))
        out.append(table('C-anchored internal stop (class A): product class', collections.Counter(
            s['product_class'] for s in S if s['stop_state'] == 'A')))
        out.append(table('stop class x product of an N-only join (A routed; B = NMD-only stop, D = unannotated: '
                         'never routed, controls)', collections.Counter(
            (s['stop_state'], s['n_only_class']) for s in S if s['stop_state'] in ('A', 'B', 'D'))))
        out.append(f'annotated start codons inside slices that are not ATG in the genome: '
                   f'{sum(s["non_atg_starts"] for s in S)}\n')
        out.append(table('anchored junctions not routed because they are frameshifts (gene | type | id | fate)',
                         collections.Counter((s['gene_symbol'], s['type'], s['rmats_id'], s['actual_fate'])
                                             for s in S if s['frameshift_blocks'] and s['slice'] == 1)))
        out.append(table('junction fate: actual -> proposed', collections.Counter(
            (j['actual_fate'], j['proposed_fate']) for j in J if j['analysed'])))
        routed = [s for s in S if s['routed']]
        out.append(table('routed slices: proposed file x product class', collections.Counter(
            (s['proposed_file'], s['product_class']) for s in routed)))
        out.append(table('routed slices: current files -> proposed file (changed only)', collections.Counter(
            (s['current_files'] or '(none)', s['proposed_file']) for s in routed
            if s['proposed_file'] not in s['current_files'].split(';'))))
        t1 = [s for s in S if s['actual_fate'] == 'SUCCESS TIER 1']
        out.append(table('baseline: slices of junctions that are T1 today, product under today\'s rule',
                         collections.Counter(s['today_product_class'] for s in t1)))
        out.append(table('routed slices: product today -> proposed (class)', collections.Counter(
            (s['today_product_class'], s['product_class']) for s in routed)))
        out.append(table('routed slices of junctions that are T1/T2 today: does the product change?', collections.Counter(
            (s['actual_fate'], s['today_product_class'], s['product_class'], s['proposed_file']) for s in routed
            if s['actual_fate'].startswith('SUCCESS'))))
    return '\n'.join(out)


if __name__ == '__main__':
    main_()
