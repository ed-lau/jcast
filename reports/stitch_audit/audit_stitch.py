# -*- coding: utf-8 -*-

"""Audit of PR 5: where a translated slice is joined ("stitched") to the canonical protein, and the anchor length.

A slice is joined by finding its first and last params.aa_stitch_length (10) residues in the canonical with str.find
(Sequence.stitch_to_canonical_aa, and Sequence._join for alternative termini). REVIEW.md MED-3: the first
occurrence is taken, an anchor occurring more than once is not detected, and an end anchor found before the start
anchor is not rejected.

The worker runs every chr15 junction through main._translate_one and records every join attempt: the peptide, the
mode (both ends, N-terminus only, C-terminus only), the canonical and the product. With the canonical translated
from the GTF, every canonical residue and every slice residue has genomic coordinates, so the correct join point
is known: the canonical residue whose codon starts where the peptide's first (last) codon starts. A peptide ending at
a stop codon, or starting at the gene's start codon, is complete at that end and is not joined there. The audit
compares the products of several joining rules (see join()) with the correct product, for anchor lengths 4 to 20.

A second run takes the canonical from the local UniProt cache (cache/uniprot-cache.db, keyed by gene id), as a
default run would; it has no coordinates, so only anchor multiplicity is measured.

Writes, into reports/stitch_audit/:

    calls.tsv         every join attempt (GTF canonical): anchors found, true anchors, product classification
    sweep.tsv         anchor length x rule: correct joins, wrong joins, missed joins, unjoined
    summary.txt       the aggregate numbers quoted in the report

Usage (from the repository root):
    python reports/stitch_audit/audit_stitch.py
"""

import argparse
import collections
import csv
import json
import logging
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']
K_RANGE = range(4, 21)
UNIPROT_CACHE = REPO / 'cache' / 'uniprot-cache.db'


# ---- worker: runs the pipeline and records every join attempt ----------------------------------------------

def coding_positions(transcript, strand):
    """Genomic position of every nucleotide of an annotated transcript's exons, in transcript orientation."""
    positions = sorted(p for s, e in transcript.exons for p in range(s, e + 1))
    return positions if strand == '+' else positions[::-1]


def true_anchor(cpos_index, cphase, cpos, span):
    """Canonical residue index whose codons are at the genomic positions span (3k nt), or -1."""
    i = cpos_index.get(span[0])
    if i is None or (i - cphase) % 3 or cpos[i:i + len(span)] != span:
        return -1
    return (i - cphase) // 3


def worker(out_dir, canonical_from):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import helpers, main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    params.use_gtf_only = canonical_from == 'gtf'

    if canonical_from == 'uniprot_cache':
        con = sqlite3.connect(UNIPROT_CACHE)
        cached = {i.split('.')[0]: s for i, s in con.execute('SELECT id, seq FROM sequences')}

        def from_cache(self, reviewed='true'):
            fasta = cached.get(self.j.gene_id.split('.')[0])
            if fasta is None or reviewed != 'true':
                return SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})
            return list(SeqIO.parse(StringIO(fasta), 'fasta'))[0]
        Sequence.get_canonical_aa_uniprot = from_cache

    # translation starts: (slice nucleotides, strand) -> {peptide: start}
    starts = collections.defaultdict(dict)
    make_pep = helpers.make_pep

    def recording_make_pep(nt, strand, phase, terminate=True):
        aa = make_pep(nt, strand, phase, terminate)
        starts[(str(nt).upper(), strand)].setdefault(aa, phase)
        return aa
    helpers.make_pep = recording_make_pep

    calls = []
    current = {}

    def record(self, k, aa, mode, product):
        nt = [self.slice1_nt, self.slice2_nt][k]
        start = starts.get((str(nt.seq).upper(), self.j.strand), {}).get(aa)
        row = {'type': current['type'], 'rmats_id': self.j.name, 'gene_id': self.j.gene_id,
               'gene_symbol': self.j.gene_symbol, 'strand': self.j.strand, 'slice': k + 1, 'mode': mode,
               'aa': aa, 'canonical': str(self.canonical_aa.seq),
               'product': None if product is None else str(product.seq), 'start': start}
        tx = getattr(self, 'gtf_canonical_transcript', None)
        if canonical_from == 'gtf' and tx is not None and start is not None:
            cpos = coding_positions(tx, self.j.strand)
            cphase = tx.starting_translation_phase
            spos = self.j.slice_positions()[k][start:start + 3 * len(aa)]
            index = {p: i for i, p in enumerate(cpos)}
            # canonical residues whose codons start, in frame, at the genomic position where the peptide's first and
            # last codons start (-1: none). A codon is placed by its first base: a splice site a few bases away (e.g.
            # NAGNAG) changes its last base but not where it is joined
            row['true_first'] = true_anchor(index, cphase, cpos, spos[:1]) if aa else -1
            row['true_last'] = true_anchor(index, cphase, cpos, spos[-3:-2]) if aa else -1
            # the coordinates must reproduce the canonical: translate it back from them
            chrom = self.j.chr if self.j.chr in self.genome_index else self.j.chr[3:]
            g = self.genome_index[chrom].seq
            pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            canonical_nt = ''.join(g[p - 1].upper() if self.j.strand == '+' else pair.get(g[p - 1].upper(), 'N')
                                   for p in cpos)
            row['coords_ok'] = make_pep(canonical_nt, '+', cphase, terminate=False) == row['canonical']
        # whether the peptide ends at a stop codon (a complete C-terminus) rather than at the end of the slice
        seq = str(nt.seq).upper() if self.j.strand == '+' else helpers.get_complementary(str(nt.seq).upper())
        row['ends_at_stop'] = start is not None and seq[start + 3 * len(aa):start + 3 * len(aa) + 3] in \
            {'TAA', 'TAG', 'TGA'}
        # the gene's start codon (the most upstream annotated start; trimming starts a slice that contains it there)
        gene_start = self.j.tx0 if self.j.strand == '+' else self.j.tx1
        row['starts_at_start'] = start is not None and self.j.slice_positions()[k][start] == gene_start
        calls.append(row)

    stitch = Sequence.stitch_to_canonical_aa

    def recording_stitch(self, slice_to_stitch, slice_has_ptc=False):
        k = slice_to_stitch - 1
        before = [self.slice1_stitched, self.slice2_stitched][k]
        current['in_stitch'] = True  # since PR 5, stitch_to_canonical_aa calls _join: record the call once
        try:
            result = stitch(self, slice_to_stitch, slice_has_ptc)
        finally:
            current['in_stitch'] = False
        after = [self.slice1_stitched, self.slice2_stitched][k]
        if len(self.canonical_aa) > 0:
            record(self, k, [self.slice1_aa, self.slice2_aa][k], 'both_or_n' if slice_has_ptc else 'both',
                   after if after is not before else None)
        return result
    Sequence.stitch_to_canonical_aa = recording_stitch

    join = Sequence._join

    def recording_join(self, aa, mode, **kwargs):
        product = join(self, aa, mode, **kwargs)
        if mode != 'none' and len(self.canonical_aa) > 0 and not current.get('in_stitch'):
            k = next((i for i, nt in enumerate([self.slice1_nt, self.slice2_nt])
                      if aa in starts.get((str(nt.seq).upper(), self.j.strand), {})), 0)
            record(self, k, aa, mode, product)
        return product
    Sequence._join = recording_join

    get_canonical = Sequence.get_canonical_aa

    def keeping_genome(self, gtf, genome_index):
        self.genome_index = genome_index
        return get_canonical(self, gtf=gtf, genome_index=genome_index)
    Sequence.get_canonical_aa = keeping_genome

    out_dir = Path(out_dir)
    work = Path(tempfile.mkdtemp(prefix='stitch_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
    for jt, rma in zip(TYPES, rmats.rmats_list):
        current['type'] = jt
        for i in range(len(rma)):
            main._translate_one(Junction(**rma.iloc[i].to_dict()), gtf, genome, args, pred_bound=args.read)
    (out_dir / 'calls.json').write_text(json.dumps(calls, default=int))


def run(out_dir, canonical_from):
    out_dir.mkdir(parents=True)
    env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND='Agg')
    subprocess.run([sys.executable, __file__, '--worker', str(out_dir), canonical_from], env=env, check=True,
                   cwd=str(out_dir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.loads((out_dir / 'calls.json').read_text())


# ---- joining rules -----------------------------------------------------------------------------------------

def occurrences(text, sub):
    """Start of every (overlapping) occurrence of sub in text."""
    out, i = [], text.find(sub)
    while i != -1:
        out.append(i)
        i = text.find(sub, i + 1)
    return out


RULES = ['first', 'unique', 'consistent', 'first+ends', 'unique+ends', 'consistent+ends']


def join(canonical, aa, mode, k, rule, ends_at_stop=False, starts_at_start=False):
    """
    The product of joining aa to canonical with anchor length k. mode 'both_or_n' joins both ends if possible, else
    the N-terminus only. Rules:
        first        current code: the first occurrence of each anchor
        unique       each anchor must occur exactly once, and the start anchor must not be after the end anchor
        consistent   as unique, but an anchor occurring more than once is resolved by the pair of occurrences (start
                     not after end) whose canonical span is closest to the peptide length; a tie is not joined
        +ends        a peptide ending at a stop codon is complete at its C-terminus and one starting at the gene's
                     start codon is complete at its N-terminus: neither is joined at that end
    """
    if len(aa) < k or not canonical:
        return None
    n_done = False
    if rule.endswith('+ends'):
        rule = rule[:-len('+ends')]
        n_done = starts_at_start and mode != 'c'
        if ends_at_stop and mode != 'c':
            mode = 'n'
    o0, o1 = occurrences(canonical, aa[:k]), occurrences(canonical, aa[-k:])
    if rule == 'first':
        m0, m1 = (o0[0] if o0 else -1), (o1[0] if o1 else -1)
    else:
        m0, m1 = (o0[0] if len(o0) == 1 else -1), (o1[0] if len(o1) == 1 else -1)
    if rule == 'consistent' and mode in ('both', 'both_or_n') and not n_done and (len(o0) > 1 or len(o1) > 1):
        pairs = sorted((abs(b + k - a - len(aa)), a, b) for a in o0 for b in o1 if a <= b)
        if pairs and (len(pairs) == 1 or pairs[0][0] < pairs[1][0]):
            m0, m1 = pairs[0][1], pairs[0][2]
    if n_done:
        before = ''  # nothing to add before the peptide
    else:
        before = canonical[:m0] if m0 != -1 else None
    both_ok = before is not None and m1 != -1 and (rule == 'first' or n_done or m0 <= m1)
    if mode in ('both', 'both_or_n') and both_ok:
        return before + aa + canonical[m1 + k:]
    if mode in ('n', 'both_or_n') and before is not None:
        return before + aa
    if mode == 'c' and m1 != -1:
        return aa + canonical[m1 + k:]
    return None


def truth(c):
    """
    The correct product: the canonical joined where the peptide's end codons lie in it. A peptide ending at a stop
    codon has a complete C-terminus, and one starting at the gene's start codon a complete N-terminus; neither is
    joined at that end.
    """
    can, aa, t0, t1 = c['canonical'], c['aa'], c['true_first'], c['true_last']
    if c['starts_at_start'] and c['mode'] != 'c':
        t0, before = 0, ''
    else:
        before = can[:t0]
    if c['ends_at_stop'] and c['mode'] != 'c':
        return before + aa if t0 != -1 else None
    if c['mode'] in ('both', 'both_or_n') and t0 != -1 and t1 != -1:
        return before + aa + can[t1 + 1:]
    if c['mode'] in ('n', 'both_or_n') and t0 != -1:
        return before + aa
    if c['mode'] == 'c' and t1 != -1:
        return aa + can[t1 + 1:]
    return None


def classify(product, truth):
    if product is None:
        return 'unjoined (correct)' if truth is None else 'missed join'
    if truth is None:
        return 'wrong join (no true anchor)'
    return 'correct join' if product == truth else 'wrong join (misplaced)'


# ---- comparison --------------------------------------------------------------------------------------------

def main_():
    sys.path.insert(0, str(REPO / 'reports' / 'utr_audit'))
    import audit_utr
    T = audit_utr.table

    work = Path(tempfile.mkdtemp(prefix='stitch_audit_'))
    gtf_calls = run(work / 'gtf', 'gtf')
    uni_calls = run(work / 'uniprot_cache', 'uniprot_cache') if UNIPROT_CACHE.exists() else []

    # deduplicate: the same peptide joined the same way to the same canonical gives the same product
    seen, calls = set(), []
    for c in gtf_calls:
        key = (c['canonical'], c['aa'], c['mode'])
        if key not in seen:
            seen.add(key)
            calls.append(c)
    usable = [c for c in calls if c.get('true_first') is not None and c.get('coords_ok')]

    rows = []
    for c in usable:
        o0, o1 = occurrences(c['canonical'], c['aa'][:10]), occurrences(c['canonical'], c['aa'][-10:])
        product = join(c['canonical'], c['aa'], c['mode'], 10, 'first', c['ends_at_stop'], c['starts_at_start'])
        assert product == c['product'], 'the replayed join must reproduce the pipeline'
        t = truth(c)
        rows.append({'type': c['type'], 'rmats_id': c['rmats_id'], 'gene_symbol': c['gene_symbol'],
                     'strand': c['strand'], 'slice': c['slice'], 'mode': c['mode'], 'aa_length': len(c['aa']),
                     'canonical_length': len(c['canonical']),
                     'start_anchor_occurrences': len(o0), 'end_anchor_occurrences': len(o1),
                     'found0': o0[0] if o0 else -1, 'found1': o1[0] if o1 else -1,
                     'true_first': c['true_first'], 'true_last': c['true_last'],
                     'ends_at_stop': c['ends_at_stop'], 'starts_at_start': c['starts_at_start'],
                     'current': classify(product, t),
                     **{rule: classify(join(c['canonical'], c['aa'], c['mode'], 10, rule, c['ends_at_stop'], c['starts_at_start']), t)
                        for rule in RULES[1:]},
                     'aa': c['aa']})
    with open(OUT / 'calls.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    sweep = []
    for k in K_RANGE:
        for rule in RULES:
            counts = collections.Counter(
                classify(join(c['canonical'], c['aa'], c['mode'], k, rule, c['ends_at_stop'], c['starts_at_start']), truth(c))
                for c in usable)
            sweep.append(dict({'k': k, 'rule': rule}, **{x: counts[x] for x in
                              ['correct join', 'wrong join (misplaced)', 'wrong join (no true anchor)',
                               'missed join', 'unjoined (correct)']}))
    with open(OUT / 'sweep.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(sweep[0]), delimiter='\t')
        w.writeheader()
        w.writerows(sweep)

    def multiplicity(cs, label):
        n = collections.Counter()
        for c in cs:
            if len(c['aa']) < 10:
                n['peptide shorter than the anchor'] += 1
                continue
            o0, o1 = occurrences(c['canonical'], c['aa'][:10]), occurrences(c['canonical'], c['aa'][-10:])
            n[f'start anchor x{min(len(o0), 2)}{"+" if len(o0) > 1 else ""}'] += 1
            n[f'end anchor x{min(len(o1), 2)}{"+" if len(o1) > 1 else ""}'] += 1
            if o0 and o1 and o1[0] < o0[0]:
                n['end anchor found before start anchor'] += 1
        return T(label, n)

    uni = []
    seen = set()
    for c in uni_calls:
        key = (c['canonical'], c['aa'], c['mode'])
        if key not in seen:
            seen.add(key)
            uni.append(c)

    out = [f'junctions: chr15 test set; join attempts recorded (GTF canonical): {len(gtf_calls)}; distinct: '
           f'{len(calls)}; with genomic truth (coordinates reproduce the canonical): {len(usable)}; '
           f'coordinates failed: {sum(1 for c in calls if c.get("true_first") is not None and not c.get("coords_ok"))}; '
           f'no translation start found: {sum(1 for c in calls if c.get("start") is None)}\n',
           T('join attempts by mode (distinct)', collections.Counter(c['mode'] for c in calls)),
           T('current rule (first occurrence), k = 10', collections.Counter(r['current'] for r in rows)),
           *[T(f'rule {rule}, k = 10', collections.Counter(r[rule] for r in rows)) for rule in RULES[1:]],
           T('current -> consistent+ends, k = 10 (changed only)', collections.Counter(
               (r['current'], r['consistent+ends']) for r in rows if r['current'] != r['consistent+ends'])),
           T("join attempts by mode | ends at a stop codon | starts at the gene's start codon", collections.Counter(
               (c['mode'], c['ends_at_stop'], c['starts_at_start']) for c in usable)),
           multiplicity(usable, 'anchor multiplicity in the GTF canonical, k = 10'),
           f'UniProt-cache run: join attempts {len(uni_calls)}; distinct {len(uni)}; genes '
           f'{len({c["gene_id"] for c in uni})}\n',
           multiplicity(uni, 'anchor multiplicity in the UniProt canonical, k = 10')]
    for k in (5, 10):
        listed = collections.Counter()
        lines = [f'consistent+ends, k = {k}: wrong and missed joins']
        for c in usable:
            t = truth(c)
            got = join(c['canonical'], c['aa'], c['mode'], k, 'consistent+ends', c['ends_at_stop'], c['starts_at_start'])
            verdict = classify(got, t)
            if verdict.startswith(('wrong', 'missed')):
                aa = c['aa']
                o0, o1 = occurrences(c['canonical'], aa[:k]), occurrences(c['canonical'], aa[-k:])
                reason = ('peptide shorter than the anchor' if len(aa) < k else
                          f'start anchor x{len(o0)}, end anchor x{len(o1)}')
                lines.append(f'  {verdict:24s} {c["type"]} {c["rmats_id"]} {c["gene_symbol"]} slice {c["slice"]} '
                             f'{c["mode"]} len {len(aa)}; {reason}; true {c["true_first"]}..{c["true_last"]}; '
                             f'found {o0[:3]} {o1[:3]}')
        out.append('\n'.join(lines) + '\n')
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())
    print(open(OUT / 'sweep.tsv').read())


if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[1] == '--worker':
        worker(*sys.argv[2:])
    else:
        main_()
