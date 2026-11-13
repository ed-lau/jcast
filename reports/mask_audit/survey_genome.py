# -*- coding: utf-8 -*-

"""Survey of soft-masked (lowercase) and ambiguous (non-ACGT) bases in a genome FASTA, and of how many annotated
coding sequences contain them.

For every sequence in the FASTA: length and counts of uppercase A/C/G/T, N, other IUPAC codes, and lowercase bases.
Every non-ACGTN base is listed with its position. Every CDS in the GTF is intersected with the non-ACGT bases
(N runs and IUPAC sites) and with lowercase runs. Writes, into --out:

    composition.tsv     per sequence: length, ACGT, N, IUPAC (non-N), lowercase, other
    iupac_sites.tsv     every non-ACGTN base: sequence, 1-based position, base, CDS transcripts overlapping it
    cds_affected.tsv    every CDS transcript with a non-ACGT or lowercase base in its CDS
    summary.txt         the aggregate numbers quoted in the report

Usage (from the repository root):
    python reports/mask_audit/survey_genome.py --genome GENOME.fa --gtf ANNOTATION.gtf --out DIR
"""

import argparse
import bisect
import collections
import csv
import re
from pathlib import Path

ACGT = frozenset(b'ACGT')
IUPAC = 'RYKMSWBDHV'
NON_ACGT_UPPER = re.compile(rb'[^ACGT]+')
LOWER = re.compile(rb'[a-z]+')
ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_fasta(path):
    """Yield (name, sequence bytes) one sequence at a time."""
    name, chunks = None, []
    with open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.rstrip(b'\r\n'))
    if name is not None:
        yield name, b''.join(chunks)


def runs(pattern, seq):
    """Sorted, non-overlapping 1-based closed intervals (start, end) matching pattern."""
    return [(m.start() + 1, m.end()) for m in pattern.finditer(seq)]


def overlap(intervals, starts, s, e):
    """Bases of intervals (sorted, non-overlapping; starts = their start positions) inside [s, e]."""
    n, i = 0, bisect.bisect_right(starts, e) - 1
    while i >= 0 and intervals[i][1] >= s:
        n += min(e, intervals[i][1]) - max(s, intervals[i][0]) + 1
        i -= 1
    return n


def read_cds(gtf):
    """CDS rows: (seqname, start, end, transcript attributes)."""
    rows = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'CDS':
                continue
            attrs = collections.defaultdict(list)
            for k, v in ATTR.findall(f[8]):
                attrs[k].append(v)
            rows.append((f[0], int(f[3]), int(f[4]), {
                'transcript_id': attrs['transcript_id'][0], 'transcript_name': attrs['transcript_name'][0],
                'transcript_type': attrs['transcript_type'][0], 'gene_name': attrs['gene_name'][0],
                'gene_type': attrs['gene_type'][0], 'tsl': (attrs.get('transcript_support_level') or ['NA'])[0],
                'mane_select': 'MANE_Select' in attrs['tag'], 'ensembl_canonical': 'Ensembl_canonical' in attrs['tag'],
                'strand': f[6]}))
    return rows


def main_(genome, gtf, out):
    out.mkdir(parents=True, exist_ok=True)
    cds = read_cds(gtf)
    cds_by_seq = collections.defaultdict(list)
    for r in cds:
        cds_by_seq[r[0]].append(r)

    comp, sites = [], []
    per_tx = {}  # transcript id -> attributes and per-class base counts within its CDS
    for name, seq in read_fasta(genome):
        counts = {c: seq.count(bytes([c])) for c in range(65, 123) if chr(c).isalpha()}
        upper = seq.upper()
        ucounts = {c: counts[c] + counts.get(c + 32, 0) for c in range(65, 91)}
        iupac = sum(ucounts[ord(c)] for c in IUPAC)
        comp.append({'sequence': name, 'length': len(seq),
                     'ACGT': sum(counts[c] for c in ACGT), 'N': ucounts[ord('N')], 'IUPAC_non_N': iupac,
                     'lowercase': sum(v for k, v in counts.items() if 97 <= k <= 122),
                     'other': len(seq) - sum(counts.values())})

        non_acgt = runs(NON_ACGT_UPPER, upper)
        lower = runs(LOWER, seq)
        iupac_pos = [(m.start() + 1, chr(upper[m.start()])) for m in re.finditer(rb'[^ACGTN]', upper)]
        na_starts, lo_starts = [s for s, _ in non_acgt], [s for s, _ in lower]
        ip = [p for p, _ in iupac_pos]

        hits_at = collections.defaultdict(set)
        for _, s, e, a in cds_by_seq.get(name, []):
            n_non_acgt, n_lower = overlap(non_acgt, na_starts, s, e), overlap(lower, lo_starts, s, e)
            lo, hi = bisect.bisect_left(ip, s), bisect.bisect_right(ip, e)
            for p in ip[lo:hi]:
                hits_at[p].add(f"{a['transcript_name']}:{a['transcript_type']}")
            t = per_tx.setdefault(a['transcript_id'], dict(a, seqname=name, cds_nt=0, non_acgt=0, iupac_non_n=0,
                                                           lowercase=0))
            t['cds_nt'] += e - s + 1
            t['non_acgt'] += n_non_acgt
            t['iupac_non_n'] += hi - lo
            t['lowercase'] += n_lower
        for p, b in iupac_pos:
            sites.append({'sequence': name, 'position': p, 'base': b, 'cds_transcripts': ';'.join(sorted(hits_at[p]))})
        print(f'{name}: {len(seq):,} nt, N {ucounts[ord("N")]:,}, IUPAC {iupac}, lowercase {comp[-1]["lowercase"]:,}',
              flush=True)

    def write(path, rows, fields=None):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)

    write(out / 'composition.tsv', comp)
    write(out / 'iupac_sites.tsv', sites, ['sequence', 'position', 'base', 'cds_transcripts'])
    affected = sorted((t for t in per_tx.values() if t['non_acgt'] or t['lowercase']),
                      key=lambda t: (t['seqname'], t['transcript_id']))
    tx_fields = ['seqname', 'transcript_id', 'transcript_name', 'transcript_type', 'gene_name', 'gene_type', 'strand',
                 'tsl', 'mane_select', 'ensembl_canonical', 'cds_nt', 'non_acgt', 'iupac_non_n', 'lowercase']
    write(out / 'cds_affected.tsv', affected, tx_fields)

    total = {k: sum(c[k] for c in comp) for k in ['length', 'ACGT', 'N', 'IUPAC_non_N', 'lowercase', 'other']}
    pc = [t for t in per_tx.values() if t['transcript_type'] == 'protein_coding']

    def describe(label, txs):
        genes = lambda ts: len({t['gene_name'] for t in ts})
        bad = [t for t in txs if t['non_acgt']]
        iu = [t for t in txs if t['iupac_non_n']]
        lo = [t for t in txs if t['lowercase']]
        return (f'{label}: {len(txs):,} transcripts ({genes(txs):,} genes); CDS with any non-ACGT base: {len(bad)} '
                f'({genes(bad)} genes), of which with a non-N IUPAC base: {len(iu)} ({genes(iu)} genes); '
                f'CDS with any lowercase base: {len(lo):,} ({genes(lo):,} genes)\n')

    iupac_by_base = collections.Counter(s['base'] for s in sites)
    out_txt = [f'genome: {genome.name}; annotation: {gtf.name}\n',
               f'sequences: {len(comp)}; total {total["length"]:,} nt\n',
               f'uppercase ACGT {total["ACGT"]:,}; N {total["N"]:,} ({total["N"] / total["length"]:.2%}); '
               f'non-N IUPAC {total["IUPAC_non_N"]}; lowercase {total["lowercase"]:,}; other {total["other"]}\n',
               f'non-N IUPAC by base: {dict(sorted(iupac_by_base.items()))}\n',
               f'non-N IUPAC sites inside any CDS: {sum(1 for s in sites if s["cds_transcripts"])}\n',
               describe('all CDS transcripts', list(per_tx.values())),
               describe('protein_coding transcripts', pc),
               describe('MANE Select', [t for t in pc if t['mane_select']]),
               describe('protein_coding, TSL 1-2', [t for t in pc if t['tsl'] in ('1', '2')])]
    (out / 'summary.txt').write_text(''.join(out_txt))
    print(''.join(out_txt))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--gtf', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    a = parser.parse_args()
    main_(a.genome, a.gtf, a.out)
