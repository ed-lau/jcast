# -*- coding: utf-8 -*-

"""Proteome-wide survey of stitch anchors: how often an anchor of k residues at an exon boundary is repeated within its
own protein, and how often a k-mer from another protein occurs in it by chance.

Joining a slice to the canonical (PR 5 audit) searches the canonical for the slice's first and last k residues. A slice
starts at an exon start and ends at an exon end, so its anchors are the first k residues after an exon boundary (start
anchor) and the last k residues before one (end anchor). For each MANE Select protein (GENCODE), and each k:

    repeated     the anchor occurs more than once in the protein (ambiguous: the proposed rule resolves it by the
                 other anchor, or does not join)
    chance       an anchor taken from another, randomly paired protein (fixed seed) occurs in this protein: the
                 rate at which an anchor that is not in the canonical would still be found

Writes, into --out: anchors.tsv (k: anchors, repeated, chance trials, chance matches) and summary.txt.

Usage (from the repository root):
    python reports/stitch_audit/survey_anchors.py --gtf gencode.v50.annotation.gtf \\
        --transcripts gencode.v50.transcripts.fa.gz --out reports/stitch_audit/gencode_v50
"""

import argparse
import collections
import csv
import gzip
import random
import re
import subprocess
from pathlib import Path

from Bio import SeqIO

K_VALUES = [4, 5, 6, 7, 8, 9, 10, 12, 15, 20]
ATTR = re.compile(r'(\S+) "?([^";]*)"?;')


def occurrences(text, sub):
    n, i = 0, text.find(sub)
    while i != -1:
        n += 1
        i = text.find(sub, i + 1)
    return n


def mane_structures(gtf):
    """MANE Select transcript id -> (strand, exons in transcript order, CDS rows in transcript order)."""
    exons, cds, strand = collections.defaultdict(list), collections.defaultdict(list), {}
    grep = subprocess.Popen(['grep', '-F', 'MANE_Select', str(gtf)], stdout=subprocess.PIPE, text=True)
    for line in grep.stdout:
        f = line.rstrip('\n').split('\t')
        if f[2] not in ('exon', 'CDS'):
            continue
        a = dict(ATTR.findall(f[8]))
        t = a['transcript_id']
        strand[t] = f[6]
        (exons if f[2] == 'exon' else cds)[t].append((int(a['exon_number']), int(f[3]), int(f[4])))
    grep.wait()
    return {t: (strand[t], [e[1:] for e in sorted(exons[t])], [c[1:] for c in sorted(cds[t])]) for t in cds}


def protein_and_boundaries(seq, strand, exons, cds):
    """Protein (translated from the transcript sequence) and the CDS nucleotide offsets of internal exon boundaries."""
    cds_start = cds[0][0] if strand == '+' else cds[0][1]
    offset = 0
    for s, e in exons:
        if s <= cds_start <= e:
            offset += cds_start - s if strand == '+' else e - cds_start
            break
        offset += e - s + 1
    lengths = [e - s + 1 for s, e in cds]
    coding = seq[offset:offset + sum(lengths)]
    protein = str(coding.translate())
    boundaries, b = [], 0
    for n in lengths[:-1]:
        b += n
        boundaries.append(b)
    return protein, boundaries


def main_(gtf, transcripts, out):
    out.mkdir(parents=True, exist_ok=True)
    structures = mane_structures(gtf)
    proteins = []
    skipped = collections.Counter()
    with gzip.open(transcripts, 'rt') as fh:
        for rec in SeqIO.parse(fh, 'fasta'):
            t = rec.id.split('|')[0]
            if t not in structures:
                continue
            protein, boundaries = protein_and_boundaries(rec.seq, *structures[t])
            if not protein.startswith('M') or '*' in protein.rstrip('*'):
                skipped['does not translate cleanly (non-ATG start, internal stop or selenoprotein)'] += 1
                continue
            proteins.append((t, protein.rstrip('*'), boundaries))

    # each protein is paired with a random other protein (fixed seed); neighbours in the file are often paralogs
    order = list(range(len(proteins)))
    random.Random(0).shuffle(order)
    partner = {order[i]: order[(i + 1) % len(order)] for i in range(len(order))}

    rows = []
    for k in K_VALUES:
        counts = collections.Counter()
        anchors = []  # per protein: list of anchors at its boundaries
        for _, p, boundaries in proteins:
            mine = []
            for b in boundaries:
                first = -(-b // 3)  # first full codon after the boundary
                last = b // 3       # residues before the boundary's codon
                if first + k <= len(p):
                    mine.append(('start', p[first:first + k]))
                if last - k >= 0:
                    mine.append(('end', p[last - k:last]))
            anchors.append(mine)
            for kind, a in mine:
                counts[kind] += 1
                if occurrences(p, a) > 1:
                    counts[f'{kind} repeated'] += 1
        for i, (_, p, _) in enumerate(proteins):
            other = anchors[partner[i]]
            for _, a in other:
                counts['chance trials'] += 1
                counts['chance matches'] += a in p
        rows.append({'k': k, 'start anchors': counts['start'], 'start repeated': counts['start repeated'],
                     'end anchors': counts['end'], 'end repeated': counts['end repeated'],
                     'repeated %': round(100 * (counts['start repeated'] + counts['end repeated']) /
                                         max(1, counts['start'] + counts['end']), 3),
                     'chance trials': counts['chance trials'], 'chance matches': counts['chance matches'],
                     'chance %': round(100 * counts['chance matches'] / max(1, counts['chance trials']), 4)})
        print(rows[-1], flush=True)

    with open(out / 'anchors.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    lengths = sorted(len(p) for _, p, _ in proteins)
    text = [f'annotation: {gtf.name}; transcripts: {transcripts.name}\n',
            f'MANE Select proteins: {len(proteins)} (median length {lengths[len(lengths) // 2]} aa); '
            f'internal exon boundaries: {sum(len(b) for _, _, b in proteins)}; skipped: {dict(skipped)}\n']
    text += [f'k={r["k"]}: anchors {r["start anchors"] + r["end anchors"]}, repeated within the protein '
             f'{r["start repeated"] + r["end repeated"]} ({r["repeated %"]}%); chance match of a foreign k-mer '
             f'{r["chance matches"]}/{r["chance trials"]} ({r["chance %"]}%)\n' for r in rows]
    (out / 'summary.txt').write_text(''.join(text))
    print(''.join(text))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--gtf', type=Path, required=True)
    parser.add_argument('--transcripts', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    a = parser.parse_args()
    main_(a.gtf, a.transcripts, a.out)
