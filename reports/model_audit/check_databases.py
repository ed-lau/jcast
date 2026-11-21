# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §19): the heart search databases, as a search engine
will read them.

Per run (--db), its jcast_uniprot_merged.fasta:
    UniProt     the first entries are the --uniprot FASTA's, header and sequence, in order
    JCAST       the rest: each accession sp|<UniProt accession>-J<n>|, joined to a Swiss-Prot canonical of the FASTA
    log         the entries against the run's logfile.log ("<n> UniProt entries and <m> JCAST records")
    accessions  distinct (the first word of the header), repeated, and any that begins with a decoy prefix (DECOY_,
                rev_) a search would add
    sequences   distinct; JCAST records identical to a UniProt entry or to each other (JCAST writes neither); letters
                other than the 20 amino acids, U and X; stops (*); the shortest and longest JCAST records
With --pair <label>,<label>: the JCAST sequences of each database not in the other's (by sequence: the -J numbers
depend on the run).

Writes --out.

Usage (from the repository root):
    python reports/model_audit/check_databases.py --db <label>=<output dir> ... --uniprot <FASTA>
        [--pair <label>,<label> ...] --out <file>
"""

import argparse
import re
from pathlib import Path

from Bio import SeqIO

AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWYUX')
DECOYS = ('DECOY_', 'rev_')
JCAST_ID = re.compile(r'^sp\|([A-Z0-9]+(?:-\d+)?)-J\d+\|')


def check(run, uniprot) -> tuple:
    """ the lines for one run, and its JCAST sequences """
    entries = list(SeqIO.parse(Path(run) / 'jcast_uniprot_merged.fasta', 'fasta'))
    n = len(uniprot)
    head = entries[:n]
    same = len(head) == n and all(a.description == b.description and str(a.seq) == str(b.seq)
                                  for a, b in zip(head, uniprot))
    jcast = entries[n:]
    canonicals = {r.id.split('|')[1] for r in uniprot if '-' not in r.id.split('|')[1]}
    joined = [JCAST_ID.match(r.id) for r in jcast]
    log = re.search(r'jcast_uniprot_merged\.fasta: (\d+) UniProt entries and (\d+) JCAST records',
                    (Path(run) / 'logfile.log').read_text())
    ids = [r.id for r in entries]
    seqs = [str(r.seq) for r in jcast]
    sp_seqs = {str(r.seq) for r in uniprot}
    odd = sum(1 for s in seqs if set(s) - AMINO_ACIDS - {'*'})
    lengths = sorted(len(s) for s in seqs)
    lines = [f'entries: {len(entries):,}; UniProt: {n:,}, the --uniprot FASTA\'s in order: {same}; JCAST: '
             f'{len(jcast):,}',
             f'log: {int(log.group(1)):,} UniProt entries and {int(log.group(2)):,} JCAST records; equal: '
             f'{(int(log.group(1)), int(log.group(2))) == (n, len(jcast))}',
             f'JCAST accessions sp|<accession>-J<n>|: {sum(bool(m) for m in joined):,} of {len(jcast):,}; joined to a '
             f'Swiss-Prot canonical of the FASTA: {sum(1 for m in joined if m and m.group(1) in canonicals):,}',
             f'accessions: {len(set(ids)):,} distinct, {len(ids) - len(set(ids))} repeated; beginning with a decoy '
             f'prefix: {sum(1 for i in ids if i.startswith(DECOYS))}',
             f'JCAST sequences: {len(set(seqs)):,} distinct, {len(seqs) - len(set(seqs))} repeated; identical to a '
             f'UniProt entry: {len(set(seqs) & sp_seqs)}; with letters other than the 20 amino acids, U and X: {odd}; '
             f'with a stop: {sum(1 for s in seqs if "*" in s)}',
             f'JCAST lengths: {lengths[0]} to {lengths[-1]}, median {lengths[len(lengths) // 2]}; under 30: '
             f'{sum(1 for x in lengths if x < 30)}']
    return lines, set(seqs)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--db', action='append', required=True, help='<label>=<output dir>')
    p.add_argument('--uniprot', required=True)
    p.add_argument('--pair', action='append', default=[], help='<label>,<label>')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    uniprot = list(SeqIO.parse(a.uniprot, 'fasta'))
    lines, seqs = [f'--uniprot: {Path(a.uniprot).name}, {len(uniprot):,} entries'], {}
    for spec in a.db:
        label, run = spec.split('=', 1)
        found, seqs[label] = check(run, uniprot)
        lines += ['', f'== {label}: {run}'] + found
    for spec in a.pair:
        x, y = spec.split(',')
        lines += ['', f'== {x} and {y}, JCAST sequences: {len(seqs[x] & seqs[y]):,} in both; '
                      f'{len(seqs[x] - seqs[y]):,} only in {x}; {len(seqs[y] - seqs[x]):,} only in {y}']
    text = '\n'.join(lines) + '\n'
    print(text)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text)


if __name__ == '__main__':
    main()
