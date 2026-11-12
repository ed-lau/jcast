# -*- coding: utf-8 -*-

""" Canonical proteins from a UniProt FASTA supplied with --uniprot, and the merged Swiss-Prot + JCAST database

See reports/2026-11-11-uniprot-fasta-evaluation.md.
"""

import collections
import csv
import gzip
import hashlib
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

from jcast.annots import gene_rows

# JCAST files whose records joined to a Swiss-Prot canonical go into the merged database
MERGE_TIERS = ['T1', 'T1_alt_cterm', 'T1_alt_nterm', 'T1_alt_both']

MERGED_FILE = 'jcast_uniprot_merged.fasta'
ISOFORMS_FILE = 'psq_isoforms.tsv'
CANONICALS_FILE = 'uniprot_canonicals.tsv'


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def gene_proteins(gtf, genome_index, gene_id) -> set:
    """
    Proteins translated from every CDS-bearing transcript of a gene in the GTF, whatever its biotype or support level.

    :param gtf: ReadAnnotations
    :param genome_index: genome in memory
    :param gene_id: gene id
    :return: set of protein sequences (without the stop)
    """
    rows = gene_rows(gtf, gene_id)
    cds = rows[rows['feature'] == 'CDS']
    proteins = set()
    for _, t in cds.groupby('transcript_id'):
        t = t.sort_values('start')
        chromosome = genome_index[str(t['seqname'].iloc[0])].seq
        nt = ''.join(str(chromosome[s - 1:e]) for s, e in zip(t['start'], t['end'])).upper()
        if t['strand'].iloc[0] == '-':
            nt, frame = str(Seq(nt).reverse_complement()), int(t['frame'].iloc[-1])
        else:
            frame = int(t['frame'].iloc[0])
        nt = nt[frame:]
        proteins.add(str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip('*'))
    return proteins


class UniProtFasta(object):

    def __init__(self, logger, path):
        """
        Read a UniProt FASTA (e.g. Swiss-Prot human with isoforms; optionally gzipped).

        :param logger: logger
        :param path: path to the FASTA
        """
        self.logger = logger
        self.path = path

        with _open(path) as fh:
            self.entries = list(SeqIO.parse(fh, 'fasta'))
        with open(path, 'rb') as fh:
            sha256 = hashlib.sha256(fh.read()).hexdigest()

        self.by_sequence = collections.defaultdict(list)     # sequence -> accessions
        self.canonicals = {'sp': collections.defaultdict(list), 'tr': collections.defaultdict(list)}  # by sequence
        self.by_gene_name = {'sp': collections.defaultdict(list), 'tr': collections.defaultdict(list)}
        self.records = {}                                     # accession -> record

        for r in self.entries:
            db, accession = r.id.split('|')[:2]
            self.records[accession] = r
            self.by_sequence[str(r.seq)].append(accession)
            # isoform accessions carry a suffix (P12345-2); canonical ones do not
            if db in self.canonicals and '-' not in accession:
                self.canonicals[db][str(r.seq)].append(accession)
                gene_name = re.search(r' GN=(\S+)', r.description)
                if gene_name:
                    self.by_gene_name[db][gene_name.group(1)].append(accession)

        n_canonical = sum(len(v) for d in self.canonicals.values() for v in d.values())
        self.logger.info(f'UniProt FASTA {os.path.basename(str(path))} (sha256 {sha256}): {len(self.entries)} entries, '
                         f'{n_canonical} canonical.')

        self.chosen = {}        # gene id -> (gene symbol, source, accession, canonical SeqRecord)

    def canonical(self, gene_id, gene_symbol, proteins, gtf_canonical):
        """
        The gene's canonical: a UniProt canonical whose sequence equals one of the gene's GTF-translated proteins,
        else the one whose gene name (GN=) is the gene symbol, in each case only if exactly one entry matches and
        reviewed (sp) before unreviewed (tr) entries; otherwise the GTF canonical. Cached per gene.

        :param gene_id: gene id
        :param gene_symbol: gene symbol
        :param proteins: callable returning the gene's GTF-translated proteins
        :param gtf_canonical: callable returning the GTF canonical SeqRecord (possibly empty)
        :return: SeqRecord
        """
        if gene_id not in self.chosen:
            self.chosen[gene_id] = self._choose(gene_id, gene_symbol, proteins, gtf_canonical)
            # tagged so that, with -p > 1, the main process logs it once per gene, as a serial run does
            self.logger.info('Canonical for {0} {1}: {2} {3}'.format(gene_id, gene_symbol, *self.chosen[gene_id][1:3]),
                             extra={'jcast_canonical_for': gene_id})
        return self.chosen[gene_id][3]

    def _choose(self, gene_id, gene_symbol, proteins, gtf_canonical):
        gene_proteins_ = proteins()
        for db in ['sp', 'tr']:
            by_sequence = sorted({a for p in gene_proteins_ for a in self.canonicals[db].get(p, [])})
            if len(by_sequence) == 1:
                return gene_symbol, f'{db} by sequence', by_sequence[0], self.records[by_sequence[0]]
            by_name = self.by_gene_name[db].get(gene_symbol, [])
            if len(by_name) == 1:
                return gene_symbol, f'{db} by gene name', by_name[0], self.records[by_name[0]]
        record = gtf_canonical()
        return gene_symbol, 'gtf' if len(record) > 0 else 'none', '', record

    def write_outputs(self, out_dir) -> None:
        """
        Write, into out_dir:
            uniprot_canonicals.tsv      the canonical chosen for each gene, and how
            psq_isoforms.tsv            every JCAST record (except psq_canonical) identical to a UniProt entry
            jcast_uniprot_merged.fasta  every UniProt entry, then the records of MERGE_TIERS joined to a Swiss-Prot
                                        canonical that are not identical to a UniProt entry or an earlier record,
                                        each renamed sp|<accession>-J<n>|<entry name>, followed by its JCAST
                                        header line

        :param out_dir: output directory
        """
        with open(os.path.join(out_dir, CANONICALS_FILE), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['gene_id', 'gene_symbol', 'source', 'accession'])
            for gene_id, (gene_symbol, source, accession, _) in sorted(self.chosen.items()):
                w.writerow([gene_id, gene_symbol, source, accession])

        files = sorted(f for f in os.listdir(out_dir)
                       if f.startswith('psq_') and f.endswith('.fasta') and f != 'psq_canonical.fasta')
        with open(os.path.join(out_dir, ISOFORMS_FILE), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['file', 'header', 'accessions', 'entry'])
            n_identical = 0
            for f in files:
                for r in SeqIO.parse(os.path.join(out_dir, f), 'fasta'):
                    accessions = self.by_sequence.get(str(r.seq))
                    if accessions:
                        n_identical += 1
                        w.writerow([f, r.id, ';'.join(accessions),
                                    'isoform' if any('-' in a for a in accessions) else 'canonical'])

        seen = set(self.by_sequence)
        merged, numbers = [], collections.Counter()
        for tier in MERGE_TIERS:
            path = os.path.join(out_dir, f'psq_{tier}.fasta')
            if not os.path.exists(path):
                continue
            for r in SeqIO.parse(path, 'fasta'):
                fields = r.id.split('|')
                if fields[0] != 'sp' or str(r.seq) in seen:
                    continue
                seen.add(str(r.seq))
                numbers[fields[1]] += 1
                r.id = f'sp|{fields[1]}-J{numbers[fields[1]]}|{fields[2]}'
                merged.append(r)
        with open(os.path.join(out_dir, MERGED_FILE), 'w') as fh:
            SeqIO.write(self.entries + merged, fh, 'fasta')

        sources = collections.Counter(c[1] for c in self.chosen.values())
        self.logger.info(f'Canonicals by source: {dict(sorted(sources.items()))}. JCAST records identical to a '
                         f'UniProt entry: {n_identical} (see {ISOFORMS_FILE}). {MERGED_FILE}: {len(self.entries)} '
                         f'UniProt entries and {len(merged)} JCAST records from {", ".join(MERGE_TIERS)}.')
