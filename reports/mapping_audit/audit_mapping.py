# -*- coding: utf-8 -*-

"""Audit of the Ensembl-to-UniProt mapping under --uniprot (reports/2026-11-12-uniprot-mapping-audit.md).

--uniprot maps no identifiers. For each gene it takes the UniProt canonical whose sequence equals a protein
translated from the gene's GTF transcripts, else the entry whose GN= is the gene symbol, else the GTF canonical
(jcast/uniprot.py). The UniProt API path instead looks the gene up by UniProt's Ensembl cross-references. This
script compares the two for every gene of a whole-genome run:

    the run's choice       uniprot_canonicals.tsv and psq_canonical.fasta of a jcast -c --uniprot run
    cross-references       UniProtKB human entries with their Ensembl genes and sequences, from the UniProt REST
                           API (downloaded by --download, with the release in the response headers)
    API answers            the local API cache (cache/uniprot-cache.db), for the genes it holds
    gene types             the GTF cache (<gtf>.cached) of the run's GTF

With --gtf and --genome, it also works out what a FASTA with TrEMBL entries as well would give the protein-coding
genes left without a canonical: the run's rules (jcast/uniprot.py) applied to the cross-referenced TrEMBL entries,
and the run's output records of those genes.

Writes into --out:
    genes.tsv      one row per gene: the run's choice, the cross-referenced entries, the API's cached answer, the
                   comparison
    summary.txt    the tables quoted in the report

Usage (from the repository root):
    python reports/mapping_audit/audit_mapping.py --jcast_out <dir> --uniprot <the FASTA it used> \\
        --gtf_cached <gtf>.cached --xrefs <file>.tsv.gz [--download] [--api_cache cache/uniprot-cache.db] \\
        [--gtf <gtf> --genome <fa>]
"""

import argparse
import collections
import csv
import gzip
import io
import logging
import re
import sqlite3
import urllib.request
from pathlib import Path

from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
XREF_URL = ('https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv'
            '&fields=accession,reviewed,gene_primary,xref_ensembl_full,sequence'
            '&query=(organism_id:9606)+AND+(database:ensembl)')


def unversioned(gene_id):
    return gene_id.split('.')[0]


def download(path):
    with urllib.request.urlopen(XREF_URL) as r:
        release = r.headers.get('X-UniProt-Release')
        Path(path).write_bytes(r.read())
    Path(str(path) + '.release').write_text(f'{release}\n{XREF_URL}\n')


def read_xrefs(path):
    """ gene -> {'sp': accessions, 'tr': accessions}; accession -> (sequence, primary gene name) """
    genes = collections.defaultdict(lambda: {'sp': set(), 'tr': set()})
    entries = {}
    with gzip.open(path, 'rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            db = 'sp' if row['Reviewed'] == 'reviewed' else 'tr'
            entries[row['Entry']] = (row['Sequence'], row['Gene Names (primary)'], db)
            for gene in set(re.findall(r'ENSG\d+', row['Ensembl'])):
                genes[gene][db].add(row['Entry'])
    return genes, entries


def read_gene_types(path):
    import pandas as pd
    types = {}
    for chunk in pd.read_table(path, usecols=['feature', 'gene_id', 'gene_type'], chunksize=2_000_000):
        g = chunk[chunk.feature == 'gene']
        types.update(zip(g.gene_id.map(unversioned), g.gene_type))
    return types


def read_api_cache(path):
    """ gene -> (db, accession, sequence) of the API's cached answer """
    out = {}
    if not Path(path).exists():
        return out
    con = sqlite3.connect(path)
    for gene, fasta in con.execute('SELECT id, seq FROM sequences'):
        r = next(SeqIO.parse(io.StringIO(fasta), 'fasta'))
        db, accession = r.id.split('|')[:2]
        out[unversioned(gene)] = (db, accession, str(r.seq))
    con.close()
    return out


def compare(source, accession, sequence, xref, entries):
    """ how the run's choice relates to the gene's cross-referenced entries """
    sp, tr = xref['sp'], xref['tr']
    if source.startswith('sp'):
        if accession in sp:
            return 'same entry'
        if sp and sequence in {entries[a][0] for a in sp}:
            return 'other entry, same sequence'
        if sp:
            return 'other entry, other sequence'
        return 'Swiss-Prot entry; cross-references: TrEMBL only' if tr else 'Swiss-Prot entry; no cross-reference'
    if sp:
        return f'{source}; cross-references: Swiss-Prot'
    if tr:
        same = sequence and sequence in {entries[a][0] for a in tr}
        return f'{source}; cross-references: TrEMBL only' + (' (same sequence)' if same else '')
    return f'{source}; no cross-reference'


def trembl_option(a, rows, entries, summary):
    """ the run's rules applied to TrEMBL entries, for protein-coding genes the run left without a canonical """
    import sys
    sys.path.insert(0, str(REPO))
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.uniprot import gene_proteins
    log = logging.getLogger('jcast-audit')
    gtf, genome = ReadAnnotations(log, a.gtf), ReadGenome(log, a.genome)
    by_sequence, by_name = collections.defaultdict(list), collections.defaultdict(list)
    for accession, (sequence, name, db) in entries.items():
        if db == 'tr':
            by_sequence[sequence].append(accession)
            by_name[name].append(accession)
    records = output_records(a.jcast_out, {r['gene_id'] for r in rows})
    found = collections.Counter()
    for r in rows:
        if r['gene_type'] != 'protein_coding' or r['source'] != 'none':
            continue
        matches = sorted({acc for p in gene_proteins(gtf, genome.genome, r['gene_id']) for acc in by_sequence.get(p, [])})
        if len(matches) == 1:
            r['trembl_option'] = f'tr by sequence {matches[0]}'
        elif len(by_name.get(r['gene_symbol'], [])) == 1:
            r['trembl_option'] = f'tr by gene name {by_name[r["gene_symbol"]][0]}'
        else:
            r['trembl_option'] = 'none'
        key = r['trembl_option'].rsplit(' ', 1)[0]
        found[key, 'genes'] += 1
        found[key, 'orphan records'] += records[r['gene_id']]['orphan']
    summary.append('protein-coding genes without a canonical in the run, with TrEMBL entries in the FASTA as well: ' +
                   ', '.join(f'{k[0]}: {n} {k[1]}' for k, n in sorted(found.items())))


def output_records(out_dir, gene_ids):
    """
    gene id -> {'joined': n, 'orphan': n} records in the run's psq files, except psq_canonical. A record's gene is the
    first of its header fields that is one of the run's gene ids (a gene without a name has its id as its symbol).
    """
    records = collections.defaultdict(collections.Counter)
    for f in sorted(Path(out_dir).glob('psq_*.fasta')):
        if f.name != 'psq_canonical.fasta':
            for r in SeqIO.parse(f, 'fasta'):
                gene = next((x for x in r.id.split('|') if x in gene_ids), None)
                if gene:
                    records[gene]['orphan' if 'orphan' in f.name else 'joined'] += 1
    return records


def main_(a):
    if a.download or not Path(a.xrefs).exists():
        download(a.xrefs)
    genes, entries = read_xrefs(a.xrefs)
    fasta = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(gzip.open(a.uniprot, 'rt')
                                                                if a.uniprot.endswith('.gz') else a.uniprot, 'fasta')}
    gtf_canonicals = {r.id.split('|')[0]: str(r.seq)
                      for r in SeqIO.parse(Path(a.jcast_out) / 'psq_canonical.fasta', 'fasta')
                      if r.id.startswith('ENSG')}
    types = read_gene_types(a.gtf_cached)
    api = read_api_cache(a.api_cache)

    rows = []
    for r in csv.DictReader(open(Path(a.jcast_out) / 'uniprot_canonicals.tsv'), delimiter='\t'):
        gene = unversioned(r['gene_id'])
        sequence = fasta.get(r['accession'], '') if r['accession'] else gtf_canonicals.get(r['gene_id'], '')
        xref = genes.get(gene, {'sp': set(), 'tr': set()})
        row = dict(gene_id=r['gene_id'], gene_symbol=r['gene_symbol'], gene_type=types.get(gene, '?'),
                   source=r['source'], accession=r['accession'],
                   xref_sp=','.join(sorted(xref['sp'])), xref_tr_count=len(xref['tr']),
                   comparison=compare(r['source'], r['accession'], sequence, xref, entries))
        if gene in api:
            db, accession, seq = api[gene]
            row['api_accession'] = f'{db}|{accession}'
            row['api_agrees'] = ('same entry' if accession == r['accession'] else
                                 'same sequence' if seq == sequence else 'differs')
        rows.append(row)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    s = []
    if a.gtf and a.genome:
        trembl_option(a, rows, entries, s)
    records = output_records(a.jcast_out, {r['gene_id'] for r in rows})
    for r in rows:
        r.update(joined_records=records[r['gene_id']]['joined'], orphan_records=records[r['gene_id']]['orphan'])
    cols = ['gene_id', 'gene_symbol', 'gene_type', 'source', 'accession', 'comparison', 'xref_sp', 'xref_tr_count',
            'api_accession', 'api_agrees', 'joined_records', 'orphan_records', 'trembl_option']
    (out / 'genes.tsv').write_text('\t'.join(cols) + '\n' + ''.join(
        '\t'.join(str(r.get(c, '')) for c in cols) + '\n' for r in rows))

    release = Path(str(a.xrefs) + '.release')
    s = s + [f'cross-references: {a.xrefs} (UniProt release {release.read_text().split()[0] if release.exists() else "?"})'
         f'; {len(entries)} human entries with an Ensembl gene, {len(genes)} Ensembl genes',
         f'genes in the run: {len(rows)}', '']
    for coding in [True, False]:
        sub = [r for r in rows if (r['gene_type'] == 'protein_coding') == coding]
        s.append(f'{"protein-coding" if coding else "other"} genes ({len(sub)}): the run\'s choice against UniProt\'s '
                 f'cross-references')
        groups = collections.defaultdict(list)
        for r in sub:
            groups[r['source'], r['comparison']].append(r)
        for (source, comparison), g in sorted(groups.items(), key=lambda x: -len(x[1])):
            s.append(f'  {len(g):6d}  {source:18s} {comparison}; records: '
                     f'{sum(r["joined_records"] for r in g)} joined, {sum(r["orphan_records"] for r in g)} orphan')
        s.append('')
    multi = sum(1 for r in rows if r['xref_sp'].count(',') >= 1)
    s.append(f'genes cross-referenced to more than one Swiss-Prot entry: {multi}')
    cached = [r for r in rows if 'api_agrees' in r]
    s.append(f'API answers in the local cache: {len(cached)} genes; the run\'s choice: '
             f'{dict(collections.Counter(r["api_agrees"] for r in cached))}')
    s.append('  differing, by the run\'s source and the API\'s database: ' + str(dict(collections.Counter(
        (r['source'], r['api_accession'].split('|')[0]) for r in cached if r['api_agrees'] == 'differs'))))
    (out / 'summary.txt').write_text('\n'.join(s) + '\n')
    print('\n'.join(s))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--jcast_out', required=True)
    p.add_argument('--uniprot', required=True)
    p.add_argument('--gtf_cached', required=True)
    p.add_argument('--xrefs', required=True, help='UniProt cross-reference TSV (.tsv.gz); downloaded if missing')
    p.add_argument('--download', action='store_true', help='download the cross-references even if the file exists')
    p.add_argument('--api_cache', default=str(REPO / 'cache' / 'uniprot-cache.db'))
    p.add_argument('--gtf', help='the run\'s GTF (its .cached file is read), for the TrEMBL option')
    p.add_argument('--genome', help='the run\'s genome, for the TrEMBL option')
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'genome'))
    main_(p.parse_args())
