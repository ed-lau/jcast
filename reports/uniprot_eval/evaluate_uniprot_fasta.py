# -*- coding: utf-8 -*-

"""Evaluation: a user-supplied UniProt FASTA instead of per-gene UniProt API calls, and joining by coordinates.

On the chr15 test set, compares three sources of each gene's canonical protein:

    api        the UniProt API answer JCAST uses today (first reviewed entry for the Ensembl gene id), taken from the
               local cache (cache/uniprot-cache.db, filled by an earlier run)
    fasta      a Swiss-Prot FASTA with isoforms (--uniprot), mapped to the gene by gene name (GN=) or by sequence
               (the entry's sequence equals a protein translated from one of the gene's GTF transcripts)
    gtf        the GTF canonical JCAST uses with use_gtf_only

and measures:
    1. mapping: how many genes map to a Swiss-Prot entry by name and by sequence, and agree with the API answer
    2. coordinates: how many canonicals equal a GTF-translated transcript of the gene, so that their residues have
       genomic coordinates (needed to join by coordinates)
    3. renaming: running every chr15 junction with the API canonicals, how many output records are identical to a
       Swiss-Prot canonical or isoform (what the post-run R merge renames or drops)

Writes, into reports/uniprot_eval/: genes.tsv, records.tsv, summary.txt.

Usage (from the repository root):
    python reports/uniprot_eval/evaluate_uniprot_fasta.py --uniprot uniprot_human_sp_isoforms.fasta.gz
"""

import argparse
import collections
import csv
import gzip
import logging
import re
import shutil
import sqlite3
import sys
import tempfile
from io import StringIO
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
CACHE = REPO / 'cache' / 'uniprot-cache.db'

sys.path.insert(0, str(REPO))
for d in ['utr_audit', 'cterm_audit']:
    sys.path.insert(0, str(REPO / 'reports' / d))


def read_uniprot(path):
    """accession -> (gene name, is canonical, sequence) for a UniProt FASTA (optionally gzipped)."""
    from Bio import SeqIO
    opener = gzip.open if str(path).endswith('.gz') else open
    out = {}
    with opener(path, 'rt') as fh:
        for r in SeqIO.parse(fh, 'fasta'):
            acc = r.id.split('|')[1]
            gn = re.search(r' GN=(\S+)', r.description)
            out[acc] = (gn.group(1) if gn else None, '-' not in acc, str(r.seq))
    return out


def main_(uniprot):
    import audit_cterm
    import audit_utr
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    T = audit_utr.table

    work = Path(tempfile.mkdtemp(prefix='uniprot_eval_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')

    sp = read_uniprot(uniprot)
    sp_by_seq = collections.defaultdict(list)
    for acc, (_, _, s) in sp.items():
        sp_by_seq[s].append(acc)
    canonical_by_gn = collections.defaultdict(list)
    for acc, (gn, is_canonical, _) in sp.items():
        if is_canonical and gn:
            canonical_by_gn[gn].append(acc)

    # proteins translated from the GTF: sequence -> transcript names; transcript name -> gene id
    proteins = audit_cterm.annotated_proteins(gtf, genome)
    tx_gene = dict(zip(gtf.annot['transcript_name'], gtf.annot['gene_id']))
    gene_proteins = collections.defaultdict(set)
    for s, names in proteins.items():
        for n in names:
            gene_proteins[tx_gene.get(n.split(':')[0])].add(s)

    con = sqlite3.connect(CACHE)
    cached = {i.split('.')[0]: f for i, f in con.execute('SELECT id, seq FROM sequences')}

    def api(gene_id):
        f = cached.get(gene_id.split('.')[0])
        if f is None:
            return None
        r = next(SeqIO.parse(StringIO(f), 'fasta'))
        return r.id.split('|')[1], str(r.seq)

    genes = {}
    for rma in rmats.rmats_list:
        for gid, sym in zip(rma['gene_id'], rma['gene_symbol']):
            genes[gid] = sym
    coding = {g for g in genes if gene_proteins.get(g)}

    rows = []
    for gid in sorted(genes):
        sym = genes[gid]
        a = api(gid)
        by_name = canonical_by_gn.get(sym, [])
        by_seq = sorted({acc for s in gene_proteins.get(gid, ()) for acc in sp_by_seq.get(s, []) if sp[acc][1]})
        rows.append({
            'gene_id': gid, 'gene_symbol': sym, 'protein_coding_in_gtf': gid in coding,
            'api_cached': a is not None, 'api_accession': a[0] if a else '',
            'api_sequence_in_current_release': (a[1] == sp[a[0].split('-')[0]][2]) if a and a[0] in sp else '',
            'api_equals_gtf_protein': (a[1] in gene_proteins.get(gid, ())) if a else '',
            'fasta_by_name': ';'.join(by_name), 'fasta_by_sequence': ';'.join(by_seq),
            'name_agrees_with_api': (by_name == [a[0]]) if a else '',
            'sequence_agrees_with_api': (by_seq == [a[0]]) if a else '',
        })
    with open(OUT / 'genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # 3. run every junction with the API canonical (as a default run), then compare records with Swiss-Prot
    params.use_gtf_only = False

    def from_cache(self, reviewed='true'):
        f = cached.get(self.j.gene_id.split('.')[0])
        if f is None or reviewed != 'true':
            return SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})
        return next(SeqIO.parse(StringIO(f), 'fasta'))
    Sequence.get_canonical_aa_uniprot = from_cache
    out_dir = work / 'out'
    out_dir.mkdir()
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
    for rma in rmats.rmats_list:
        for i in range(len(rma)):
            main._translate_one(Junction(**rma.iloc[i].to_dict()), gtf, genome, args, pred_bound=args.read)
    records = []
    for p in sorted(out_dir.glob('psq_*.fasta')):
        for r in SeqIO.parse(p, 'fasta'):
            accs = sp_by_seq.get(str(r.seq), [])
            records.append({'file': p.name, 'header': r.id, 'length': len(r.seq),
                            'swissprot': ';'.join(accs),
                            'match': ('isoform' if any('-' in a for a in accs) else
                                      'canonical' if accs else 'none')})
    with open(OUT / 'records.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0]), delimiter='\t')
        w.writeheader()
        w.writerows(records)

    cached_rows = [r for r in rows if r['api_cached']]
    coding_rows = [r for r in rows if r['protein_coding_in_gtf']]
    out = [f'UniProt FASTA: {Path(uniprot).name}: {len(sp)} entries, '
           f'{sum(1 for v in sp.values() if v[1])} canonical; GTF: Ensembl 89 chr15\n',
           f'chr15 rMATS genes: {len(rows)}; with a protein-coding transcript in the GTF: {len(coding_rows)}; '
           f'with a cached API answer: {len(cached_rows)}\n',
           T('protein-coding genes: Swiss-Prot canonicals matching the gene name (GN=)', collections.Counter(
               min(len(r['fasta_by_name'].split(';')) if r['fasta_by_name'] else 0, 2) for r in coding_rows)),
           T('protein-coding genes: Swiss-Prot canonicals equal to one of the gene\'s GTF proteins', collections.Counter(
               min(len(r['fasta_by_sequence'].split(';')) if r['fasta_by_sequence'] else 0, 2) for r in coding_rows)),
           T('cached API answers: name mapping agrees | sequence mapping agrees', collections.Counter(
               (r['name_agrees_with_api'], r['sequence_agrees_with_api']) for r in cached_rows)),
           T('cached API answers: the cached sequence equals the current release', collections.Counter(
               r['api_sequence_in_current_release'] for r in cached_rows)),
           T('cached API answers: the canonical equals a GTF-translated transcript (coordinates available)',
             collections.Counter(r['api_equals_gtf_protein'] for r in cached_rows)),
           T('output records (API canonicals) identical to a Swiss-Prot entry: file | match', collections.Counter(
               (r['file'], r['match']) for r in records)),
           T('output records: match', collections.Counter(r['match'] for r in records))]
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--uniprot', type=Path, required=True, help='UniProt FASTA (Swiss-Prot, with isoforms)')
    main_(parser.parse_args().uniprot)
