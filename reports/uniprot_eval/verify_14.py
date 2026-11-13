# -*- coding: utf-8 -*-

"""Verification of PR 14 (--uniprot; reports/2026-11-11-uniprot-fasta-evaluation.md §8).

1. Canonicals: the canonical chosen for each chr15 rMATS gene equals the evaluation's "sequence, then name" mapping
   (genes.tsv), and is compared with the cached API answer.
2. Every junction is run through main._translate_one twice: with the cached UniProt API answers (as the evaluation,
   a default run) and with --uniprot. The UniProt API is patched to raise in the --uniprot run, so it makes no call.
   Every fate change and every FASTA record present in one run only is attributed to its gene's canonical: changed
   (or absent from the cache) vs the same.
3. The merged database equals the R merge (combine_JCAST_comprehensive_canonical.R) applied to the same psq_T1.fasta,
   for the records from psq_T1: drop tr| records, keep those not identical to a UniProt entry, unique sequences,
   named sp|<acc>-J<n>|<entry>. Compared by sequence, accession and entry name: the R script numbers every record
   before dropping those identical to a UniProt entry, so its -J numbers have gaps; JCAST numbers only the records
   it keeps.

Writes final/summary.txt, final/fates.tsv, final/fasta_diff.tsv.

Usage (from the repository root):
    python reports/uniprot_eval/verify_14.py --uniprot <Swiss-Prot human FASTA with isoforms>
"""

import argparse
import collections
import csv
import logging
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
for d in ['utr_audit']:
    sys.path.insert(0, str(REPO / 'reports' / d))


def main_(uniprot_path):
    import audit_utr
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence
    from jcast.uniprot import UniProtFasta, MERGED_FILE, CANONICALS_FILE

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    T = audit_utr.table

    work = Path(tempfile.mkdtemp(prefix='verify_14_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    params.use_gtf_only = False

    con = sqlite3.connect(CACHE)
    cached = {i.split('.')[0]: f for i, f in con.execute('SELECT id, seq FROM sequences')}

    def from_cache(self, reviewed='true'):
        f = cached.get(self.j.gene_id.split('.')[0])
        if f is None or reviewed != 'true':
            return SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})
        return next(SeqIO.parse(StringIO(f), 'fasta'))

    def no_api(self, reviewed='true'):
        raise RuntimeError('UniProt API called with --uniprot')

    def run(out_dir, uniprot):
        out_dir.mkdir()
        args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
        fates, genes = {}, {}
        for splice_type, rma in zip(['MXE', 'SE', 'RI', 'A5SS', 'A3SS'], rmats.rmats_list):
            for i in range(len(rma)):
                jx = Junction(**rma.iloc[i].to_dict())
                fates[(splice_type, jx.name)] = main._translate_one(jx, gtf, genome, args, pred_bound=1,
                                                                    uniprot=uniprot).split('.')[0]
                genes[(splice_type, jx.name)] = (jx.gene_id, jx.gene_symbol)
        if uniprot is not None:
            uniprot.write_outputs(out_dir)
        return fates, genes

    Sequence.get_canonical_aa_uniprot = from_cache
    api_fates, genes = run(work / 'api', None)
    Sequence.get_canonical_aa_uniprot = no_api
    uniprot = UniProtFasta(logger=logger, path=uniprot_path)
    fasta_fates, _ = run(work / 'fasta', uniprot)

    # 1. canonicals
    evaluation = {r['gene_id']: r for r in csv.DictReader(open(OUT / 'genes.tsv'), delimiter='\t')}
    chosen = {r['gene_id']: r for r in csv.DictReader(open(work / 'fasta' / CANONICALS_FILE), delimiter='\t')}

    def expected(e):
        for column in ['fasta_by_sequence', 'fasta_by_name']:
            if e[column] and ';' not in e[column]:
                return e[column]
        return ''
    mapping_mismatch = [g for g in chosen if chosen[g]['accession'] != expected(evaluation[g])]
    api_accession = {g: next(SeqIO.parse(StringIO(f), 'fasta')).id.split('|')[1] for g, f in cached.items()}
    api_accession = {g: api_accession.get(g.split('.')[0]) for g in chosen if g.split('.')[0] in api_accession}
    same_canonical = {g for g in chosen if chosen[g]['accession'] and chosen[g]['accession'] == api_accession.get(g)}

    # 2. fates and records
    final = OUT / 'final'
    final.mkdir(exist_ok=True)
    changed = sorted(k for k in api_fates if api_fates[k] != fasta_fates[k])
    with open(final / 'fates.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['type', 'rmats_id', 'gene_symbol', 'canonical', 'api_cache', 'uniprot_fasta'])
        for k in changed:
            g, sym = genes[k]
            w.writerow([*k, sym, 'same' if g in same_canonical else 'changed', api_fates[k], fasta_fates[k]])

    fa, ff = audit_utr.read_fastas(work / 'api'), audit_utr.read_fastas(work / 'fasta')
    diff = []
    for name in sorted(set(fa) | set(ff)):
        a, f = fa.get(name, {}), ff.get(name, {})
        for change, seqs, other in [('api only', a, f), ('uniprot only', f, a)]:
            for seq in sorted(set(seqs) - set(other)):
                gene = next((g for g in chosen if g in seqs[seq]), '')
                diff.append({'file': name, 'change': change,
                             'canonical': 'same' if gene in same_canonical else 'changed',
                             'header': seqs[seq], 'length': len(seq)})
    with open(final / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'change', 'canonical', 'header', 'length'], delimiter='\t')
        w.writeheader()
        w.writerows(diff)

    # 3. the R merge on the same psq_T1
    sp = list(uniprot.entries)
    sp_seqs = {str(r.seq) for r in sp}
    r_merged = {}
    for r in SeqIO.parse(work / 'fasta' / 'psq_T1.fasta', 'fasta'):
        if r.id.startswith('tr|'):
            continue
        if str(r.seq) not in sp_seqs and str(r.seq) not in r_merged:
            r_merged[str(r.seq)] = tuple(r.id.split('|')[1:3])
    t1_headers = {r.id for r in SeqIO.parse(work / 'fasta' / 'psq_T1.fasta', 'fasta')}
    ours = {str(r.seq): (r.id.split('|')[1].split('-J')[0], r.id.split('|')[2])
            for r in SeqIO.parse(work / 'fasta' / MERGED_FILE, 'fasta')
            if '-J' in r.id and r.description.split(' ', 1)[1].split(' ')[0] in t1_headers}
    merged_all = list(SeqIO.parse(work / 'fasta' / MERGED_FILE, 'fasta'))

    out = [f'UniProt FASTA: {Path(uniprot_path).name}; chr15 test set; API answers from the local cache\n',
           f'1. genes: {len(chosen)}; canonical differing from the evaluation\'s sequence-then-name mapping: '
           f'{len(mapping_mismatch)} {mapping_mismatch}\n',
           T('canonical source', collections.Counter(r['source'] for r in chosen.values())),
           T('canonical source | gene has a cached API answer | same accession as the API',
             collections.Counter((r['source'], g in api_accession, g in same_canonical) for g, r in chosen.items())),
           f'2. UniProt API calls in the --uniprot run: 0 (patched to raise)\n',
           f'junctions: {len(api_fates)}; fate changed: {len(changed)}\n',
           T('fate changes (canonical | API cache | --uniprot)', collections.Counter(
               ('same' if genes[k][0] in same_canonical else 'changed', api_fates[k], fasta_fates[k])
               for k in changed)),
           T('FASTA records (file | API cache | --uniprot)', collections.Counter(
               {(n, len(fa.get(n, {})), len(ff.get(n, {}))): 1 for n in sorted(set(fa) | set(ff))})),
           T('records in one run only (file | change | canonical)', collections.Counter(
               (d['file'], d['change'], d['canonical']) for d in diff)),
           f'3. merged database: {len(merged_all)} records = {len(sp)} UniProt entries + '
           f'{len(merged_all) - len(sp)} JCAST records\n',
           f'records from psq_T1: R merge {len(r_merged)}, JCAST {len(ours)}; identical sequences, accessions and entry names: '
           f'{r_merged == ours}\n']
    (final / 'summary.txt').write_text('\n'.join(out))
    print((final / 'summary.txt').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--uniprot', type=Path, required=True, help='UniProt FASTA (Swiss-Prot, with isoforms)')
    main_(parser.parse_args().uniprot)
