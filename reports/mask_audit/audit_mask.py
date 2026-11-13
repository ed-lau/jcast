# -*- coding: utf-8 -*-

"""Audit of PR 3: soft-masked (lowercase) and ambiguous bases (REVIEW.md HIGH-2/HIGH-3).

Runs every chr15 junction through main._translate_one with the unmasked test genome (Ensembl 89 `dna`), with the
soft-masked genome of the same release (`dna_sm`, identical apart from case) and with a hard-masked genome derived from
it (every soft-masked base replaced by N). A real run stops at the first
exception; here each junction's exception is caught and recorded, so every affected junction is counted.

    audit:   both genomes at BASELINE_COMMIT (or --baseline); the soft-masked genome also without a GTF canonical
             (as in a default run, which takes the canonical from UniProt), so that only the slices are translated
    verify:  the unmasked and soft-masked genomes in the working tree must each reproduce the baseline unmasked run
             exactly (fates and FASTA files); the hard-masked genome in the working tree must raise nothing

Also counts, from the chr15 GTF, how many protein-coding CDS contain soft-masked bases. Writes, into --out:

    junctions.tsv     per junction: fate (or exception) for every run
    summary.txt       the aggregate numbers quoted in the report

Usage (from the repository root):
    python reports/mask_audit/audit_mask.py --softmasked Homo_sapiens.GRCh38.dna_sm.chromosome.15.fa.gz [--verify]

The soft-masked genome is at
https://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.15.fa.gz
"""

import argparse
import collections
import csv
import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
UNMASKED = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
BASELINE_COMMIT = 'a784030'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']


# ---- worker: runs inside one code version ------------------------------------------------------------------

def worker(out_dir, genome_path, canonical_from='gtf'):
    sys.path.insert(0, str(REPO / 'reports' / 'utr_audit'))
    from audit_utr import short_fate
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    # 'gtf': canonical translated from the genome; 'uniprot': as a default run, which translates no canonical
    # (UniProt is stubbed to return nothing, so only the slices are translated)
    params.use_gtf_only = canonical_from == 'gtf'
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    out_dir = Path(out_dir)
    work = Path(tempfile.mkdtemp(prefix='mask_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=genome_path)
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask='NnXx', qvalue=[0, 1])

    rows = []
    for jt, rma in zip(TYPES, rmats.rmats_list):
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            try:
                fate, where = short_fate(main._translate_one(jx, gtf, genome, args, pred_bound=args.read)), ''
            except Exception as e:
                frames = [f for f in traceback.extract_tb(e.__traceback__) if '/jcast/' in f.filename]
                fate = f'EXCEPTION {type(e).__name__}: {e}'
                where = ' <- '.join(f'{Path(f.filename).name}:{f.name}' for f in reversed(frames[-3:]))
            rows.append({'type': jt, 'rmats_id': jx.name, 'gene_symbol': jx.gene_symbol, 'strand': jx.strand,
                         'fate': fate, 'where': where})
    (out_dir / 'junctions.json').write_text(json.dumps(rows, default=int))


def run(repo, out_dir, genome_path, canonical_from='gtf'):
    out_dir.mkdir(parents=True)
    env = dict(os.environ, PYTHONPATH=str(repo), MPLBACKEND='Agg')
    subprocess.run([sys.executable, __file__, '--worker', str(out_dir), str(genome_path), canonical_from], env=env, check=True,
                   cwd=str(out_dir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return {(r['type'], r['rmats_id']): r for r in json.loads((out_dir / 'junctions.json').read_text())}


def hardmask(softmasked, path):
    """Write the soft-masked genome with every lowercase base replaced by N (a hard-masked genome)."""
    from Bio import SeqIO
    from Bio.Seq import Seq

    with gzip.open(softmasked, 'rt') as fh:
        rec = next(SeqIO.parse(fh, 'fasta'))
    seq = str(rec.seq)
    rec.seq = Seq(seq.translate(str.maketrans('acgtn', 'NNNNN')))
    with gzip.open(path, 'wt') as fh:
        SeqIO.write(rec, fh, 'fasta')


# ---- soft-masked CDS on chr15 ------------------------------------------------------------------------------

def softmasked_cds(softmasked):
    """Protein-coding transcripts whose CDS (incl. stop codon) contains soft-masked bases, from the chr15 GTF."""
    from Bio import SeqIO
    import pandas as pd
    from gtfparse import read_gtf

    with gzip.open(softmasked, 'rt') as fh:
        seq = str(next(SeqIO.parse(fh, 'fasta')).seq)
    work = Path(tempfile.mkdtemp(prefix='mask_audit_gtf_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    g = read_gtf(str(gtf_copy))
    g = g[(g['feature'].isin(['CDS', 'stop_codon'])) & (g['transcript_biotype'] == 'protein_coding')]
    lower = collections.defaultdict(int)
    total = collections.defaultdict(int)
    for r in g.itertuples():
        s = seq[r.start - 1:r.end]
        total[r.transcript_id] += len(s)
        lower[r.transcript_id] += sum(c.islower() for c in s)
    tx = g.drop_duplicates('transcript_id').set_index('transcript_id')
    tsl = tx['transcript_support_level'].astype(str).str.split(' ').str[0]
    out = {}
    for label, ids in [('protein_coding', list(tx.index)), ('protein_coding, TSL 1-2', list(tx.index[tsl.isin(['1', '2'])]))]:
        hit = [t for t in ids if lower[t]]
        out[label] = (len(ids), tx.loc[ids, 'gene_id'].nunique(), len(hit), tx.loc[hit, 'gene_id'].nunique(),
                      sum(lower[t] for t in ids), sum(total[t] for t in ids))
    return out


# ---- comparison --------------------------------------------------------------------------------------------

def main_(softmasked, baseline_commit, verify):
    sys.path.insert(0, str(REPO / 'reports' / 'utr_audit'))
    import audit_utr

    work = Path(tempfile.mkdtemp(prefix='mask_audit_'))
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), baseline_commit],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hardmasked = work / 'hardmasked.chromosome.15.fa.gz'
    hardmask(softmasked, hardmasked)
    runs = {'baseline_unmasked': (baseline_repo, UNMASKED, 'gtf'), 'baseline_softmasked': (baseline_repo, softmasked, 'gtf'),
            'baseline_unmasked_no_gtf_canonical': (baseline_repo, UNMASKED, 'uniprot'),
            'baseline_softmasked_no_gtf_canonical': (baseline_repo, softmasked, 'uniprot'),
            'baseline_hardmasked': (baseline_repo, hardmasked, 'gtf'),
            'baseline_hardmasked_no_gtf_canonical': (baseline_repo, hardmasked, 'uniprot')}
    if verify:
        runs.update({'candidate_unmasked': (REPO, UNMASKED, 'gtf'), 'candidate_softmasked': (REPO, softmasked, 'gtf'),
                     'candidate_hardmasked': (REPO, hardmasked, 'gtf'),
                     'candidate_hardmasked_no_gtf_canonical': (REPO, hardmasked, 'uniprot')})
    try:
        res = {k: run(repo, work / k, g, c) for k, (repo, g, c) in runs.items()}
    finally:
        subprocess.run(['git', '-C', str(REPO), 'worktree', 'remove', '--force', str(baseline_repo)], check=True)

    ref = res['baseline_unmasked']
    rows = [dict({'type': k[0], 'rmats_id': k[1], 'gene_symbol': ref[k]['gene_symbol'], 'strand': ref[k]['strand']},
                 **{name: r[k]['fate'] for name, r in res.items()},
                 **{f'{name}_where': r[k]['where'] for name, r in res.items() if any(r[j]['where'] for j in r)})
            for k in ref]
    out_dir = OUT / ('final' if verify else 'baseline')
    out_dir.mkdir(exist_ok=True)
    with open(out_dir / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    T = audit_utr.table
    fastas = {k: audit_utr.read_fastas(work / k) for k in res}
    out = [f'baseline commit: {baseline_commit}; junctions: {len(ref)}; soft-masked genome: {Path(softmasked).name}\n']
    for name, r in res.items():
        exc = [k for k in r if r[k]['fate'].startswith('EXCEPTION')]
        out.append(f'{name}: junctions raising an exception: {len(exc)} '
                   f'({len({r[k]["gene_symbol"] for k in exc})} genes); first in run order: '
                   f'{exc[0] if exc else "-"}\n')
    soft = res['baseline_softmasked']
    exc = [k for k in soft if soft[k]['fate'].startswith('EXCEPTION')]
    slices = res['baseline_softmasked_no_gtf_canonical']
    out.append(T('baseline_softmasked_no_gtf_canonical exceptions: exception type and location', collections.Counter(
        (slices[k]['fate'].split(':')[0], slices[k]['where']) for k in slices if slices[k]['where'])))
    out.append(T('baseline_softmasked exceptions: unmasked fate of the junction',
                 collections.Counter(ref[k]['fate'] for k in exc)))
    out.append(T('baseline_softmasked exceptions: exception type and location',
                 collections.Counter((soft[k]['fate'].split(':')[0], soft[k]['where']) for k in exc)))
    out.append(T('baseline_softmasked: fate of junctions not raising, compared with unmasked',
                 collections.Counter((ref[k]['fate'], soft[k]['fate']) for k in soft if k not in exc)))
    for name, reference in [('baseline_hardmasked', 'baseline_unmasked'),
                            ('baseline_hardmasked_no_gtf_canonical', 'baseline_unmasked_no_gtf_canonical')]:
        hard, r0 = res[name], res[reference]
        out.append(T(f'{name}: fate transitions from {reference} (changed only)', collections.Counter(
            (r0[k]['fate'], hard[k]['fate'].split(':')[0]) for k in hard if hard[k]['fate'] != r0[k]['fate'])))
    for label, (n, ng, hit, nghit, lo, tot) in softmasked_cds(softmasked).items():
        out.append(f'chr15 {label}: {n} transcripts ({ng} genes); CDS with any soft-masked base: {hit} '
                   f'({nghit} genes); soft-masked CDS bases: {lo:,} of {tot:,} ({lo / tot:.2%})\n')
    if verify:
        for name in ['candidate_unmasked', 'candidate_softmasked']:
            same_fates = all(res[name][k]['fate'] == ref[k]['fate'] for k in ref)
            out.append(f'{name}: fates identical to baseline_unmasked: {same_fates}; FASTA files identical: '
                       f'{fastas[name] == fastas["baseline_unmasked"]}\n')
        for name, reference in [('candidate_hardmasked', 'baseline_hardmasked'),
                                ('candidate_hardmasked_no_gtf_canonical', 'baseline_hardmasked_no_gtf_canonical')]:
            out.append(T(f'{name}: fate transitions from {reference} (changed only)', collections.Counter(
                (res[reference][k]['fate'].split(':')[0], res[name][k]['fate'])
                for k in ref if res[name][k]['fate'] != res[reference][k]['fate'])))
        out.append(T('FASTA records (file | baseline_unmasked | candidate_softmasked)', collections.Counter(
            {(n, len(fastas['baseline_unmasked'].get(n, {})), len(fastas['candidate_softmasked'].get(n, {}))): 1
             for n in sorted(set(fastas['baseline_unmasked']) | set(fastas['candidate_softmasked']))})))
    (out_dir / 'summary.txt').write_text('\n'.join(out))
    print((out_dir / 'summary.txt').read_text())


if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[1] == '--worker':
        worker(*sys.argv[2:])
    else:
        parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
        parser.add_argument('--softmasked', type=Path, required=True, help='Ensembl 89 dna_sm chromosome 15 FASTA')
        parser.add_argument('--baseline', default=BASELINE_COMMIT, help='commit to compare against')
        parser.add_argument('--verify', action='store_true', help='also run the working tree (after the fix)')
        a = parser.parse_args()
        main_(a.softmasked.resolve(), a.baseline, a.verify)
