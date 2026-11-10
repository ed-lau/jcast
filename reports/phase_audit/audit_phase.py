# -*- coding: utf-8 -*-

"""Audit of PR 4b: per-slice reading frame from the first retained coding exon, and the three-frame fallback.

Runs every junction in tests/data/rmats through main._translate_one in three configurations, each in its own
process:

    baseline     a git worktree of BASELINE_COMMIT (PR 4a2: one phase from a fixed exon matched exactly)
    fallback     this repository's working tree, params.three_frame_fallback = True
    one_frame    this repository's working tree, params.three_frame_fallback = False
    no_startless as one_frame, and params.frame_for_start_less_slice = False: a slice lacking the gene's start
                 codon that the other slice contains is given no reading frame

and writes:

    junctions.tsv     one row per junction: phases, whether three-frame translation was used, and fate under each
    fasta_diff.tsv    every output sequence that differs between baseline and either candidate, with the
                      annotated transcripts (name:biotype) it matches exactly
    summary.txt       the aggregate numbers quoted in the report

Canonical sequences come from the GTF only and UniProt is stubbed out, so the run is deterministic and needs
no network.

Usage (from the repository root):
    python reports/phase_audit/audit_phase.py [--baseline COMMIT]
"""

import argparse
import collections
import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
BASELINE_COMMIT = 'bc7d894'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']
VARIANTS = ['baseline', 'fallback', 'one_frame', 'no_startless']


# ---- worker: runs inside one code version ------------------------------------------------------------------

def worker(out_dir, fallback, startless='True'):
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
    params.use_gtf_only = True
    params.three_frame_fallback = fallback == 'True'
    params.frame_for_start_less_slice = startless == 'True'
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    # record whether any translation of the junction guessed the frame
    original = Sequence.translate

    def recording(self, use_phase=True, log=True):
        phases = getattr(self.j, 'slice_phases', (self.j.phase, self.j.phase))
        one_frame = use_phase and (all(p in [0, 1, 2] for p in phases) or not getattr(params, 'three_frame_fallback', True))
        if not one_frame:
            self.j.three_frame_used = True
        return original(self, use_phase=use_phase, log=log)
    Sequence.translate = recording

    out_dir = Path(out_dir)
    work = Path(tempfile.mkdtemp(prefix='phase_audit_'))
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask='NnXx', qvalue=[0, 1])

    rows = []
    for jt, rma in zip(TYPES, rmats.rmats_list):
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            fate = main._translate_one(jx, gtf, genome, args, pred_bound=args.read)
            rows.append({'type': jt, 'rmats_id': jx.name, 'gene_id': jx.gene_id, 'gene_symbol': jx.gene_symbol,
                         'strand': jx.strand, 'sum_sjc': int(jx.sum_sjc), 'fate': short_fate(fate),
                         'phase': jx.phase, 'slice_phases': list(getattr(jx, 'slice_phases', (jx.phase, jx.phase))),
                         'phase_ambiguous': getattr(jx, 'phase_ambiguous', None),
                         'three_frame_used': getattr(jx, 'three_frame_used', False),
                         'tx0': jx.tx0, 'tx1': jx.tx1,
                         'slice_exons': [list(map(list, e)) for e in jx.slice_exons()] if hasattr(jx, 'slice_exons') else None})
    (out_dir / 'junctions.json').write_text(json.dumps(rows, default=int))


def run_variant(repo, out_dir, fallback, startless=True):
    out_dir.mkdir(parents=True)
    env = dict(os.environ, PYTHONPATH=str(repo), MPLBACKEND='Agg')
    subprocess.run([sys.executable, __file__, '--worker', str(out_dir), str(fallback), str(startless)], env=env, check=True,
                   cwd=str(out_dir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return {(r['type'], r['rmats_id']): r for r in json.loads((out_dir / 'junctions.json').read_text())}


# ---- comparison --------------------------------------------------------------------------------------------

def main_(baseline_commit):
    sys.path.insert(0, str(REPO))
    for d in ['utr_audit', 'tx1_audit', 'cterm_audit']:
        sys.path.insert(0, str(REPO / 'reports' / d))
    import audit_utr
    import audit_cterm
    from jcast.annots import ReadAnnotations, ReadGenome

    work = Path(tempfile.mkdtemp(prefix='phase_audit_'))
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), baseline_commit],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        res = {'baseline': run_variant(baseline_repo, work / 'out_baseline', True),
               'fallback': run_variant(REPO, work / 'out_fallback', True),
               'one_frame': run_variant(REPO, work / 'out_one_frame', False),
               'no_startless': run_variant(REPO, work / 'out_no_startless', False, False)}
    finally:
        subprocess.run(['git', '-C', str(REPO), 'worktree', 'remove', '--force', str(baseline_repo)], check=True)

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    proteins = audit_cterm.annotated_proteins(gtf, genome)
    coding_genes = set(gtf.annot.query('feature == "CDS" & transcript_biotype == "protein_coding"')['gene_id'])

    def start_codon_cover(r):
        """For each slice, whether it contains the gene's annotated (most upstream) start codon."""
        start = r['tx0'] if r['strand'] == '+' else r['tx1']
        if start is None or r['slice_exons'] is None:
            return [None, None]
        return [any(es <= start <= ee for es, ee in exons) for exons in r['slice_exons']]

    rows = []
    for key, b in res['baseline'].items():
        f, o, n = res['fallback'][key], res['one_frame'][key], res['no_startless'][key]
        has_start = start_codon_cover(f)
        rows.append({
            'type': key[0], 'rmats_id': key[1], 'gene_symbol': b['gene_symbol'], 'gene_id': b['gene_id'],
            'strand': b['strand'], 'sum_sjc': b['sum_sjc'], 'gene_has_protein_coding_cds': b['gene_id'] in coding_genes,
            'baseline_phase': b['phase'], 'candidate_slice1_phase': f['slice_phases'][0],
            'candidate_slice2_phase': f['slice_phases'][1], 'candidate_phase_ambiguous': f['phase_ambiguous'],
            'slice1_has_start_codon': has_start[0], 'slice2_has_start_codon': has_start[1],
            'baseline_three_frame': b['three_frame_used'], 'fallback_three_frame': f['three_frame_used'],
            'baseline_fate': b['fate'], 'fallback_fate': f['fate'], 'one_frame_fate': o['fate'],
            'no_startless_fate': n['fate'],
        })
    with open(OUT / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    fastas = {v: audit_utr.read_fastas(work / f'out_{v}') for v in VARIANTS}

    def annotation(seq):
        m = proteins.get(seq, [])
        if not m:
            return 'no annotated match'
        return 'annotated protein_coding' if any(x.endswith(':protein_coding') for x in m) else \
            'annotated ' + '/'.join(sorted({x.split(':')[1] for x in m}))

    diff = []
    for v in ['fallback', 'one_frame', 'no_startless']:
        for name in sorted(set(fastas['baseline']) | set(fastas[v])):
            b_seqs, c_seqs = fastas['baseline'].get(name, {}), fastas[v].get(name, {})
            for change, seqs, other in [('lost', b_seqs, c_seqs), ('gained', c_seqs, b_seqs)]:
                for seq in sorted(set(seqs) - set(other)):
                    diff.append({'variant': v, 'file': name, 'change': change, 'header': seqs[seq],
                                 'length': len(seq), 'annotation': annotation(seq),
                                 'matches': ';'.join(proteins.get(seq, [])), 'sequence': seq})
    with open(OUT / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(diff[0]), delimiter='\t')
        w.writeheader()
        w.writerows(diff)

    T = audit_utr.table
    out = [f'baseline commit: {baseline_commit}; junctions: {len(rows)}\n']
    for v in VARIANTS:
        out.append(T(f'fate counts: {v}', collections.Counter(r[f'{v}_fate'] for r in rows)))
    for v in ['fallback', 'one_frame', 'no_startless']:
        out.append(T(f'fate transitions baseline -> {v}', collections.Counter(
            (r['baseline_fate'], r[f'{v}_fate']) for r in rows if r['baseline_fate'] != r[f'{v}_fate'])))
    out.append(T('fate transitions fallback -> one_frame', collections.Counter(
        (r['fallback_fate'], r['one_frame_fate']) for r in rows if r['fallback_fate'] != r['one_frame_fate'])))
    out.append(T('fate transitions one_frame -> no_startless', collections.Counter(
        (r['one_frame_fate'], r['no_startless_fate']) for r in rows if r['one_frame_fate'] != r['no_startless_fate'])))

    live = [r for r in rows if not r['baseline_fate'].startswith('SKIPPED') or not r['fallback_fate'].startswith('SKIPPED')]
    def phase_state(r):
        p1, p2 = r['candidate_slice1_phase'], r['candidate_slice2_phase']
        if p1 is None and p2 is None:
            return 'neither slice'
        if p1 is None or p2 is None:
            return 'one slice'
        return 'both, same' if p1 == p2 else 'both, different'
    out.append(T('candidate phase availability (junctions not skipped) | baseline had a phase', collections.Counter(
        (phase_state(r), r['baseline_phase'] is not None) for r in live)))
    agree = [r for r in live if r['baseline_phase'] is not None and r['candidate_slice1_phase'] is not None]
    out.append(f'baseline phase available and candidate slice-1 phase available: {len(agree)}; '
               f'they differ: {sum(r["baseline_phase"] != r["candidate_slice1_phase"] for r in agree)} '
               f'({", ".join(r["gene_symbol"] + " " + r["type"] + " " + r["rmats_id"] for r in agree if r["baseline_phase"] != r["candidate_slice1_phase"])})\n')
    out.append(f'candidate phase ambiguous (transcripts disagree): {sum(bool(r["candidate_phase_ambiguous"]) for r in live)}\n')
    out.append(T('slice with no annotated frame: does the gene have any protein-coding CDS?', collections.Counter(
        (phase_state(r), r['gene_has_protein_coding_cds']) for r in live if phase_state(r) in ('neither slice', 'one slice'))))
    def start_state(r):
        a, b = r['slice1_has_start_codon'], r['slice2_has_start_codon']
        if a is None:
            return 'no start codon annotated'
        return {(True, True): 'both slices', (False, False): 'neither slice'}.get((a, b), 'one slice only')
    out.append(T('annotated start codon contained in (junctions not skipped) | baseline fate | one_frame fate',
                 collections.Counter((start_state(r), r['baseline_fate'], r['one_frame_fate']) for r in live
                                     if r['baseline_fate'] != r['one_frame_fate'])))
    out.append(T('annotated start codon contained in, all junctions not skipped (one_frame fate)', collections.Counter(
        (start_state(r), r['one_frame_fate']) for r in live)))
    out.append(T('three-frame translation used (baseline, fallback)', collections.Counter(
        (r['baseline_three_frame'], r['fallback_three_frame']) for r in live)))
    out.append(T('fallback: fate of junctions that used three-frame translation', collections.Counter(
        r['fallback_fate'] for r in live if r['fallback_three_frame'])))
    out.append(T('FASTA records (file | baseline | fallback | one_frame)', collections.Counter(
        {(n, *[len(fastas[v].get(n, {})) for v in VARIANTS]): 1
         for n in sorted(set().union(*[set(fastas[v]) for v in VARIANTS]))})))
    out.append(T('FASTA diff vs baseline (variant | file | change | annotation)', collections.Counter(
        (d['variant'], d['file'], d['change'], d['annotation']) for d in diff)))
    for v in VARIANTS:
        out.append(T(f'all records by annotation, T1/T2 files: {v}', collections.Counter(
            (n, annotation(s)) for n, seqs in fastas[v].items() if n in ('psq_T1.fasta', 'psq_T2.fasta')
            for s in seqs)))
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())


if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[1] == '--worker':
        worker(*sys.argv[2:])
    else:
        parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
        parser.add_argument('--baseline', default=BASELINE_COMMIT, help='commit to compare against')
        main_(parser.parse_args().baseline)
