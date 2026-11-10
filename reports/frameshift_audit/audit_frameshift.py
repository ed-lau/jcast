# -*- coding: utf-8 -*-

"""Audit of PR 4b2: the T1/T2 frameshift test accounts for each slice's own reading frame.

Predicts, from the pipeline's own slices and frames, which T1/T2 junctions the old test (raw slice lengths) and the
new test (coding lengths: slice length minus its frame) disagree on. Then runs every junction through
main._translate_one at BASELINE_COMMIT and in the working tree (harness of reports/phase_audit/audit_phase.py) and
checks that exactly the predicted junctions change fate, in the predicted direction. Writes:

    junctions.tsv     the predicted junctions: frames, slice lengths, old and new test, fates before and after
    fasta_diff.tsv    every output sequence present under one version but not the other, with exact annotation matches
    summary.txt       the aggregate numbers quoted in the report

Usage (from the repository root):
    python reports/frameshift_audit/audit_frameshift.py [--baseline COMMIT]
"""

import argparse
import collections
import csv
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
BASELINE_COMMIT = '5c0a697'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']

sys.path.insert(0, str(REPO))
for d in ['phase_audit', 'utr_audit', 'tx1_audit', 'cterm_audit']:
    sys.path.insert(0, str(REPO / 'reports' / d))


def predict(gtf, genome, rmats):
    """T1/T2 junctions whose slices have different frames, with both tests evaluated on the actual slices."""
    from jcast.junctions import Junction
    from jcast.sequences import Sequence

    out = {}
    for jt, rma in zip(TYPES, rmats.rmats_list):
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            if jx.sum_sjc <= 1:
                continue
            jx.trim_cds(gtf)
            jx.get_translated_phase(gtf)
            p1, p2 = jx.slice_phases
            if p1 is None or p2 is None or p1 == p2 or jx.event_outside_cds:
                continue
            sq = Sequence(jx)
            sq.make_slice_localgenome(genome.genome)
            if any(m in sq.slice1_nt or m in sq.slice2_nt for m in 'NnXx'):
                continue
            sq.translate(log=False)
            if not (sq.slice1_aa and sq.slice2_aa):
                continue  # only junctions where both slices translate reach the T1/T2 test
            n1, n2 = len(sq.slice1_nt), len(sq.slice2_nt)
            raw, coding = (n1 - n2) % 3 != 0, ((n1 - p1) - (n2 - p2)) % 3 != 0
            out[(jt, jx.name)] = {'gene_symbol': jx.gene_symbol, 'strand': jx.strand, 'phases': (p1, p2),
                                  'slice_nt': (n1, n2), 'raw_frameshift': raw, 'coding_frameshift': coding}
    return out


def main_(baseline_commit):
    import audit_phase
    import audit_utr
    import audit_cterm
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from jcast import params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import RmatsResults
    from jcast.sequences import Sequence

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))

    work = Path(tempfile.mkdtemp(prefix='frameshift_audit_'))
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), baseline_commit],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        base = audit_phase.run_variant(baseline_repo, work / 'out_baseline', False)
        cand = audit_phase.run_variant(REPO, work / 'out_candidate', False)
    finally:
        subprocess.run(['git', '-C', str(REPO), 'worktree', 'remove', '--force', str(baseline_repo)], check=True)

    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(logger=logger, rmats_dir=DATA / 'rmats')
    proteins = audit_cterm.annotated_proteins(gtf, genome)
    pred = predict(gtf, genome, rmats)

    tier = {True: 'SUCCESS TIER 2', False: 'SUCCESS TIER 1'}
    disagree = {k: v for k, v in pred.items() if v['raw_frameshift'] != v['coding_frameshift']}
    rows = []
    for k, v in sorted(disagree.items()):
        rows.append({'type': k[0], 'rmats_id': k[1], 'gene_symbol': v['gene_symbol'], 'strand': v['strand'],
                     'slice1_phase': v['phases'][0], 'slice2_phase': v['phases'][1],
                     'slice1_nt': v['slice_nt'][0], 'slice2_nt': v['slice_nt'][1],
                     'old_test': tier[v['raw_frameshift']], 'new_test': tier[v['coding_frameshift']],
                     'baseline_fate': base[k]['fate'], 'candidate_fate': cand[k]['fate']})
    with open(OUT / 'junctions.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    fb, fc = audit_utr.read_fastas(work / 'out_baseline'), audit_utr.read_fastas(work / 'out_candidate')
    diff = []
    for name in sorted(set(fb) | set(fc)):
        b_seqs, c_seqs = fb.get(name, {}), fc.get(name, {})
        for change, seqs, other in [('lost', b_seqs, c_seqs), ('gained', c_seqs, b_seqs)]:
            for seq in sorted(set(seqs) - set(other)):
                diff.append({'file': name, 'change': change, 'header': seqs[seq], 'length': len(seq),
                             'matches': ';'.join(proteins.get(seq, [])), 'sequence': seq})
    with open(OUT / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'change', 'header', 'length', 'matches', 'sequence'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(diff)

    changed = {k for k in base if base[k]['fate'] != cand[k]['fate']}
    as_predicted = all(cand[k]['fate'] == tier[disagree[k]['coding_frameshift']] and
                       base[k]['fate'] == tier[disagree[k]['raw_frameshift']] for k in disagree)
    all_lost = {d['sequence'] for d in diff if d['change'] == 'lost'}
    all_gained = {d['sequence'] for d in diff if d['change'] == 'gained'}
    T = audit_utr.table
    out = [f'baseline commit: {baseline_commit}; junctions: {len(base)}\n',
           f'T1/T2-eligible junctions whose slices have different frames: {len(pred)}; '
           f'old and new frameshift test disagree: {len(disagree)}\n',
           f'junctions whose fate changed: {len(changed)}; identical to the predicted set: {changed == set(disagree)}; '
           f'every change in the predicted direction: {as_predicted}\n',
           f'sequences lost from every file: {len(all_lost - all_gained)}; '
           f'new sequences in no file before: {len(all_gained - all_lost)}\n']
    out.append(T('fate transitions (baseline | candidate)', collections.Counter(
        (base[k]['fate'], cand[k]['fate']) for k in changed)))
    out.append(T('FASTA records (file | baseline | candidate)', collections.Counter(
        {(n, len(fb.get(n, {})), len(fc.get(n, {}))): 1 for n in sorted(set(fb) | set(fc))})))
    out.append(T('FASTA diff (file | change | annotated protein_coding)', collections.Counter(
        (d['file'], d['change'], ':protein_coding' in d['matches']) for d in diff)))
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--baseline', default=BASELINE_COMMIT, help='commit to compare against')
    main_(parser.parse_args().baseline)
