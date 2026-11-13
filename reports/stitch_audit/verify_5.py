# -*- coding: utf-8 -*-

"""Verification of the PR 5 implementation (reports/2026-11-11-stitch-anchor-audit.md).

1. Every join the pipeline makes in the working tree equals the audited rule "unique + complete ends" at the
   working tree's anchor length (replayed by audit_stitch.join), and is classified against the genomic truth.
2. Every junction is run through main._translate_one at BASELINE_COMMIT and in the working tree (harness of
   reports/phase_audit/audit_phase.py): fate changes, and every FASTA record present in one but not the other, with
   its exact matches to annotated proteins.

Writes final/summary.txt, final/fates.tsv and final/fasta_diff.tsv.

Usage (from the repository root, at 2142711):
    python reports/stitch_audit/verify_5.py
"""

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
BASELINE_COMMIT = 'f1db483'

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(OUT))
for d in ['phase_audit', 'utr_audit', 'tx1_audit', 'cterm_audit']:
    sys.path.insert(0, str(REPO / 'reports' / d))


def main_():
    import audit_phase
    import audit_stitch
    import audit_utr
    import audit_cterm
    from jcast import params
    from jcast.annots import ReadAnnotations, ReadGenome

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    T = audit_utr.table
    k = params.aa_stitch_length

    work = Path(tempfile.mkdtemp(prefix='verify_5_'))

    # 1. the pipeline's joins equal the audited rule
    calls = audit_stitch.run(work / 'calls', 'gtf')
    seen, distinct = set(), []
    for c in calls:
        key = (c['canonical'], c['aa'], c['mode'])
        if key not in seen:
            seen.add(key)
            distinct.append(c)
    replay_mismatch = [c for c in distinct if c['product'] !=
                       audit_stitch.join(c['canonical'], c['aa'], c['mode'], k, 'unique+ends', c['ends_at_stop'],
                                         c['starts_at_start'])]
    verdicts = collections.Counter(audit_stitch.classify(c['product'], audit_stitch.truth(c)) for c in distinct)

    # 2. baseline vs working tree
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), BASELINE_COMMIT],
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
    proteins = audit_cterm.annotated_proteins(gtf, genome)

    final = OUT / 'final'
    final.mkdir(exist_ok=True)
    changed = sorted(key for key in base if base[key]['fate'] != cand[key]['fate'])
    with open(final / 'fates.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['type', 'rmats_id', 'gene_symbol', 'baseline', 'candidate'])
        for key in changed:
            w.writerow([*key, base[key]['gene_symbol'], base[key]['fate'], cand[key]['fate']])

    fb, fc = audit_utr.read_fastas(work / 'out_baseline'), audit_utr.read_fastas(work / 'out_candidate')
    diff = []
    for name in sorted(set(fb) | set(fc)):
        b_seqs, c_seqs = fb.get(name, {}), fc.get(name, {})
        for change, seqs, other in [('lost', b_seqs, c_seqs), ('gained', c_seqs, b_seqs)]:
            for seq in sorted(set(seqs) - set(other)):
                diff.append({'file': name, 'change': change, 'header': seqs[seq], 'length': len(seq),
                             'matches': ';'.join(proteins.get(seq, [])), 'sequence': seq})
    with open(final / 'fasta_diff.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'change', 'header', 'length', 'matches', 'sequence'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(diff)
    everywhere_b = {s for seqs in fb.values() for s in seqs}
    everywhere_c = {s for seqs in fc.values() for s in seqs}

    out = [f'baseline commit: {BASELINE_COMMIT}; candidate: working tree; anchor length {k}\n',
           f'1. distinct joins: {len(distinct)}; differing from the audited rule (unique + complete ends): '
           f'{len(replay_mismatch)}\n',
           T('pipeline joins against the genomic truth', verdicts),
           f'2. junctions whose fate changed: {len(changed)}\n',
           T('fate transitions (baseline | candidate)', collections.Counter(
               (base[key]['fate'], cand[key]['fate']) for key in changed)),
           T('FASTA records (file | baseline | candidate)', collections.Counter(
               {(n, len(fb.get(n, {})), len(fc.get(n, {}))): 1 for n in sorted(set(fb) | set(fc))})),
           T('FASTA diff (file | change | matches an annotated protein_coding protein)', collections.Counter(
               (d['file'], d['change'], ':protein_coding' in d['matches']) for d in diff)),
           f'annotated proteins in some file before but in none now: '
           f'{len([s for s in everywhere_b - everywhere_c if s in proteins])}; '
           f'now in some file but in none before: {len([s for s in everywhere_c - everywhere_b if s in proteins])}\n']
    (final / 'summary.txt').write_text('\n'.join(out))
    print((final / 'summary.txt').read_text())


if __name__ == '__main__':
    main_()
