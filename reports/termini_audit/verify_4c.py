# -*- coding: utf-8 -*-

"""Verification of the PR 4c implementation against the alternative-termini simulation.

Runs every junction through main._translate_one at BASELINE_COMMIT (PR 4b2) and in the working tree (PR 4c), using
the harness of reports/phase_audit/audit_phase.py, and checks against expected_after_4b2/ (the simulation re-run on
PR 4b2, policy tsl2):

    1. every junction's fate equals the simulated fate;
    2. every routed slice is written where the simulation put it, and its record matches the same annotated
       transcripts. An anchored slice that cannot be joined goes to the orphan file of its anchor (the simulation
       calls it psq_T1_alt_orphan) as its unjoined peptide, compared with the simulated peptide. A product already
       written by a sibling event is not written again by the FASTA writer; it is looked up in the expected file;
    3. across all files, no sequence that matches an annotated protein is lost.

Writes final/summary.txt, final/slices.tsv and final/fasta_diff.tsv.

Usage (from the repository root, at the PR 4c commit):
    python reports/termini_audit/verify_4c.py
"""

import collections
import csv
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = REPO / 'tests' / 'data'
BASELINE_COMMIT = 'a0f3338'
HEADER = re.compile(r'\|(MXE|SE|RI|A5SS|A3SS)([12])\|(\d+)\|')

sys.path.insert(0, str(REPO))
for d in ['phase_audit', 'utr_audit', 'tx1_audit', 'cterm_audit']:
    sys.path.insert(0, str(REPO / 'reports' / d))


def records_by_slice(out_dir):
    """(type, rmats id, slice) -> list of (file stem, sequence) for every record carrying a junction id."""
    out = collections.defaultdict(list)
    for p in sorted(out_dir.glob('psq_*.fasta')):
        for r in SeqIO.parse(p, 'fasta'):
            m = HEADER.search(r.id)
            if m:
                out[(m.group(1), m.group(3), m.group(2))].append((p.stem, str(r.seq)))
    return out


def main_():
    import audit_phase
    import audit_utr
    import audit_cterm
    from jcast.annots import ReadAnnotations, ReadGenome

    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False

    work = Path(tempfile.mkdtemp(prefix='verify_4c_'))
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

    expected_dir = OUT / 'expected_after_4b2'
    with open(expected_dir / 'junctions.tsv') as fh:
        exp_j = {(r['type'], r['rmats_id']): r for r in csv.DictReader(fh, delimiter='\t')
                 if r['policy'] in ('tsl2', '-')}
    with open(expected_dir / 'slices.tsv') as fh:
        exp_s = [r for r in csv.DictReader(fh, delimiter='\t') if r['policy'] == 'tsl2' and r['routed'] == 'True']

    # 1. fates
    def expected_fate(k):
        r = exp_j[(k[0], str(k[1]))]
        return r['proposed_fate'] if r['analysed'] == 'True' else r['actual_fate']
    fate_mismatch = [k for k in cand if cand[k]['fate'].split(':')[0].split('.')[0] != expected_fate(k).split('.')[0]]

    # 2. routed slices
    fb, fc = audit_utr.read_fastas(work / 'out_baseline'), audit_utr.read_fastas(work / 'out_candidate')
    recs = records_by_slice(work / 'out_candidate')
    rows, slice_mismatch = [], []
    for e in exp_s:
        got = recs.get((e['type'], e['rmats_id'], e['slice']), [])
        files = sorted({f for f, _ in got})
        matches = sorted({m for _, s in got for m in proteins.get(s, [])})
        expected_file = e['proposed_file']
        expected_matches = sorted(m for m in e['product_matches'].split(';') if m)
        how = 'own record'
        if expected_file == 'psq_canonical':
            ok, how = not got, 'written as the canonical record (no junction id)'
        elif expected_file == 'psq_T1_alt_orphan':
            # orphans are written as the unjoined slice peptide: compare with the simulated peptide
            anchor = 'nterm' if e['n_anchored'] == 'True' else 'cterm'
            ok = files == [f'psq_T1_alt_{anchor}_orphan'] and [s for _, s in got] == [e['aa']]
        elif got:
            ok = files == [expected_file] and matches == expected_matches
        else:
            # the FASTA writer skips a sequence already in the file (written by a sibling event): look for it there
            in_file = {m for s in fc.get(expected_file + '.fasta', {}) for m in proteins.get(s, [])}
            ok, how = bool(expected_matches) and set(expected_matches) <= in_file, 'deduplicated (sibling event)'
        rows.append({'type': e['type'], 'rmats_id': e['rmats_id'], 'slice': e['slice'], 'gene_symbol': e['gene_symbol'],
                     'expected_file': expected_file, 'written_files': ';'.join(files) or '-', 'how': how,
                     'expected_matches': e['product_matches'], 'written_matches': ';'.join(matches), 'agrees': ok})
        if not ok:
            slice_mismatch.append(rows[-1])

    final = OUT / 'final'
    final.mkdir(exist_ok=True)
    with open(final / 'slices.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # 3. output files
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
    annotated_lost = [s for s in everywhere_b - everywhere_c if s in proteins]

    T = audit_utr.table
    out = [f'baseline commit: {BASELINE_COMMIT} (PR 4b2); candidate: working tree (PR 4c); junctions: {len(cand)}\n',
           f'1. junctions whose fate differs from the simulation: {len(fate_mismatch)}\n',
           f'2. routed slices: {len(rows)}; written where the simulation put them, with the same annotation: '
           f'{len(rows) - len(slice_mismatch)}\n',
           f'3. sequences matching an annotated protein present before but in no file now: {len(annotated_lost)}\n']
    out.append(T('routed slices: how the written record was found', collections.Counter(r['how'] for r in rows)))
    out.append(T('fate counts: candidate', collections.Counter(c['fate'] for c in cand.values())))
    out.append(T('fate transitions baseline -> candidate', collections.Counter(
        (base[k]['fate'], cand[k]['fate']) for k in cand if base[k]['fate'] != cand[k]['fate'])))
    out.append(T('FASTA records (file | baseline | candidate)', collections.Counter(
        {(n, len(fb.get(n, {})), len(fc.get(n, {}))): 1 for n in sorted(set(fb) | set(fc))})))
    out.append(T('FASTA diff (file | change | annotated protein_coding)', collections.Counter(
        (d['file'], d['change'], ':protein_coding' in d['matches']) for d in diff)))
    if fate_mismatch or slice_mismatch:
        out.append('MISMATCHES\n' + '\n'.join(map(str, fate_mismatch + slice_mismatch)) + '\n')
    (final / 'summary.txt').write_text('\n'.join(out))
    print((final / 'summary.txt').read_text())


if __name__ == '__main__':
    main_()
