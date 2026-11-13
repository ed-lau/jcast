# -*- coding: utf-8 -*-

"""Audit of PR 4a2: building the RI retained-intron slice from the retained exon alone (REVIEW.md HIGH-11).

Runs every junction in tests/data/rmats through main._translate_one twice, each in its own process, using the
harness of reports/utr_audit/audit_utr.py:

    baseline    a git worktree of BASELINE_COMMIT (PR 4a: slice 1 = anc + alt1 + down)
    candidate   this repository's working tree (PR 4a2: slice 1 = alt1)

and writes:

    junctions.tsv     one row per RI junction: trimmed coordinates, slice lengths and fate under both versions
    fasta_diff.tsv    every output sequence present under one version but not the other, with the annotated
                      transcripts (name:biotype) it matches exactly
    summary.txt       the aggregate numbers quoted in the report, including a check that no non-RI junction changed

Usage (from the repository root):
    python reports/ri_audit/audit_ri.py [--baseline COMMIT]
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
BASELINE_COMMIT = '89e1a82'

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'reports' / 'utr_audit'))
sys.path.insert(0, str(REPO / 'reports' / 'tx1_audit'))
sys.path.insert(0, str(REPO / 'reports' / 'cterm_audit'))

import audit_utr  # noqa: E402  (worker harness)


def length(es, ee):
    return ee - es + 1 if es > 0 and ee > 0 else 0


def main_(baseline_commit):
    work = Path(tempfile.mkdtemp(prefix='ri_audit_'))
    baseline_repo = work / 'baseline'
    subprocess.run(['git', '-C', str(REPO), 'worktree', 'add', '--detach', str(baseline_repo), baseline_commit],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        base = audit_utr.run_variant(baseline_repo, work / 'out_baseline')
        cand = audit_utr.run_variant(REPO, work / 'out_candidate')
    finally:
        subprocess.run(['git', '-C', str(REPO), 'worktree', 'remove', '--force', str(baseline_repo)], check=True)

    import audit_cterm
    from jcast.annots import ReadAnnotations, ReadGenome
    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    gtf_copy = work / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
    shutil.copy(DATA / 'genome' / gtf_copy.name, gtf_copy)
    gtf = ReadAnnotations(logger=logger, path=gtf_copy)
    genome = ReadGenome(logger=logger, path=DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    proteins = audit_cterm.annotated_proteins(gtf, genome)

    # non-RI junctions must be untouched
    other_changed = [k for k, b in base.items() if k[0] != 'RI' and
                     (b['fate'], b['after'], b['phase']) != (cand[k]['fate'], cand[k]['after'], cand[k]['phase'])]

    rows = []
    for key, b in base.items():
        if key[0] != 'RI':
            continue
        c = cand[key]
        a = c['after']
        rows.append({
            'rmats_id': key[1], 'gene_symbol': c['gene_symbol'], 'gene_id': c['gene_id'], 'strand': c['strand'],
            'sum_sjc': c['sum_sjc'], 'tx0': c['tx0'], 'tx1': c['tx1'], 'phase': c['phase'],
            'anc_len': length(a['anc_es'], a['anc_ee']), 'retained_len': length(a['alt1_es'], a['alt1_ee']),
            'down_len': length(a['down_es'], a['down_ee']),
            'baseline_slice1_len': sum(length(a[f'{e}_es'], a[f'{e}_ee']) for e in ['anc', 'alt1', 'down']),
            'candidate_slice1_len': length(a['alt1_es'], a['alt1_ee']),
            'slice2_len': length(a['anc_es'], a['anc_ee']) + length(a['down_es'], a['down_ee']),
            'baseline_fate': b['fate'], 'candidate_fate': c['fate'],
            'coords_identical': b['after'] == a,
            **{f'{k}': v for k, v in a.items()},
        })
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

    def match_class(d):
        if not d['matches']:
            return 'no annotated match'
        return ('annotated protein_coding' if ':protein_coding' in d['matches']
                else 'annotated ' + '/'.join(sorted({m.split(':')[1] for m in d['matches'].split(';')})))

    def is_ri(d):
        return '|RI1|' in d['header'] or '|RI2|' in d['header']

    T = audit_utr.table
    out = [f'baseline commit: {baseline_commit}; RI junctions: {len(rows)}; '
           f'non-RI junctions whose fate, coordinates or phase changed: {len(other_changed)}\n']
    out.append(T('RI fate: baseline', collections.Counter(r['baseline_fate'] for r in rows)))
    out.append(T('RI fate: candidate', collections.Counter(r['candidate_fate'] for r in rows)))
    out.append(T('RI fate transitions (baseline | candidate)', collections.Counter(
        (r['baseline_fate'], r['candidate_fate']) for r in rows if r['baseline_fate'] != r['candidate_fate'])))
    out.append(f'RI junctions whose trimmed coordinates differ: {sum(not r["coords_identical"] for r in rows)}\n')
    out.append(f'retained exon == anc_es..down_ee before trimming: see reports/2026-11-10-ri-slice-audit.md §2\n')
    analysed = [r for r in rows if not r['candidate_fate'].startswith('SKIPPED')]
    out.append(T('slice 1 length, analysed RI junctions: baseline anc+retained+down vs candidate retained',
                 collections.Counter({'baseline total nt': sum(r['baseline_slice1_len'] for r in analysed),
                                      'candidate total nt': sum(r['candidate_slice1_len'] for r in analysed)})))
    out.append(T('intron length mod 3 (frameshift test), analysed RI junctions', collections.Counter(
        (r['retained_len'] - r['slice2_len']) % 3 for r in analysed)))
    out.append(T('FASTA records (file | baseline | candidate)', collections.Counter(
        {(n, len(fb.get(n, {})), len(fc.get(n, {}))): 1 for n in sorted(set(fb) | set(fc))})))
    out.append(T('FASTA diff (file | change | slice | annotation)', collections.Counter(
        (d['file'], d['change'], 'RI' if is_ri(d) else 'other', match_class(d)) for d in diff)))
    out.append(T('RI records in every output file, by annotation: baseline', collections.Counter(
        (n, match_class({'matches': ';'.join(proteins.get(s, []))})) for n, seqs in fb.items()
        for s, h in seqs.items() if '|RI1|' in h or '|RI2|' in h)))
    out.append(T('RI records in every output file, by annotation: candidate', collections.Counter(
        (n, match_class({'matches': ';'.join(proteins.get(s, []))})) for n, seqs in fc.items()
        for s, h in seqs.items() if '|RI1|' in h or '|RI2|' in h)))
    (OUT / 'summary.txt').write_text('\n'.join(out))
    print((OUT / 'summary.txt').read_text())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--baseline', default=BASELINE_COMMIT, help='commit to compare against')
    main_(parser.parse_args().baseline)
