# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §17.4): the maintainer's runs on folders made with
--statoff, against their runs on the same samples' folders made with statistics.

For each data set, two rMATS folders and the JCAST run on each (the same options, -m):
    rMATS       the same b1.txt and b2.txt; every event of the folder made with statistics in the --statoff folder,
                every column equal but PValue and FDR; every FDR NA
    JCAST       the -m cutoff, both runs; every junction of the earlier run in the new one, with the same read count;
                at or above the earlier cutoff, and below the new one, every column equal but the FDR and new_records;
                between the two, the junctions skipped before and what they give now; the junctions new to the folder;
                the events the new log says were left out, against the folder's events less its junctions
    FASTA       each output FASTA's sequences, and whether it keeps every earlier one; the Swiss-Prot isoforms
                (accession -N, less any sequence that is also a canonical) in the psq_ files, as cutoff_recovery.py
                counts them

Writes <out>/reruns.txt.

Usage (from the repository root):
    python reports/model_audit/check_reruns.py --set <label>=<rMATS --statoff>,<JCAST on it>,<rMATS with statistics>,
        <JCAST on it>,<Swiss-Prot FASTA> [--set ...] [--out <dir>]
"""

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

TYPES = ['SE', 'MXE', 'RI', 'A5SS', 'A3SS']
KEY = ['junction_type', 'junction_id']


def sha1(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()


def sequences(path) -> set:
    return {sha1(str(r.seq)) for r in SeqIO.parse(path, 'fasta')}


def logged(log, pattern) -> list:
    return re.findall(pattern, Path(log).read_text())


def rmats(new, old) -> tuple:
    """ the lines, and the --statoff folder's events """
    same = all((new / f).read_text().strip() == (old / f).read_text().strip() for f in ('b1.txt', 'b2.txt'))
    lines, events = [f'b1.txt and b2.txt as before: {same}'], 0
    for t in TYPES:
        a = pd.read_csv(new / f'{t}.MATS.JC.txt', sep='\t', dtype=str)
        b = pd.read_csv(old / f'{t}.MATS.JC.txt', sep='\t', dtype=str)
        cols = [c for c in b.columns if c not in ('PValue', 'FDR')]
        m = b[cols].merge(a[cols], on='ID', how='left', suffixes=('', '_new'), indicator=True)
        differ = sum(int((m[c] != m[c + '_new']).sum()) for c in cols if c != 'ID')
        lines.append(f"{t}: {len(a)} events, {len(b)} with statistics; all of those in the --statoff folder: "
                     f"{(m['_merge'] == 'both').all()}, every column but PValue and FDR equal: {differ == 0}; "
                     f"FDR all NA: {a['FDR'].isna().all()}")
        events += len(a)
    return lines, events


def jcast(new, old, events) -> list:
    c_new, c_old = (int(logged(d / 'logfile.log', r'modeled read count cutoff (\d+)')[0]) for d in (new, old))
    a = pd.read_csv(new / 'junctions.tsv', sep='\t', dtype=str, keep_default_na=False)
    b = pd.read_csv(old / 'junctions.tsv', sep='\t', dtype=str, keep_default_na=False)
    left = sum(int(n) for n in logged(new / 'logfile.log', r'rMATS \w+: (\d+) events left out'))
    records = [logged(d / 'logfile.log', r'and (\d+) JCAST records')[0] for d in (old, new)]
    lines = [f'-m: {c_old} -> {c_new}; junctions: {len(b)} -> {len(a)} (+{len(a) - len(b)}); FDR all NA: '
             f"{(a['fdr'] == 'NA').all()}; left out (log) {left}, the folder's events less its junctions "
             f'{events - len(a)}; database records (log): {records[0]} -> {records[1]}']
    m = b.merge(a, on=KEY, how='left', suffixes=('', '_new'), indicator=True)
    lines.append(f"every earlier junction in the new run: {(m['_merge'] == 'both').all()}; with the same read count: "
                 f"{(m['min_jc'] == m['min_jc_new']).all()}")
    n = m['min_jc'].astype(int)
    cols = [c for c in b.columns if c not in KEY + ['fdr', 'new_records', 'min_jc']]
    for label, rows in [(f'at or above {c_old}', m[n >= c_old]), (f'below {c_new}', m[n < c_new])]:
        differ = [c for c in cols if not (rows[c] == rows[c + '_new']).all()]
        lines.append(f"  {label}: {len(rows)} junctions, every column but fdr and new_records equal: "
                     f"{'yes' if not differ else 'no: ' + ', '.join(differ)}")
    band = m[(n >= c_new) & (n < c_old)]
    lines.append(f"  {c_new} to {c_old - 1}: {len(band)} junctions, all skipped_low_count before: "
                 f"{(band['fate'] == 'skipped_low_count').all()}; now " +
                 ', '.join(f'{f} {k}' for f, k in band['fate_new'].value_counts().items()))
    added = a.merge(b[KEY], on=KEY, how='left', indicator=True)
    added = added[added['_merge'] == 'left_only']
    above = added[added['min_jc'].astype(int) >= c_new]
    lines.append(f'  new to the folder: {len(added)} junctions; at or above {c_new}: {len(above)} (' +
                 ', '.join(f'{f} {k}' for f, k in above['fate'].value_counts().items()) + ')')
    return lines


def fastas(new, old, uniprot) -> list:
    iso, can = set(), set()
    for r in SeqIO.parse(uniprot, 'fasta'):
        (iso if '-' in r.id.split('|')[1] else can).add(sha1(str(r.seq)))
    iso -= can
    found = {}
    for d in (old, new):
        found[d] = set().union(*(sequences(f) & iso for f in d.glob('psq_*.fasta') if f.stem != 'psq_canonical'))
    lines = [f'Swiss-Prot isoforms in the psq_ files: {len(found[old])} -> {len(found[new])}, all earlier ones kept: '
             f'{found[old] <= found[new]}']
    for name in sorted(p.name for p in new.glob('*.fasta')):
        s_new, s_old = sequences(new / name), sequences(old / name)
        lines.append(f'  {name}: {len(s_old)} -> {len(s_new)} sequences, all earlier ones kept: {s_old <= s_new}')
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set', action='append', required=True,
                        help='<label>=<rMATS --statoff>,<JCAST>,<rMATS with statistics>,<JCAST>,<Swiss-Prot FASTA>')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_24')
    args = parser.parse_args()

    lines = []
    for s in args.set:
        label, paths = s.split('=', 1)
        rmats_new, jcast_new, rmats_old, jcast_old, uniprot = (Path(p) for p in paths.split(','))
        found, events = rmats(rmats_new, rmats_old)
        block = [f'== {label}', f'--statoff: {rmats_new}, {jcast_new}', f'statistics: {rmats_old}, {jcast_old}',
                 *found, *jcast(jcast_new, jcast_old, events), *fastas(jcast_new, jcast_old, uniprot), '']
        print('\n'.join(block), flush=True)
        lines += block

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'reruns.txt').write_text('\n'.join(lines).rstrip() + '\n')


if __name__ == '__main__':
    main()
