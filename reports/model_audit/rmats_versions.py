# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §5): the same BAMs through rMATS 4.1.0 and 4.4.0.

4.1.1 added --allow-clipping and fixed "novel event classification" (the release notes say no more); both runs here
used --novelSS and --allow-clipping. Per splice type: the events, the events rMATS calls novel, and, for the events
both versions report (matched by chromosome, strand and every coordinate column), how the read count JCAST uses
(Junction.min_jc: the smaller of the summed inclusion and skipping junction counts) and the summed junction reads
changed, and how the novel calls moved. (Every event in *.MATS.JC.txt has at least one read in each form, in both
versions: min_jc >= 1.)

Writes <out>/rmats_versions.tsv and prints it.

Usage:
    python reports/model_audit/rmats_versions.py --old <rMATS 4.1.0 dir> --new <rMATS 4.4.0 dir> --out <dir>
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

SPLICE_TYPES = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']


def summed(col: pd.Series) -> pd.Series:
    """ a column of comma-separated replicate counts, summed per row (anything not a number counts 0) """
    return pd.to_numeric(col.astype(str).str.split(',').explode(), errors='coerce').fillna(0).groupby(level=0).sum()


def events(rmats_dir, st) -> pd.DataFrame:
    m = pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t')
    coords = list(m.columns[3:list(m.columns).index('ID.1')])        # chr, strand and the coordinates
    ijc = summed(m['IJC_SAMPLE_1']) + summed(m['IJC_SAMPLE_2'])
    sjc = summed(m['SJC_SAMPLE_1']) + summed(m['SJC_SAMPLE_2'])
    out = pd.DataFrame({'key': m[coords].astype(str).agg(':'.join, axis=1), 'ID': m['ID'].astype(str),
                        'min_jc': np.minimum(ijc, sjc).astype(int), 'reads': (ijc + sjc).astype(int)})
    for kind in ['novelJunction', 'novelSpliceSite']:
        f = Path(rmats_dir) / f'fromGTF.{kind}.{st}.txt'
        ids = set(pd.read_csv(f, sep='\t', usecols=['ID'])['ID'].astype(str)) if f.exists() else set()
        out[kind] = out['ID'].isin(ids)
    out['call'] = np.select([out['novelJunction'], out['novelSpliceSite']], ['novelJunction', 'novelSpliceSite'],
                            'annotated')
    return out


def compare(old_dir, new_dir) -> pd.DataFrame:
    rows = []
    for st in SPLICE_TYPES:
        o, n = events(old_dir, st), events(new_dir, st)
        both = o.merge(n, on='key', suffixes=('_old', '_new'))
        moved = pd.crosstab(both['call_old'], both['call_new'])
        nz = both[both['min_jc_old'] > 0]
        row = {'splice_type': st, 'events 4.1.0': len(o), 'events 4.4.0': len(n), 'in both': len(both),
               'only 4.1.0': len(o) - len(both), 'only 4.4.0': len(n) - len(both)}
        for v, d in [('4.1.0', o), ('4.4.0', n)]:
            row[f'novelJunction {v}'] = int(d['novelJunction'].sum())
            row[f'novelSpliceSite {v}'] = int(d['novelSpliceSite'].sum())
            row[f'median min_jc, {v}'] = float(d['min_jc'].median())
        row['in both: novelJunction 4.1.0 -> annotated 4.4.0'] = int(moved.get('annotated', {}).get('novelJunction', 0))
        row['in both: annotated 4.1.0 -> novel 4.4.0'] = int(
            both[(both['call_old'] == 'annotated') & (both['call_new'] != 'annotated')].shape[0])
        row['in both: same call'] = int((both['call_old'] == both['call_new']).sum())
        row['in both: summed reads, new/old (total)'] = round(both['reads_new'].sum() / both['reads_old'].sum(), 4)
        row['in both: min_jc, new/old (total)'] = round(both['min_jc_new'].sum() / both['min_jc_old'].sum(), 4)
        row['in both, min_jc >= 1 in 4.1.0: share with a higher min_jc in 4.4.0'] = round(
            (nz['min_jc_new'] > nz['min_jc_old']).mean(), 4)
        row['in both, min_jc >= 1 in 4.1.0: share with a lower min_jc in 4.4.0'] = round(
            (nz['min_jc_new'] < nz['min_jc_old']).mean(), 4)
        rows.append(row)
    t = pd.DataFrame(rows)
    total = t.drop(columns='splice_type').sum()
    for c in t.columns:
        if 'median' in c or '/' in c or 'share' in c:
            total[c] = np.nan
    return pd.concat([t, pd.DataFrame([{'splice_type': 'all', **total.to_dict()}])], ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--old', required=True)
    p.add_argument('--new', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t = compare(a.old, a.new)
    t.to_csv(out / 'rmats_versions.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(t.set_index('splice_type').T.to_string())


if __name__ == '__main__':
    main()
