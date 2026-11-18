# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §10): do premature-stop junctions with few reads respond
to the NMD inhibitor as those with many do? If they do, they come from transcripts, not from noise.

The human cells are DMSO (rMATS sample 1, b1.txt) and NMDI14 (sample 2). Per event and condition, the inclusion form's
share (PSI) from the summed replicate counts, each form per junction position (rMATS's IncFormLen and SkipFormLen):
(I / IncFormLen) / (I / IncFormLen + S / SkipFormLen); an event with no read in one condition has no PSI there. A T4
junction's slice k hit a premature stop (junctions.tsv's fate_detail, "slice k"), and slice 1 is the inclusion form in
every splice type (sequences.make_slice_localgenome). The stop form's change is PSI(NMDI14) - PSI(DMSO) for slice 1,
and its negative for slice 2: positive when the form gains under the inhibitor. T1 junctions (no stop, no frameshift)
are the baseline: for them, the inclusion form's change, and its negative, average the same if NMDI14 favours neither.

Writes, into --out:
    nmd.tsv             per group (T4 with the stop in the inclusion form, in the skipping form, either; T1), rMATS's
                        call (annotated, novel, both) and read-count bin (Junction.min_jc, summed over all six
                        samples): the events with a PSI in both conditions, the mean change of the stop form (T1: the
                        inclusion form) with a 95% interval, and the share that rise of those that change; of the
                        events rMATS finds changed (FDR < 0.05), how many and the share that rise
    nmd_genes.tsv       T4 junctions by whether their gene has a transcript in the folder's DESeq2 list of NMD targets
                        (--targets, R/01_out/DTE_significant_NMD_targets.csv)
    nmd_gencode.tsv     with --gtf, a positive control: the events whose one form is only in GENCODE's
                        nonsense_mediated_decay transcripts, the other in a protein_coding one (each of a form's
                        junctions an intron of the transcript; RI left out), and the change of the NMD form
    nmd_transcripts.tsv with --dte (R/01_out/DTE_results_NMD_vs_DMSO.csv), the transcripts by GENCODE type (baseMean
                        >= 10): whether NMD transcripts as a class rise under NMDI14 more than coding ones
    nmd.png             the mean change by read count: JCAST's T4 stop forms, GENCODE's NMD forms, T1 (baseline)

Usage:
    python reports/model_audit/nmd_response.py --rmats <rMATS 4.4.0 dir> --run <-r 1 output dir>
        [--targets <DTE_significant_NMD_targets.csv>] [--dte <DTE_results_NMD_vs_DMSO.csv>] [--gtf <gtf>
        --cache_dir <dir>] --out <dir>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cutoff_recovery as cr  # noqa: E402
from rmats_versions import summed  # noqa: E402

BINS = [1, 3, 8, 16, 32, 128, 1024, np.inf]
LABELS = ['1-2', '3-7', '8-15', '16-31', '32-127', '128-1023', '>=1024']


def events(rmats_dir) -> pd.DataFrame:
    """ per event: PSI in each condition, rMATS's FDR, and its call """
    out = []
    for st in cr.SPLICE_TYPES:
        m = pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t')
        novel = set()
        for kind in ['novelJunction', 'novelSpliceSite']:
            f = Path(rmats_dir) / f'fromGTF.{kind}.{st}.txt'
            if f.exists():
                novel |= set(pd.read_csv(f, sep='\t', usecols=['ID'])['ID'].astype(str))
        ids = m['ID'].astype(str)
        psi = {}
        for k, s in [('DMSO', 1), ('NMDI14', 2)]:
            inc = summed(m[f'IJC_SAMPLE_{s}']) / m['IncFormLen']
            skp = summed(m[f'SJC_SAMPLE_{s}']) / m['SkipFormLen']
            psi[k] = (inc / (inc + skp)).where(inc + skp > 0)
        out.append(pd.DataFrame({'junction_type': st, 'junction_id': ids, 'call': np.where(ids.isin(novel), 'novel',
                                                                                             'annotated'),
                                 'psi_dmso': psi['DMSO'], 'psi_nmdi14': psi['NMDI14'], 'rmats_fdr': m['FDR']}))
    return pd.concat(out, ignore_index=True)


def table(jx) -> pd.DataFrame:
    """ the junctions to compare, each with its group and the change of the form that counts """
    t4 = jx[jx['fate'] == 'T4'].copy()
    t4['group'] = np.where(t4['fate_detail'] == 'slice 1', 'T4, stop in the inclusion form',
                           'T4, stop in the skipping form')
    t4['change'] = np.where(t4['fate_detail'] == 'slice 1', 1, -1) * (t4['psi_nmdi14'] - t4['psi_dmso'])
    t1 = jx[jx['fate'] == 'T1'].copy()
    t1['group'] = 'T1 (the inclusion form)'
    t1['change'] = t1['psi_nmdi14'] - t1['psi_dmso']
    t = pd.concat([t4, t4.assign(group='T4, either form'), t1], ignore_index=True)
    t = t[t['change'].notna()]
    t['bin'] = pd.cut(t['min_jc'], BINS, right=False, labels=LABELS)
    return t


def summarise(g) -> dict:
    ch = g['change']
    moved = ch[ch != 0]
    sig = g[g['rmats_fdr'] < 0.05]['change']
    half = 1.96 * ch.std(ddof=1) / np.sqrt(len(ch)) if len(ch) > 1 else np.nan
    return {'events': len(ch), 'mean change': round(ch.mean(), 4), '95% from': round(ch.mean() - half, 4),
            '95% to': round(ch.mean() + half, 4), 'share rising': round((moved > 0).mean(), 4) if len(moved) else None,
            'FDR < 0.05': len(sig), 'FDR < 0.05: share rising': round((sig > 0).mean(), 4) if len(sig) else None}


def by_bin(t, by) -> pd.DataFrame:
    rows = []
    for keys, g in t.groupby(by, observed=True):
        rows.append({**dict(zip(by, keys)), 'bin': 'all', **summarise(g)})
        for b, gb in g.groupby('bin', observed=True):
            rows.append({**dict(zip(by, keys)), 'bin': b, **summarise(gb)})
    return pd.DataFrame(rows)


def form_junctions(m, st) -> tuple:
    """ per event, the splice junctions of each form, (left, right) in rMATS's coordinates; RI's inclusion form has
    none, so RI is left out """
    z = lambda *cols: zip(*(m[c] for c in cols))
    if st == 'SE':
        return ([[(u, s), (e, d)] for u, s, e, d in z('upstreamEE', 'exonStart_0base', 'exonEnd', 'downstreamES')],
                [[(u, d)] for u, d in z('upstreamEE', 'downstreamES')])
    if st == 'MXE':
        cols = ('upstreamEE', '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd', 'downstreamES')
        rows = list(z(*cols))
        return ([[(u, s1), (e1, d)] for u, s1, e1, s2, e2, d in rows],
                [[(u, s2), (e2, d)] for u, s1, e1, s2, e2, d in rows])
    right = m['flankingES'] >= m['longExonEnd']          # A5SS, A3SS: the flanking exon after the alternative one
    rows = list(z('longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE'))
    long_ = [[(le, fs)] if r else [(fe, ls)] for r, (ls, le, ss, se, fs, fe) in zip(right, rows)]
    short = [[(se, fs)] if r else [(fe, ss)] for r, (ls, le, ss, se, fs, fe) in zip(right, rows)]
    return long_, short


def poison_forms(rmats_dir, gtf, cache_dir) -> pd.DataFrame:
    """ the events whose one form is in GENCODE's nonsense_mediated_decay transcripts only, and the other in a
    protein_coding transcript: a form is in a transcript if each of its junctions is one of the transcript's introns """
    import logging
    from jcast import cache
    from jcast.annots import ReadAnnotations
    cache.set_directory(cache_dir)
    annot = ReadAnnotations(logging.getLogger('jcast-audit'), gtf).annot
    ex = annot.loc[annot['feature'] == 'exon', ['seqname', 'strand', 'transcript_id', 'transcript_type', 'start',
                                                'end']].sort_values(['transcript_id', 'start'])
    del annot
    nxt = ex.groupby('transcript_id')['start'].shift(-1)
    has = nxt.notna()
    introns = {}
    for c, s, t, a, b in zip(ex.loc[has, 'seqname'], ex.loc[has, 'strand'], ex.loc[has, 'transcript_id'],
                             ex.loc[has, 'end'], nxt[has]):
        introns.setdefault((c, s, int(a), int(b) - 1), set()).add(t)
    tx_type = dict(zip(ex['transcript_id'], ex['transcript_type']))
    del ex

    def types(c, s, js):
        sets = [introns.get((c, s, int(a), int(b)), set()) for a, b in js]
        return {tx_type[t] for t in set.intersection(*sets)} if sets else set()

    out = []
    for st in ['SE', 'MXE', 'A5SS', 'A3SS']:
        m = pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t')
        inc, skp = form_junctions(m, st)
        for i, (c, s, ji, js) in enumerate(zip(m['chr'], m['strand'], inc, skp)):
            ti, ts = types(c, s, ji), types(c, s, js)
            nmd_i = 'nonsense_mediated_decay' in ti and 'protein_coding' not in ti and 'protein_coding' in ts
            nmd_s = 'nonsense_mediated_decay' in ts and 'protein_coding' not in ts and 'protein_coding' in ti
            if nmd_i or nmd_s:
                out.append({'junction_type': st, 'junction_id': str(m['ID'].iloc[i]),
                            'poison': 'inclusion' if nmd_i else 'skipping'})
    return pd.DataFrame(out)


def plot(series, path):
    """ series: (label, junctions with a change and a bin) """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.4), facecolor=cr.SURFACE)
    cr._style(ax)
    x = np.arange(len(LABELS))
    for (label, s), color in zip(series, cr.SERIES):
        rows = [summarise(s[s['bin'] == b]) if (s['bin'] == b).sum() > 1 else None for b in LABELS]
        mean = np.array([r['mean change'] if r else np.nan for r in rows])
        lo = np.array([r['95% from'] if r else np.nan for r in rows])
        hi = np.array([r['95% to'] if r else np.nan for r in rows])
        ax.plot(x, 100 * mean, color=color, linewidth=2, marker='o', markersize=4, label=label)
        ax.fill_between(x, 100 * lo, 100 * hi, color=color, alpha=0.15, linewidth=0)
    ax.axhline(0, color=cr.INK2, linewidth=0.8, linestyle=(0, (3, 3)))
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, fontsize=8)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('read count of the junction, all six samples summed', color=cr.INK, fontsize=9)
    ax.set_ylabel("mean change in the form's share, NMDI14 - DMSO (points)", color=cr.INK, fontsize=9)
    ax.set_title('Forms with a premature stop under the NMD inhibitor (95% intervals)', color=cr.INK, fontsize=11,
                 loc='left')
    ax.legend(frameon=False, fontsize=8, loc='upper left', labelcolor=cr.INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=cr.SURFACE)
    plt.close(fig)


def transcripts(path) -> pd.DataFrame:
    d = pd.read_csv(path)
    d = d[d['baseMean'] >= 10]
    rows = []
    for t in ['protein_coding', 'nonsense_mediated_decay', 'retained_intron']:
        s = d[d['transcript_type'] == t]
        lfc = s['log2FoldChange'].dropna()
        rows.append({'transcript_type': t, 'transcripts': len(s), 'median log2FC': round(lfc.median(), 4),
                     'share up': round((lfc > 0).mean(), 4),
                     'padj < 0.05, up': int(((s['padj'] < 0.05) & (s['log2FoldChange'] > 0)).sum()),
                     'padj < 0.05, down': int(((s['padj'] < 0.05) & (s['log2FoldChange'] < 0)).sum())})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--rmats', required=True)
    p.add_argument('--run', required=True)
    p.add_argument('--targets', default=None)
    p.add_argument('--dte', default=None)
    p.add_argument('--gtf', default=None)
    p.add_argument('--cache_dir', default=None)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    jx = pd.read_csv(Path(a.run) / 'junctions.tsv', sep='\t', dtype={'junction_id': str})
    jx = jx.merge(events(a.rmats), on=['junction_type', 'junction_id'], how='left')
    t = table(jx)
    res = pd.concat([by_bin(t, ['group', 'call']), by_bin(t.assign(call='both'), ['group', 'call'])],
                    ignore_index=True)
    res.to_csv(out / 'nmd.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(res.to_string(index=False))
    if a.targets:
        targets = set(pd.read_csv(a.targets)['gene_id'].astype(str).str.split('.').str[0])
        t4 = t[t['group'] == 'T4, either form']
        t4 = t4.assign(target=np.where(cr.unversioned(t4['gene_id']).isin(targets), 'gene with an NMD target',
                                       'other genes'))
        tg = by_bin(t4, ['target', 'call'])
        tg.to_csv(out / 'nmd_genes.tsv', sep='\t', index=False)
        with pd.option_context('display.width', 250, 'display.max_columns', 40):
            print(tg.to_string(index=False))
    if a.gtf:
        pf = poison_forms(a.rmats, a.gtf, a.cache_dir)
        g = jx.merge(pf, on=['junction_type', 'junction_id'])
        g = g.assign(change=np.where(g['poison'] == 'inclusion', 1, -1) * (g['psi_nmdi14'] - g['psi_dmso']),
                     group=np.where(g['fate'] == 'T4', 'JCAST T4', 'JCAST other fates'))
        g = g[g['change'].notna()]
        g['bin'] = pd.cut(g['min_jc'], BINS, right=False, labels=LABELS)
        pg = pd.concat([by_bin(g.assign(group='all'), ['group', 'call']), by_bin(g, ['group', 'call'])],
                       ignore_index=True)
        pg.to_csv(out / 'nmd_gencode.tsv', sep='\t', index=False)
        with pd.option_context('display.width', 250, 'display.max_columns', 40):
            print(pg.to_string(index=False))
    if a.dte:
        tt = transcripts(a.dte)
        tt.to_csv(out / 'nmd_transcripts.tsv', sep='\t', index=False)
        print(tt.to_string(index=False))
    series = [('JCAST T4: the form with the stop', t[t['group'] == 'T4, either form'])]
    if a.gtf:
        series.append(("GENCODE: the form only in NMD transcripts (control)", g))
    series.append(('JCAST T1: the inclusion form (baseline)', t[t['group'] == 'T1 (the inclusion form)']))
    plot(series, out / 'nmd.png')


if __name__ == '__main__':
    main()
