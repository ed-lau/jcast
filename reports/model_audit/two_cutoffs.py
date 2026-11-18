# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §9): a read-count cutoff for annotated events and another
for novel events, from the records of one run at -r 1 (cutoff_recovery.py's inputs; the same known sets).

An event is novel if rMATS lists it in fromGTF.novelJunction.<type>.txt or fromGTF.novelSpliceSite.<type>.txt. At the
pair (a, n), a junction is kept if it is annotated with a read count >= a, or novel with a read count >= n; a sequence
is kept if a kept junction gives it. a = n is the current rule; a = 1 is the cutoff on novel events only (§§4.5, 5.6).

Writes, into --out:
    pairs.tsv       per pair: junctions kept, database records (the merged FASTA's JCAST records), Swiss-Prot isoforms
                    and GENCODE alternative proteins of at least 1 TPM reproduced, and records in neither
    frontier.tsv    at the database size of the current rule at each of a few cutoffs: the current rule, and the pair
                    that reproduces the most Swiss-Prot isoforms within that size
    gtf_calls.tsv   per splice type, rMATS's call against one made from the GTF alone: an event is annotated if each of
                    its splice junctions is an intron of a GTF transcript (and, for RI, a GTF exon also covers the
                    intron: the retained form has no junction), a novel junction if not but each splice site is one,
                    and a novel splice site otherwise. Whether JCAST could tell novel events itself, from any rMATS
                    version (with --gtf)
    pairs.png       Swiss-Prot isoforms (left) and expressed GENCODE isoforms (right) against database records

Usage (from the repository root):
    python reports/model_audit/two_cutoffs.py --records <records.tsv> --run <-r 1 output dir> --rmats <dir>
        --uniprot <FASTA> --gencode <pc_translations.fa.gz> --gtf <gtf> [--cache_dir <dir>] [--salmon <dir>]
        [--model_cutoff <n>] --out <dir>
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cutoff_recovery as cr  # noqa: E402

ANNOTATED_CUTOFFS = [1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 100]  # main adds every cutoff of the current rule
NOVEL_CUTOFFS = cr.CUTOFFS + [10 ** 9]          # the last: no novel event


def best_by_call(rec, novel) -> pd.DataFrame:
    """ per sequence: the highest read count of an annotated junction giving it, and of a novel one (0 if none) """
    key = pd.MultiIndex.from_arrays([rec['splice_type'], rec['junction_id']])
    is_novel = novel.reindex(key).fillna(True).to_numpy(dtype=bool)
    r = rec.assign(annotated=np.where(is_novel, 0, rec['min_jc']), novel=np.where(is_novel, rec['min_jc'], 0))
    return r.groupby('sha1')[['annotated', 'novel']].max()


def pairs(rec, jx, known, sp_all, tpm) -> pd.DataFrame:
    novel = jx.set_index(['junction_type', 'junction_id'])['novel']
    every = best_by_call(rec, novel)
    db = best_by_call(rec[rec['suffix'].isin(cr.MERGE_TIERS) & rec['header'].str.startswith('sp|')
                          & ~rec['sha1'].isin(sp_all)], novel)
    sets = {'Swiss-Prot': every.reindex(list(known['Swiss-Prot'])).dropna()}
    if tpm is not None:
        sets['GENCODE TPM>=1'] = every.reindex(list(tpm[tpm >= 1].index)).dropna()
    neither = db.drop(list(known['GENCODE']), errors='ignore')
    j = jx[jx['fate'].isin(cr.TRANSLATED + ['skipped_low_count'])]
    j_novel = j['novel'].fillna(True).to_numpy(dtype=bool)
    rows = []
    for a in ANNOTATED_CUTOFFS:
        for n in NOVEL_CUTOFFS:
            kept = lambda t: int(((t['annotated'] >= a) | (t['novel'] >= n)).sum())
            row = {'annotated cutoff': a, 'novel cutoff': n if n < 10 ** 9 else 'none kept',
                   'junctions kept': int(np.where(j_novel, j['min_jc'] >= n, j['min_jc'] >= a).sum()),
                   'database records': kept(db)}
            row.update({k: kept(v) for k, v in sets.items()})
            row['database records in neither'] = kept(neither)
            rows.append(row)
    return pd.DataFrame(rows)


def frontier(pt, targets) -> pd.DataFrame:
    """ at the size of the current rule at each target cutoff: that rule, and the pairs within the size with the most
    Swiss-Prot isoforms (of all pairs, and of those keeping some novel events), and with the most expressed GENCODE
    isoforms """
    now = pt[pt['annotated cutoff'] == pt['novel cutoff']]
    rows = []
    for c in targets:
        cur = now[now['annotated cutoff'] == c]
        if cur.empty:
            continue
        cur = cur.iloc[0]
        within = pt[pt['database records'] <= cur['database records']]
        best = lambda d, col: d.sort_values([col, 'database records'], ascending=[False, True]).iloc[0]
        picks = [(f'now, -r {c}', cur), ('most Swiss-Prot isoforms', best(within, 'Swiss-Prot')),
                 ('the same, keeping novel events', best(within[within['novel cutoff'] != 'none kept'], 'Swiss-Prot'))]
        if 'GENCODE TPM>=1' in pt:
            picks.append(('most expressed GENCODE isoforms', best(within, 'GENCODE TPM>=1')))
        for rule, r in picks:
            rows.append({'size of': f'-r {c}', 'rule': rule, **r.to_dict()})
    return pd.DataFrame(rows)


# --- novel events, from the GTF -------------------------------------------------------------------------------------

def gtf_introns(gtf, cache_dir):
    """ the GTF's introns as (chr, strand, last base of the exon before, 0-based start of the exon after): rMATS's
    coordinates; and its splice sites, each side apart """
    from jcast import cache
    from jcast.annots import ReadAnnotations
    cache.set_directory(cache_dir)
    annot = ReadAnnotations(logging.getLogger('jcast-audit'), gtf).annot
    ex = annot.loc[annot['feature'] == 'exon', ['seqname', 'strand', 'transcript_id', 'start', 'end']]
    del annot
    ex = ex.sort_values(['transcript_id', 'start'])
    nxt = ex.groupby('transcript_id').shift(-1)
    has = nxt['start'].notna()
    introns = pd.DataFrame({'chr': ex.loc[has, 'seqname'].to_numpy(), 'strand': ex.loc[has, 'strand'].to_numpy(),
                            'left': ex.loc[has, 'end'].to_numpy(int),
                            'right': (nxt.loc[has, 'start'] - 1).to_numpy(int)}).drop_duplicates()
    key = lambda d, *cols: set(zip(*(d[c] for c in ('chr', 'strand') + cols)))
    spans = {}                                   # per chr and strand: exon starts, sorted, and the furthest end yet
    for (c, s), g in ex.groupby(['seqname', 'strand']):
        g = g.sort_values('start')
        spans[(c, s)] = (g['start'].to_numpy(), np.maximum.accumulate(g['end'].to_numpy()))
    return key(introns, 'left', 'right'), key(introns, 'left'), key(introns, 'right'), spans


def spanned(spans, c, s, first, last) -> bool:
    """ whether a GTF exon covers the 1-based bases first to last """
    if (c, s) not in spans:
        return False
    starts, ends = spans[(c, s)]
    i = np.searchsorted(starts, first, side='right') - 1
    return i >= 0 and ends[i] >= last


def event_junctions(m, st) -> list:
    """ per event, its splice junctions as (left, right) pairs in rMATS's coordinates """
    if st == 'SE':
        cols = [('upstreamEE', 'exonStart_0base'), ('exonEnd', 'downstreamES'), ('upstreamEE', 'downstreamES')]
        return [[(r[a], r[b]) for a, b in cols] for _, r in m.iterrows()]
    if st == 'MXE':
        cols = [('upstreamEE', '1stExonStart_0base'), ('1stExonEnd', 'downstreamES'),
                ('upstreamEE', '2ndExonStart_0base'), ('2ndExonEnd', 'downstreamES')]
        return [[(r[a], r[b]) for a, b in cols] for _, r in m.iterrows()]
    if st in ('A5SS', 'A3SS'):
        out = []
        for _, r in m.iterrows():
            if r['flankingES'] >= r['longExonEnd']:
                out.append([(r['longExonEnd'], r['flankingES']), (r['shortEE'], r['flankingES'])])
            else:
                out.append([(r['flankingEE'], r['longExonStart_0base']), (r['flankingEE'], r['shortES'])])
        return out
    return [[(r['upstreamEE'], r['downstreamES'])] for _, r in m.iterrows()]


def gtf_calls(rmats_dir, gtf, cache_dir) -> pd.DataFrame:
    introns, lefts, rights, spans = gtf_introns(gtf, cache_dir)
    rows = []
    for st in cr.SPLICE_TYPES:
        m = pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t')
        call = pd.Series('annotated', index=m.index)
        for kind in ['novelSpliceSite', 'novelJunction']:
            f = Path(rmats_dir) / f'fromGTF.{kind}.{st}.txt'
            if f.exists():
                ids = set(pd.read_csv(f, sep='\t', usecols=['ID'])['ID'])
                call[m['ID'].isin(ids)] = kind
        mine = []
        for (_, r), js in zip(m.iterrows(), event_junctions(m, st)):
            c, s = r['chr'], r['strand']
            retained = st != 'RI' or spanned(spans, c, s, r['upstreamEE'] + 1, r['downstreamES'])
            if retained and all((c, s, a, b) in introns for a, b in js):
                mine.append('annotated')
            elif all((c, s, a) in lefts and (c, s, b) in rights for a, b in js):
                mine.append('novelJunction')
            else:
                mine.append('novelSpliceSite')
        t = pd.crosstab(call.rename('rMATS'), pd.Series(mine, index=m.index, name='GTF'))
        for (a, b), v in t.stack().items():
            rows.append({'splice_type': st, 'rMATS': a, 'GTF': b, 'events': int(v)})
    return pd.DataFrame(rows)


# --- the chart ------------------------------------------------------------------------------------------------------

def plot_pairs(pt, path, measures, marks):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(measures), figsize=(5.4 * len(measures), 4.6), facecolor=SURFACE)
    axes = np.atleast_1d(axes)
    now = pt[pt['annotated cutoff'] == pt['novel cutoff']]
    numeric = pt[pt['novel cutoff'] != 'none kept']
    for ax, (col, title) in zip(axes, measures):
        cr._style(ax)
        ax.plot(now['database records'], now[col], color=cr.SERIES[1], linewidth=2, marker='o', markersize=4,
                label='one cutoff for every event (now)')
        for a, color in zip([1, 3, 7, 20], cr.QUARTILES):
            s = numeric[(numeric['annotated cutoff'] == a) & (numeric['novel cutoff'].astype(float) >= a)]
            ax.plot(s['database records'], s[col], color=color, linewidth=2, marker='o', markersize=3,
                    label=f'annotated >= {a}, novel cutoff varied' + (' (novel only)' if a == 1 else ''))
        for c, text in marks:
            r = now[now['annotated cutoff'] == c]
            if not r.empty:
                model = text.startswith('-m')            # above and left, clear of a -r label nearby
                ax.annotate(text, (r['database records'].iloc[0], r[col].iloc[0]),
                            xytext=(-6, 6) if model else (5, -12), ha='right' if model else 'left',
                            textcoords='offset points', fontsize=8, color=cr.INK2)
        ax.set_xlabel('JCAST records added to the search database', color=cr.INK, fontsize=9)
        ax.set_title(title, color=cr.INK, fontsize=10, loc='left')
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
    axes[0].legend(frameon=False, fontsize=8, loc='lower right', labelcolor=cr.INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


SURFACE = cr.SURFACE


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for arg in ['--records', '--run', '--rmats', '--uniprot', '--gencode', '--gtf', '--out']:
        p.add_argument(arg, required=True)
    p.add_argument('--cache_dir', default=None)
    p.add_argument('--salmon', default=None)
    p.add_argument('--model_cutoff', type=int, default=40)
    a = p.parse_args()
    if a.model_cutoff not in cr.CUTOFFS:
        cr.CUTOFFS.append(a.model_cutoff)
        cr.CUTOFFS.sort()
        NOVEL_CUTOFFS[:] = cr.CUTOFFS + [10 ** 9]
    ANNOTATED_CUTOFFS[:] = sorted(set(ANNOTATED_CUTOFFS) | set(cr.CUTOFFS))    # the diagonal: the current rule
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    rec = cr.read_records(a.records)
    jx = pd.read_csv(Path(a.run) / 'junctions.tsv', sep='\t', dtype={'junction_id': str})
    jx = jx.merge(cr.frame_kept(a.rmats), on=['junction_type', 'junction_id'], how='left')
    sp_iso, sp_can, sp_all = cr.swissprot(a.uniprot)
    gc, translations = cr.gencode(a.gencode, a.gtf, a.cache_dir, sp_can)
    known = {'Swiss-Prot': sp_iso, **gc}
    tpm = cr.expression(a.salmon, translations, gc['GENCODE']) if a.salmon else None

    pt = pairs(rec, jx, known, sp_all, tpm)
    pt.to_csv(out / 'pairs.tsv', sep='\t', index=False)
    fr = frontier(pt, [3, 7, 20, a.model_cutoff, 100])
    fr.to_csv(out / 'frontier.tsv', sep='\t', index=False)
    gcalls = gtf_calls(a.rmats, a.gtf, a.cache_dir)
    gcalls.to_csv(out / 'gtf_calls.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(fr.to_string(index=False))
        print(pt[pt['annotated cutoff'].isin([1, 2, 3, 5, 7])].to_string(index=False))
        print(gcalls.pivot_table(index=['splice_type', 'rMATS'], columns='GTF', values='events', fill_value=0))
    measures = [('Swiss-Prot', 'Swiss-Prot isoforms reproduced')]
    if tpm is not None:
        measures.append(('GENCODE TPM>=1', 'GENCODE alternative proteins, at least 1 TPM'))
    plot_pairs(pt, out / 'pairs.png', measures,
               [(3, '-r 3'), (7, '-r 7'), (20, '-r 20'), (a.model_cutoff, f'-m ({a.model_cutoff})')])


if __name__ == '__main__':
    main()
