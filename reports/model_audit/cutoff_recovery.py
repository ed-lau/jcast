# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md): which known protein-coding isoforms a read-count
cutoff loses, measured on one run at -r 1.

A run writes a sequence once per file, under the first junction that gives it, so its output alone cannot say what a
higher cutoff would keep. The run is therefore made with an audit-only switch (audit_records.diff, in this folder):
JCAST_AUDIT_RECORDS=<file> makes main._finish write a line for every record of every junction, written or not, with
the junction's read count and the SHA-1 of the sequence. A sequence is in the output at cutoff c if any junction with
a read count >= c gives it. --check compares that with full runs at their own cutoffs.

The known isoforms (each a set of distinct sequences):
    Swiss-Prot      the isoforms (accession -N) of the --uniprot FASTA, less any sequence that is also a canonical
    GENCODE         the translations (--gencode, gencode.vNN.pc_translations.fa.gz) of protein_coding transcripts,
                    less any sequence of an Ensembl_canonical transcript or a Swiss-Prot canonical; with subsets
                    MANE (MANE_Select or MANE_Plus_Clinical, not canonical) and TSL <= 2
    GENCODE NMD     translations of nonsense_mediated_decay transcripts (a reference for T4 records)

Also, by read count, from the rMATS tables: whether the alternative region's length is a multiple of 3 (both forms
keep the frame), for the events JCAST translated (fate T1, T1_alt_terminus, T2 or T4); and whether rMATS calls the
event novel (in fromGTF.novelJunction.<type>.txt or fromGTF.novelSpliceSite.<type>.txt).

With --salmon, the expression of each GENCODE isoform: the TPM of the transcripts that encode its sequence, summed,
each averaged over the samples (<salmon>/<sample>/quant.sf). recovery.tsv then also counts the GENCODE isoforms with
at least 1 and 5 TPM, and expression.tsv gives, for the GENCODE isoforms reproduced at -r 1, their TPM by the read
count of the best junction giving them.

recovery_novel_only.tsv is recovery.tsv with the cutoff applied only to the events rMATS calls novel: an annotated
event is kept at any read count. tradeoff.tsv compares the two rules at equal database size: each row is the
novel-only rule at one cutoff, against the current rule interpolated (linearly, in database records) to the same size.

With --paxdb (a PaxDb integrated abundance file for the tissue), genes.tsv: per cutoff, the genes JCAST gives an
alternative sequence (a junction with a read count >= the cutoff and fate T1, T1_alt_terminus, T2 or T4), of all genes
and of the proteins PaxDb finds in the tissue, by abundance quartile; and, for those proteins, how many have a junction
with at least that many reads in either form (the gene's best junction: the highest, over its rMATS events, of the
larger of the summed inclusion and skipping counts). That is a positive control: a protein PaxDb finds is in the
tissue, so its best junction's read count is what a protein known to be there can have. paxdb.tsv summarises it by
quartile. paxdb_by_count.tsv turns it round: of the protein-coding genes (a GENCODE translation) with an rMATS event,
by the read count of the best junction, the share whose protein PaxDb finds. A PaxDb protein is matched to its gene
by the Ensembl protein of its STRING id (unversioned, in the --gencode translations), else by the gene name. Only the
proteins whose gene has an rMATS event are counted. With -c, JCAST writes every such gene's canonical whatever the
cutoff, so the canonicals themselves do not change with it.

Writes, into --out: recovery.tsv, recovery_novel_only.tsv (per cutoff) and tradeoff.tsv; bins.tsv (per read-count bin
and splice type); expression.tsv (with --salmon); genes.tsv, paxdb.tsv, paxdb_by_count.tsv and paxdb.png (with
--paxdb); summary.txt (with the --check results); recovery.png, bins.png, novel.png and tradeoff.png. The -m cutoff
(--model_cutoff) is added to the cutoffs if it is not one of them.

Usage (from the repository root; the whole-genome inputs are not bundled):
    python reports/model_audit/cutoff_recovery.py --records <records.tsv> --run <-r 1 output dir> --rmats <dir>
        --uniprot <Swiss-Prot isoform FASTA> --gencode <pc_translations.fa.gz> --gtf <gtf> [--cache_dir <dir>]
        [--salmon <dir>] [--paxdb <PaxDb file>] [--check <output dir>:<cutoff> ...] [--model_cutoff <n>] --out <dir>
"""

import argparse
import csv
import gzip
import hashlib
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rmats_versions import summed  # noqa: E402

CUTOFFS = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]
TICKS, TICKS_SPARSE = CUTOFFS[::2], CUTOFFS[::3]     # main adds the -m cutoff to CUTOFFS; the axes keep these
BINS = [1, 2, 3, 5, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, np.inf]   # [lo, next lo)
TRANSLATED = ['T1', 'T1_alt_terminus', 'T2', 'T4']
MERGE_TIERS = ['T1', 'T1_alt_cterm', 'T1_alt_nterm', 'T1_alt_both']      # uniprot.MERGE_TIERS
GROUPS = {'T1': MERGE_TIERS, 'T2': ['T2'], 'T4': ['T4'],
          'orphans': ['T1_orphan', 'T1_alt_cterm_orphan', 'T1_alt_nterm_orphan', 'T1_alt_both_orphan',
                      'T2_orphan', 'T4_orphan']}
SPLICE_TYPES = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']

# the reference palette's light steps (dataviz skill, references/palette.md)
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
QUARTILES = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']    # ordinal blue, steps 250 to 700
LABEL_BOX = dict(facecolor=SURFACE, edgecolor='none', pad=1)   # a cutoff's label, over its dashed line


def sha1(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()


def read_records(path) -> pd.DataFrame:
    cols = ['splice_type', 'junction_id', 'gene_id', 'min_jc', 'suffix', 'sha1', 'length', 'written', 'header']
    rec = pd.read_csv(path, sep='\t', names=cols, header=None, quoting=csv.QUOTE_NONE,
                      dtype={'junction_id': str, 'min_jc': int, 'written': int})
    return rec[rec['suffix'] != 'canonical']


def swissprot(path):
    """ (isoform sequences less canonicals, canonical sequences, every Swiss-Prot sequence) as SHA-1 sets """
    iso, can = set(), set()
    for r in SeqIO.parse(path, 'fasta'):
        (iso if '-' in r.id.split('|')[1] else can).add(sha1(str(r.seq)))
    return iso - can, can, iso | can


def gencode(path, gtf, cache_dir, sp_canonical):
    """ the GENCODE sets (see the module's docstring), and a table of the translations """
    from jcast import cache
    from jcast.annots import ReadAnnotations
    cache.set_directory(cache_dir)
    annot = ReadAnnotations(logging.getLogger('jcast-audit'), gtf).annot
    tx = annot.loc[annot['feature'] == 'transcript',
                   ['transcript_id', 'transcript_type', 'transcript_support_level', 'tag']].set_index('transcript_id')
    del annot
    rows = []
    with gzip.open(path, 'rt') as fh:
        for r in SeqIO.parse(fh, 'fasta'):
            f = r.id.split('|')
            rows.append((f[0], f[1], f[2], f[6], sha1(str(r.seq).rstrip('*'))))
    t = pd.DataFrame(rows, columns=['protein_id', 'transcript_id', 'gene_id', 'gene_name', 'sha1']
                     ).join(tx, on='transcript_id')
    tags = t['tag'].fillna('').str.split(',').map(set)
    canonical = set(t.loc[tags.map(lambda s: 'Ensembl_canonical' in s), 'sha1']) | sp_canonical
    pc = t['transcript_type'] == 'protein_coding'
    sets = {
        'GENCODE': set(t.loc[pc, 'sha1']) - canonical,
        'GENCODE MANE': set(t.loc[pc & tags.map(lambda s: bool(s & {'MANE_Select', 'MANE_Plus_Clinical'})), 'sha1'])
        - canonical,
        'GENCODE TSL<=2': set(t.loc[pc & (t['transcript_support_level'] <= 2), 'sha1']) - canonical,
        'GENCODE NMD': set(t.loc[t['transcript_type'] == 'nonsense_mediated_decay', 'sha1']) - canonical,
    }
    return sets, t


def expression(salmon_dir, translations, seqs) -> pd.Series:
    """ per sequence of seqs: the summed TPM of the transcripts encoding it, each averaged over the samples """
    quants = sorted(Path(salmon_dir).glob('*/quant.sf'))
    tpm = sum(pd.read_csv(q, sep='\t', usecols=['Name', 'TPM'], index_col='Name')['TPM'] for q in quants) / len(quants)
    t = translations[translations['sha1'].isin(seqs)]
    return t.assign(tpm=t['transcript_id'].map(tpm).fillna(0)).groupby('sha1')['tpm'].sum()


def frame_kept(rmats_dir) -> pd.DataFrame:
    """ per rMATS event: whether the alternative region's length is a multiple of 3; whether rMATS calls it novel; and
    the read count of its more used form (the larger of the summed inclusion and skipping junction counts) """
    out = []
    for st in SPLICE_TYPES:
        m = pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t')
        if st == 'SE':
            d = m['exonEnd'] - m['exonStart_0base']
        elif st == 'MXE':
            d = (m['1stExonEnd'] - m['1stExonStart_0base']) - (m['2ndExonEnd'] - m['2ndExonStart_0base'])
        elif st in ('A5SS', 'A3SS'):
            d = (m['longExonEnd'] - m['longExonStart_0base']) - (m['shortEE'] - m['shortES'])
        else:
            d = m['downstreamES'] - m['upstreamEE']
        novel = set()
        for kind in ['novelJunction', 'novelSpliceSite']:
            f = Path(rmats_dir) / f'fromGTF.{kind}.{st}.txt'
            if f.exists():
                novel |= set(pd.read_csv(f, sep='\t', usecols=['ID'])['ID'].astype(str))
        ids = m['ID'].astype(str)
        major = np.maximum(summed(m['IJC_SAMPLE_1']) + summed(m['IJC_SAMPLE_2']),
                           summed(m['SJC_SAMPLE_1']) + summed(m['SJC_SAMPLE_2']))
        out.append(pd.DataFrame({'junction_type': st, 'junction_id': ids, 'frame_kept': d % 3 == 0,
                                 'novel': ids.isin(novel), 'major_jc': major.astype(int)}))
    return pd.concat(out)


def max_count(rec, suffixes=None) -> pd.Series:
    """ per distinct sequence (within the files), the highest read count of a junction giving it """
    r = rec if suffixes is None else rec[rec['suffix'].isin(suffixes)]
    return r.groupby('sha1')['min_jc'].max()


def recovery(rec, jx, known, sp_all, tpm=None) -> pd.DataFrame:
    best = max_count(rec)
    if tpm is not None:
        known = {**known, 'GENCODE TPM>=1': set(tpm[tpm >= 1].index), 'GENCODE TPM>=5': set(tpm[tpm >= 5].index)}
    by_group = {g: max_count(rec, s) for g, s in GROUPS.items()}
    sp_rec = rec[rec['suffix'].isin(MERGE_TIERS) & rec['header'].str.startswith('sp|') & ~rec['sha1'].isin(sp_all)]
    merged = sp_rec.groupby('sha1')['min_jc'].max()
    rows = []
    for c in CUTOFFS:
        kept = jx[jx['min_jc'] >= c]
        row = {'cutoff': c, 'junctions': len(kept)}
        for f in TRANSLATED:
            row[f'{f} junctions'] = int((kept['fate'] == f).sum())
        for g, s in by_group.items():
            row[f'{g} records'] = int((s >= c).sum())
        row['merged JCAST records'] = int((merged >= c).sum())
        row['merged, not GENCODE'] = int((merged.drop(list(known['GENCODE']), errors='ignore') >= c).sum())
        for name, seqs in known.items():
            b = best.reindex(list(seqs)).dropna()
            row[name] = int((b >= c).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def equal_size(rc, rn, cutoffs=(3, 7, 20, 40, 100)) -> pd.DataFrame:
    """ the novel-only rule at each cutoff, and the current rule interpolated to the same number of database records """
    now = rc.sort_values('merged JCAST records')
    rows = []
    for c in cutoffs:
        r = rn[rn['cutoff'] == c].iloc[0]
        row = {'novel-only cutoff': c, 'database records': int(r['merged JCAST records'])}
        for col in ['Swiss-Prot', 'GENCODE TPM>=1', 'merged, not GENCODE']:
            if col in rn:
                row[f'{col}: novel-only'] = int(r[col])
                row[f'{col}: now'] = int(round(np.interp(r['merged JCAST records'], now['merged JCAST records'],
                                                         now[col])))
        rows.append(row)
    return pd.DataFrame(rows)


def novel_only(rec, jx):
    """ the records and junctions with the cutoff applied to novel events only: an annotated event is always kept """
    novel = jx.set_index(['junction_type', 'junction_id'])['novel']
    key = pd.MultiIndex.from_arrays([rec['splice_type'], rec['junction_id']])
    always = 10 ** 9
    rec = rec.assign(min_jc=np.where(novel.reindex(key).fillna(True).to_numpy(dtype=bool), rec['min_jc'], always))
    jx = jx.assign(min_jc=np.where(jx['novel'].fillna(True).to_numpy(dtype=bool), jx['min_jc'], always))
    return rec, jx


def bins(rec, jx, known) -> pd.DataFrame:
    """ per read-count bin and splice type, for the translated junctions: frame kept, and a known isoform given """
    union = known['Swiss-Prot'] | known['GENCODE']
    hit = rec[rec['sha1'].isin(union)].groupby(['splice_type', 'junction_id']).size()
    t = jx[jx['fate'].isin(TRANSLATED)].copy()
    t['known'] = pd.MultiIndex.from_arrays([t['junction_type'], t['junction_id']]).isin(hit.index)
    t['bin'] = pd.cut(t['min_jc'], BINS, right=False, labels=[f'{a:g}' for a in BINS[:-1]])
    rows = []
    for st in ['all'] + SPLICE_TYPES:
        s = t if st == 'all' else t[t['junction_type'] == st]
        for b, g in s.groupby('bin', observed=True):
            rows.append({'splice_type': st, 'bin': b, 'junctions': len(g),
                         'T1 fraction': round((g['fate'].isin(['T1', 'T1_alt_terminus'])).mean(), 4),
                         'frame kept fraction': round(g['frame_kept'].mean(), 4),
                         'known isoform fraction': round(g['known'].mean(), 4),
                         'novel event fraction': round(g['novel'].mean(), 4),
                         'known fraction, annotated events': round(g.loc[~g['novel'], 'known'].mean(), 4),
                         'known fraction, novel events': round(g.loc[g['novel'], 'known'].mean(), 4),
                         'frame kept, annotated events': round(g.loc[~g['novel'], 'frame_kept'].mean(), 4),
                         'frame kept, novel events': round(g.loc[g['novel'], 'frame_kept'].mean(), 4)})
    return pd.DataFrame(rows)


def expression_by_count(rec, tpm) -> pd.DataFrame:
    """ the GENCODE isoforms reproduced at -r 1, by the read count of the best junction giving them """
    best = max_count(rec).reindex(tpm.index).dropna()
    e = pd.DataFrame({'best': best, 'tpm': tpm.reindex(best.index)})
    e['bin'] = pd.cut(e['best'], BINS, right=False, labels=[f'{a:g}' for a in BINS[:-1]])
    return e.groupby('bin', observed=True)['tpm'].agg(
        isoforms='size', median_tpm='median', tpm_ge_1=lambda s: round((s >= 1).mean(), 4),
        tpm_ge_5=lambda s: round((s >= 5).mean(), 4)).reset_index()


def unversioned(ids: pd.Series) -> pd.Series:
    return ids.astype(str).str.split('.').str[0]


def read_paxdb(path, translations) -> pd.DataFrame:
    """ PaxDb's proteins with an abundance above 0, each with its gene (see the module's docstring) and quartile """
    pax = pd.read_csv(path, sep='\t', comment='#', names=['gene_name', 'string_id', 'abundance'])
    pax = pax[pax['abundance'] > 0].copy()
    by_protein = translations.assign(p=unversioned(translations['protein_id']), g=unversioned(translations['gene_id']))
    by_name = by_protein.groupby('gene_name')['g'].unique()
    by_name = by_name[by_name.map(len) == 1].map(lambda g: g[0])
    pax['gene'] = pax['string_id'].str.split('.', n=1).str[1].map(by_protein.drop_duplicates('p').set_index('p')['g'])
    pax['gene'] = pax['gene'].fillna(pax['gene_name'].map(by_name))
    pax['quartile'] = pd.qcut(pax['abundance'].rank(method='first'), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
    return pax


def genes(jx, pax, translations) -> tuple:
    """ per cutoff: the genes given an alternative sequence, and the PaxDb proteins, by quartile; by the best
    junction's read count, the protein-coding genes PaxDb finds (module docstring) """
    g = jx.assign(gene=unversioned(jx['gene_id']))
    best_min = g[g['fate'].isin(TRANSLATED)].groupby('gene')['min_jc'].max()
    best_major = g.groupby('gene')['major_jc'].max()
    p = pax[pax['gene'].isin(best_major.index)].assign(
        best_major=lambda d: d['gene'].map(best_major), best_min=lambda d: d['gene'].map(best_min))
    rows = []
    for c in CUTOFFS:
        row = {'cutoff': c, 'genes given an alternative sequence': int((best_min >= c).sum())}
        for q in ['all', 'Q1', 'Q2', 'Q3', 'Q4']:
            s = p if q == 'all' else p[p['quartile'] == q]
            row[f'PaxDb {q}: best junction >= cutoff'] = int((s['best_major'] >= c).sum())
            row[f'PaxDb {q}: given an alternative sequence'] = int((s['best_min'] >= c).sum())
        rows.append(row)
    summary = []
    for q in ['all', 'Q1', 'Q2', 'Q3', 'Q4']:
        s, sp = (pax, p) if q == 'all' else (pax[pax['quartile'] == q], p[p['quartile'] == q])
        summary.append({'quartile': q, 'proteins': len(s), 'abundance from (ppm)': s['abundance'].min(),
                        'gene matched': int(s['gene'].notna().sum()), 'gene has an rMATS event': len(sp),
                        'best junction, 10th percentile': sp['best_major'].quantile(0.1),
                        'best junction, median': sp['best_major'].median(),
                        **{f'best junction < {c}': round((sp['best_major'] < c).mean(), 4) for c in (7, 20, 40)},
                        'given an alternative sequence at -r 1': round(sp['best_min'].notna().mean(), 4)})
    coding = set(unversioned(translations.loc[translations['transcript_type'] == 'protein_coding', 'gene_id']))
    d = best_major[best_major.index.isin(coding)].to_frame('best_major')
    d['in PaxDb'] = d.index.isin(set(pax['gene'].dropna()))
    d['bin'] = pd.cut(d['best_major'], BINS, right=False, labels=[f'{a:g}' for a in BINS[:-1]])
    by_count = d.groupby('bin', observed=True)['in PaxDb'].agg(genes='size', in_paxdb='sum', share='mean')
    by_count['share'] = by_count['share'].round(4)
    return pd.DataFrame(rows), pd.DataFrame(summary), by_count.reset_index()


def check(rec, jx, known, out_dir, cutoff, sp_iso) -> list:
    """ a full run at its own cutoff, against what the -r 1 run's records give at that cutoff """
    lines = [f'== {out_dir} at {cutoff}']
    other = pd.read_csv(Path(out_dir) / 'junctions.tsv', sep='\t', dtype={'junction_id': str})
    a = jx[jx['min_jc'] >= cutoff].set_index(['junction_type', 'junction_id'])['fate']
    b = other[other['fate'] != 'skipped_low_count'].set_index(['junction_type', 'junction_id'])['fate']
    same = a.index.equals(b.index) or set(a.index) == set(b.index)
    differ = (a.reindex(b.index) != b).sum() if same else None
    lines.append(f'junctions kept: {len(a)} derived, {len(b)} in the run; same junctions: {same}; '
                 f'fates differing: {differ}')
    for f in sorted(Path(out_dir).glob('psq_*.fasta')):
        suffix = f.stem[4:]
        if suffix == 'canonical':
            continue
        run = {sha1(str(r.seq)) for r in SeqIO.parse(f, 'fasta')}
        r = rec[rec['suffix'] == suffix]
        derived = set(r.loc[r['min_jc'] >= cutoff, 'sha1'])
        lines.append(f'{f.name}: {len(derived)} derived, {len(run)} in the run, equal: {derived == run}')
    run_iso = set()
    for f in Path(out_dir).glob('psq_*.fasta'):
        if f.stem != 'psq_canonical':
            run_iso |= {sha1(str(r.seq)) for r in SeqIO.parse(f, 'fasta')} & sp_iso
    best = max_count(rec).reindex(list(sp_iso)).dropna()
    derived_iso = set(best[best >= cutoff].index)
    lines.append(f'Swiss-Prot isoform sequences: {len(derived_iso)} derived, {len(run_iso)} in the run, '
                 f'equal: {derived_iso == run_iso}')
    return lines


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_recovery(rc, path, marks):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    series = [('Swiss-Prot', 'Swiss-Prot isoforms'), ('GENCODE', 'GENCODE alternative proteins'),
              ('GENCODE TPM>=1', 'GENCODE alternative proteins, at least 1 TPM'),
              ('merged, not GENCODE', 'database records in neither (unannotated)')]
    series = [s for s in series if s[0] in rc]
    fig, ax = plt.subplots(figsize=(7.5, 4.4), facecolor=SURFACE)
    _style(ax)
    for (col, label), color in zip(series, SERIES):
        y = 100 * rc[col] / rc[col].iloc[0]
        ax.plot(rc['cutoff'], y, color=color, linewidth=2, marker='o', markersize=4,
                label=f'{label} ({rc[col].iloc[0]:,} at -r 1)')
    for x, text in marks:
        ax.axvline(x, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)))
        ax.text(x, 101, text, fontsize=8, color=INK2, ha='center', va='bottom', bbox=LABEL_BOX)
    ax.set_xscale('log')
    ax.set_xticks(TICKS)
    ax.set_xticklabels([str(c) for c in TICKS])
    ax.set_ylim(0, 108)
    ax.set_xlabel('read-count cutoff (-r): junctions kept have at least this many reads', color=INK, fontsize=9)
    ax.set_ylabel('% of the value at -r 1', color=INK, fontsize=9)
    ax.set_title('Sequences reproduced at each cutoff, relative to no cutoff (-r 1)', color=INK, fontsize=11,
                 loc='left')
    ax.legend(frameon=False, fontsize=8.5, loc='lower left', labelcolor=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_tradeoff(rc, rn, path, marks):
    """ Swiss-Prot isoforms reproduced against database records, for the two rules """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.6), facecolor=SURFACE)
    _style(ax)
    for r, label, color, offset in [(rc, 'cutoff on every junction (now)', SERIES[0], (4, -13)),
                                    (rn, 'cutoff on novel events only', SERIES[1], (-4, 7))]:
        ax.plot(r['merged JCAST records'], r['Swiss-Prot'], color=color, linewidth=2, marker='o', markersize=4,
                label=label)
        for c in marks:
            row = r[r['cutoff'] == c].iloc[0]
            ax.annotate(f'{c}', (row['merged JCAST records'], row['Swiss-Prot']), xytext=offset,
                        textcoords='offset points', fontsize=8, color=INK2, ha='left' if offset[0] > 0 else 'right')
    ax.set_xlabel('JCAST records added to the search database (labels: the cutoff in reads)', color=INK, fontsize=9)
    ax.set_ylabel('Swiss-Prot isoforms reproduced', color=INK, fontsize=9)
    ax.set_title('Known isoforms kept for the database size, two rules', color=INK, fontsize=11, loc='left')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax.legend(frameon=False, fontsize=8.5, loc='lower right', labelcolor=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_paxdb(bc, gn, ps, path, marks):
    """ left: the protein-coding genes PaxDb finds, by the best junction's read count; right: the PaxDb proteins with
    an rMATS event given an alternative sequence at each cutoff, by abundance quartile """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.4), facecolor=SURFACE, sharey=True)
    _style(left)
    x = np.arange(len(bc))
    left.plot(x, 100 * bc['share'], color=SERIES[0], linewidth=2, marker='o', markersize=4)
    left.set_xticks(x)
    left.set_xticklabels(list(bc['bin']), fontsize=8)
    left.set_xlabel("read count of the gene's best junction, either form (bin starting at)", color=INK, fontsize=9)
    left.set_ylabel('%', color=INK, fontsize=9)
    left.set_title('Protein-coding genes whose protein PaxDb finds', color=INK, fontsize=10, loc='left')
    _style(right)
    n = ps.set_index('quartile')['gene has an rMATS event']
    for q, color, label in zip(['Q1', 'Q2', 'Q3', 'Q4'], QUARTILES,
                               ['least abundant quarter', 'second', 'third', 'most abundant quarter']):
        right.plot(gn['cutoff'], 100 * gn[f'PaxDb {q}: given an alternative sequence'] / n[q], color=color,
                   linewidth=2, marker='o', markersize=4, label=f'{label} ({n[q]:,})')
    for xc, text in marks:
        right.axvline(xc, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)))
        right.text(xc, 101, text, fontsize=8, color=INK2, ha='center', va='bottom', bbox=LABEL_BOX)
    right.set_xscale('log')
    right.set_xticks(TICKS_SPARSE)
    right.set_xticklabels([str(c) for c in TICKS_SPARSE])
    right.set_xlabel('read-count cutoff (-r)', color=INK, fontsize=9)
    right.set_title('PaxDb proteins given an alternative sequence', color=INK, fontsize=10, loc='left')
    right.legend(frameon=False, fontsize=8.5, loc='lower left', labelcolor=INK, title='PaxDb abundance',
                 title_fontsize=8.5)
    left.set_ylim(0, 108)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_novel(bn, path):
    """ by read count, annotated against novel events: the frame kept (left) and a known isoform given (right) """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a = bn[bn['splice_type'] == 'all']
    x = np.arange(len(a))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), facecolor=SURFACE, sharey=True)
    for ax, measure, title in [(axes[0], 'frame kept', 'Alternative region a multiple of 3 nt'),
                               (axes[1], 'known fraction', 'Gives a known protein-coding isoform')]:
        _style(ax)
        for group, color in [('annotated', SERIES[0]), ('novel', SERIES[1])]:
            ax.plot(x, 100 * a[f'{measure}, {group} events'], color=color, linewidth=2, marker='o', markersize=4,
                    label=f'{group} events')
        ax.set_xticks(x)
        ax.set_xticklabels(list(a['bin']), fontsize=8)
        ax.set_xlabel('read count of the junction (bin starting at)', color=INK, fontsize=9)
        ax.set_title(title, color=INK, fontsize=10, loc='left')
    axes[0].axhline(100 / 3, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)), label='one in three: a random length')
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel('% of translated junctions in the bin', color=INK, fontsize=9)
    axes[0].legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_bins(bn, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a = bn[bn['splice_type'] == 'all']
    x = np.arange(len(a))
    fig, ax = plt.subplots(figsize=(7.5, 4.4), facecolor=SURFACE)
    _style(ax)
    for (col, label), color in zip([('known isoform fraction', 'gives a known protein-coding isoform'),
                                    ('frame kept fraction', 'alternative region a multiple of 3 nt'),
                                    ('T1 fraction', 'translated as T1 (no frameshift, no stop)')], SERIES):
        ax.plot(x, 100 * a[col], color=color, linewidth=2, marker='o', markersize=4, label=label)
    ax.axhline(100 / 3, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)), label='one in three: a random length')
    ax.set_xticks(x)
    ax.set_xticklabels(list(a['bin']), fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_xlabel('read count of the junction (bin starting at)', color=INK, fontsize=9)
    ax.set_ylabel('% of translated junctions in the bin', color=INK, fontsize=9)
    ax.set_title('Translated junctions by read count', color=INK, fontsize=11, loc='left')
    ax.legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--records', required=True)
    p.add_argument('--run', required=True, help='the -r 1 output directory the records came from')
    p.add_argument('--rmats', required=True)
    p.add_argument('--uniprot', required=True)
    p.add_argument('--gencode', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--cache_dir', default=None)
    p.add_argument('--salmon', default=None, help='a folder of Salmon outputs, <sample>/quant.sf')
    p.add_argument('--paxdb', default=None, help='a PaxDb integrated abundance file for the tissue')
    p.add_argument('--check', action='append', default=[], help='<output dir>:<cutoff> of a full run')
    p.add_argument('--model_cutoff', type=int, default=40, help='the -m cutoff, marked on the plot')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if a.model_cutoff not in CUTOFFS:
        CUTOFFS.append(a.model_cutoff)
        CUTOFFS.sort()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    rec = read_records(a.records)
    jx = pd.read_csv(Path(a.run) / 'junctions.tsv', sep='\t', dtype={'junction_id': str})
    jx = jx.merge(frame_kept(a.rmats), on=['junction_type', 'junction_id'], how='left')
    sp_iso, sp_can, sp_all = swissprot(a.uniprot)
    gc, translations = gencode(a.gencode, a.gtf, a.cache_dir, sp_can)
    known = {'Swiss-Prot': sp_iso, **gc}
    known['Swiss-Prot or GENCODE'] = sp_iso | gc['GENCODE']
    tpm = expression(a.salmon, translations, gc['GENCODE']) if a.salmon else None

    rc = recovery(rec, jx, known, sp_all, tpm)
    rc.to_csv(out / 'recovery.tsv', sep='\t', index=False)
    rn = recovery(*novel_only(rec, jx), known, sp_all, tpm)
    rn.to_csv(out / 'recovery_novel_only.tsv', sep='\t', index=False)
    equal_size(rc, rn).to_csv(out / 'tradeoff.tsv', sep='\t', index=False)
    bn = bins(rec, jx, known)
    bn.to_csv(out / 'bins.tsv', sep='\t', index=False)
    if tpm is not None:
        ex = expression_by_count(rec, tpm)
        ex.to_csv(out / 'expression.tsv', sep='\t', index=False)
    if a.paxdb:
        gn, ps, bc = genes(jx, read_paxdb(a.paxdb, translations), translations)
        gn.to_csv(out / 'genes.tsv', sep='\t', index=False)
        ps.to_csv(out / 'paxdb.tsv', sep='\t', index=False)
        bc.to_csv(out / 'paxdb_by_count.tsv', sep='\t', index=False)

    lines = [f'records: {len(rec)} (not canonical), {rec["sha1"].nunique()} distinct sequences; '
             f'junctions.tsv: {len(jx)} rows']
    lines += [f'known set {k}: {len(v)} sequences' for k, v in known.items()]
    for c in a.check:
        d, cutoff = c.rsplit(':', 1)
        lines += check(rec, jx, known, d, int(cutoff), sp_iso)
    (out / 'summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(rc.to_string(index=False))
        print('cutoff on novel events only:')
        print(rn.to_string(index=False))
        print(bn[bn['splice_type'] == 'all'].to_string(index=False))
        if tpm is not None:
            print(ex.to_string(index=False))
        if a.paxdb:
            print(gn.to_string(index=False))
            print(ps.to_string(index=False))
            print(bc.to_string(index=False))

    plot_recovery(rc, out / 'recovery.png', [(7, '-r 7'), (a.model_cutoff, f'-m ({a.model_cutoff})')])
    plot_bins(bn, out / 'bins.png')
    plot_novel(bn, out / 'novel.png')
    plot_tradeoff(rc, rn, out / 'tradeoff.png', [1, 3, 7, 20, a.model_cutoff])
    if a.paxdb:
        plot_paxdb(bc, gn, ps, out / 'paxdb.png', [(7, '-r 7'), (a.model_cutoff, f'-m ({a.model_cutoff})')])


if __name__ == '__main__':
    main()
