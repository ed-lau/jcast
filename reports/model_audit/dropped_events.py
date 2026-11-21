# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §15): the events rMATS drops before JCAST reads them.

With its statistics on (the default), rMATS's post step keeps an event in JC.raw.input.<type>.txt, and so in
<type>.MATS.JC.txt, only if each group has a read, the event has an inclusion read and a skipping read, and both forms
have a length (rmats.py, filter_countfile, rMATS 4.4.0). An event whose reads all come from one group is dropped; it
is not given an NA p-value. With --statoff nothing is dropped, and every PValue and FDR is NA. This compares a folder
made with statistics against a post step with --statoff on the same prep files (the same events and IDs).

Without --rmats (§19, the heart folders, made with --statoff only), the events the statistics would keep are those
that pass filter_countfile's conditions, from the --statoff folder's counts and form lengths. With --rmats, the
conditions are checked against it: they keep exactly its events.

The read count is JCAST's (Junction.min_jc): the smaller of the inclusion and skipping counts, each summed over both
groups and their replicates.

Two steps:
    prepare     count the dropped events, and write those with a read in both forms as an rMATS folder (--dropped),
                FDR and PValue set to 0 so that -q keeps them, for a JCAST run at -r 1 with the audit switch
    measure     (with --records, --dropped_records and --dropped_run) what the dropped events would add to the run

Writes, into --out:
    dropped.tsv         per splice type: events rMATS found, kept and dropped; of those dropped, why; with a read in
                        both forms, at or above the -m cutoff, and novel (rMATS's call)
    dropped_bins.tsv    per read-count bin: events kept and dropped
    dropped_genes.tsv   per gene: dropped events with a read in both forms, the highest read count, the group with none
    dropped_model.txt   the -m cutoff with and without the dropped events, and on the --statoff folder as it is
    dropped_known.tsv   (measure) what the dropped events' records add that no kept event gives, at -r 1, at the
                        -m cutoff and at the cutoff fitted with the dropped events (a --statoff folder's counts of at
                        least 1): junctions by fate, database records, Swiss-Prot isoforms, GENCODE alternative
                        proteins (and those of at least 1 TPM, with --salmon)

Usage (from the repository root):
    python reports/model_audit/dropped_events.py --rmats <dir> --statoff <dir> --dropped <dir> --out <dir>
    python reports/model_audit/dropped_events.py --rmats <dir> --statoff <dir> --out <dir>
        --records <records.tsv> --dropped_records <records.tsv> --dropped_run <-r 1 output dir of --dropped>
        --uniprot <FASTA> --gencode <pc_translations.fa.gz> --gtf <gtf> [--cache_dir <dir>] [--salmon <dir>]
    python reports/model_audit/dropped_events.py --statoff <dir> --out <dir>
        --records <records.tsv of a -r 1 run on --statoff> --dropped_run <that run's output dir>
        --uniprot <FASTA> --gencode <pc_translations.fa.gz> --gtf <gtf> [--cache_dir <dir>] [--salmon <dir>]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cutoff_recovery as cr  # noqa: E402
from rmats_versions import summed  # noqa: E402

COUNTS = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2']
REASONS = ['no reads', 'one form only', 'no reads in group 1', 'no reads in group 2', 'other']


def events(rmats_dir, statoff_dir) -> pd.DataFrame:
    """ every event of the --statoff folder: its counts, whether rMATS kept it with statistics on (in --rmats, or
    without it by filter_countfile's conditions), whether those conditions keep it, why not, and whether rMATS calls
    it novel """
    out = []
    for st in cr.SPLICE_TYPES:
        m = pd.read_csv(Path(statoff_dir) / f'{st}.MATS.JC.txt', sep='\t',
                        usecols=['ID', 'GeneID', 'geneSymbol', *COUNTS, 'IncFormLen', 'SkipFormLen'])
        i1, s1, i2, s2 = (summed(m[c]) for c in COUNTS)
        passes = ((i1 + s1 > 0) & (i2 + s2 > 0) & (i1 + i2 > 0) & (s1 + s2 > 0) & (m['IncFormLen'] > 0)
                  & (m['SkipFormLen'] > 0))
        if rmats_dir is None:
            kept = set(m.loc[passes, 'ID'])
        else:
            kept = set(pd.read_csv(Path(rmats_dir) / f'{st}.MATS.JC.txt', sep='\t', usecols=['ID'])['ID'])
            missing = kept - set(m['ID'])
            if missing:
                raise ValueError(f'{st}: {len(missing)} events kept with statistics are not in the --statoff folder')
        novel = set()
        for kind in ['novelJunction', 'novelSpliceSite']:
            novel |= set(pd.read_csv(Path(statoff_dir) / f'fromGTF.{kind}.{st}.txt', sep='\t', usecols=['ID'])['ID'])
        e = pd.DataFrame({'splice_type': st, 'ID': m['ID'], 'gene_id': m['GeneID'], 'gene': m['geneSymbol'],
                          'group1': (i1 + s1).astype(int), 'group2': (i2 + s2).astype(int),
                          'min_jc': np.minimum(i1 + i2, s1 + s2).astype(int), 'kept': m['ID'].isin(kept),
                          'passes': passes.to_numpy(), 'novel': m['ID'].isin(novel)})
        e['reason'] = np.select([e['kept'], e['group1'] + e['group2'] == 0, e['min_jc'] == 0, e['group1'] == 0,
                                 e['group2'] == 0],
                                ['', 'no reads', 'one form only', 'no reads in group 1', 'no reads in group 2'],
                                default='other')
        out.append(e)
    return pd.concat(out, ignore_index=True)


def by_type(ev, cutoff) -> pd.DataFrame:
    rows = []
    for st in cr.SPLICE_TYPES + ['all']:
        e = ev if st == 'all' else ev[ev['splice_type'] == st]
        d = e[~e['kept']]
        read = d[d['min_jc'] >= 1]
        row = {'splice_type': st, 'events': len(e), 'kept': int(e['kept'].sum()), 'dropped': len(d)}
        row.update({r: int((d['reason'] == r).sum()) for r in REASONS})
        row.update({'dropped, read in both forms': len(read),
                    'the same, as a share of kept events with one': round(len(read) / int(
                        (e['kept'] & (e['min_jc'] >= 1)).sum()), 4),
                    f'dropped, >= {cutoff} reads': int((read['min_jc'] >= cutoff).sum()),
                    f'kept, >= {cutoff} reads': int((e['kept'] & (e['min_jc'] >= cutoff)).sum()),
                    'dropped with a read in both forms: novel': round(read['novel'].mean(), 4),
                    'kept: novel': round(e.loc[e['kept'], 'novel'].mean(), 4)})
        rows.append(row)
    return pd.DataFrame(rows)


def by_bin(ev) -> pd.DataFrame:
    e = ev[ev['min_jc'] >= 1].copy()
    e['bin'] = pd.cut(e['min_jc'], cr.BINS, right=False, labels=[f'{a:g}' for a in cr.BINS[:-1]])
    t = e.groupby('bin', observed=True).agg(kept=('kept', 'sum'), events=('kept', 'size'))
    t['dropped'] = t['events'] - t['kept']
    t['dropped share'] = (t['dropped'] / t['events']).round(4)
    return t.reset_index()


def by_gene(ev) -> pd.DataFrame:
    d = ev[~ev['kept'] & (ev['min_jc'] >= 1)]
    return (d.groupby(['gene_id', 'gene'], dropna=False)
            .agg(events=('ID', 'size'), max_reads=('min_jc', 'max'),
                 no_reads_in_group_1=('reason', lambda s: int((s == 'no reads in group 1').sum())),
                 no_reads_in_group_2=('reason', lambda s: int((s == 'no reads in group 2').sum())))
            .reset_index().sort_values(['max_reads', 'events'], ascending=False))


def model_cutoffs(ev) -> list:
    """ the -m cutoff fitted to the kept events (as JCAST fits it), with the dropped events that have a read in both
    forms, and to every event (as JCAST would read the --statoff folder: with the counts of 0) """
    from jcast import model
    lines = []
    for name, e in [('kept events (as JCAST reads the folder)', ev[ev['kept']]),
                    ('kept, and dropped with a read in both forms', ev[ev['kept'] | (ev['min_jc'] >= 1)]),
                    ('every event (as JCAST would read the --statoff folder)', ev)]:
        _, params, cutoff, why = model.fit(e['min_jc'].to_numpy() + 1)
        lines.append(f'-m cutoff, {name}: {cutoff if why is None else why} ({len(e)} events)')
    return lines


def write_dropped(ev, statoff_dir, dropped_dir):
    """ the dropped events with a read in both forms, as an rMATS folder; FDR and PValue 0 so that -q keeps them """
    Path(dropped_dir).mkdir(parents=True, exist_ok=True)
    for st in cr.SPLICE_TYPES:
        ids = set(ev.loc[(ev['splice_type'] == st) & ~ev['kept'] & (ev['min_jc'] >= 1), 'ID'])
        m = pd.read_csv(Path(statoff_dir) / f'{st}.MATS.JC.txt', sep='\t', dtype=str, keep_default_na=False)
        m = m[m['ID'].astype(int).isin(ids)].assign(PValue='0', FDR='0')
        m.columns = [c.split('.')[0] for c in m.columns]                # ID.1 back to ID
        m.to_csv(Path(dropped_dir) / f'{st}.MATS.JC.txt', sep='\t', index=False)


def measure(a, cutoffs, ev) -> pd.DataFrame:
    jx = pd.read_csv(Path(a.dropped_run) / 'junctions.tsv', sep='\t', dtype={'junction_id': str})
    if a.dropped_records:
        kept = cr.read_records(a.records)
        drop = cr.read_records(a.dropped_records)
    else:                       # one -r 1 run on --statoff: its records and junctions, split by event
        dropped = set(ev.loc[~ev['kept'], 'splice_type'] + ':' + ev.loc[~ev['kept'], 'ID'].astype(str))
        rec = cr.read_records(a.records)
        of_dropped = (rec['splice_type'] + ':' + rec['junction_id']).isin(dropped)
        kept, drop = rec[~of_dropped], rec[of_dropped]
        jx = jx[(jx['junction_type'] + ':' + jx['junction_id']).isin(dropped)]
    sp_iso, sp_can, sp_all = cr.swissprot(a.uniprot)
    gc, translations = cr.gencode(a.gencode, a.gtf, a.cache_dir, sp_can)
    known = {'Swiss-Prot isoforms': sp_iso, 'GENCODE alternative proteins': gc['GENCODE']}
    if a.salmon:
        tpm = cr.expression(a.salmon, translations, gc['GENCODE'])
        known['the same, >= 1 TPM'] = set(tpm[tpm >= 1].index)

    def database(rec):          # the merged FASTA's JCAST records (cutoff_recovery.recovery)
        return set(rec.loc[rec['suffix'].isin(cr.MERGE_TIERS) & rec['header'].str.startswith('sp|')
                           & ~rec['sha1'].isin(sp_all), 'sha1'])

    rows = []
    for c in [1, *cutoffs]:
        k, d = kept[kept['min_jc'] >= c], drop[drop['min_jc'] >= c]
        j = jx[jx['min_jc'] >= c]
        row = {'cutoff': c, 'dropped junctions': len(j)}
        for f in cr.TRANSLATED:
            row[f'{f} junctions'] = int((j['fate'] == f).sum())
        ks, ds = set(k['sha1']), set(d['sha1'])
        row['records: kept events'] = len(ks)
        row['records: dropped events, not given by a kept one'] = len(ds - ks)
        kd, dd = database(k), database(d)
        row['database records: kept events'] = len(kd)
        row['database records: added'] = len(dd - kd)
        for name, seqs in known.items():
            row[f'{name}: kept events'] = len(ks & seqs)
            row[f'{name}: added'] = len((ds - ks) & seqs)
        rows.append(row)
    return pd.DataFrame(rows).set_index('cutoff').T.reset_index().rename(columns={'index': 'measure'})


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--rmats', default=None, help='the rMATS folder, statistics on, as JCAST reads it; without it, '
                                                 "filter_countfile's conditions on --statoff")
    p.add_argument('--statoff', required=True, help="a post step with --statoff on the same samples' prep files")
    p.add_argument('--dropped', default=None, help='prepare: write the dropped events here, as an rMATS folder')
    p.add_argument('--records', default=None, help='measure: the audit records of the -r 1 run on --rmats (without '
                                                   '--dropped_records: on --statoff)')
    p.add_argument('--dropped_records', default=None, help='measure: the audit records of the -r 1 run on --dropped')
    p.add_argument('--dropped_run', default=None, help='measure: that run\'s output dir (without --dropped_records: '
                                                       'the run on --statoff)')
    p.add_argument('--uniprot')
    p.add_argument('--gencode')
    p.add_argument('--gtf')
    p.add_argument('--cache_dir', default=None)
    p.add_argument('--salmon', default=None)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ev = events(a.rmats, a.statoff)
    lines = model_cutoffs(ev)
    if a.rmats:
        lines.append(f"filter_countfile's conditions on the --statoff counts keep exactly the events of --rmats: "
                     f"{bool((ev['kept'] == ev['passes']).all())} ({int(ev['passes'].sum())} events pass)")
    cutoff, statoff_cutoff = (int(line.split(': ')[1].split(' ')[0]) for line in lines[:2])
    bt = by_type(ev, cutoff)
    bt.to_csv(out / 'dropped.tsv', sep='\t', index=False)
    bb = by_bin(ev)
    bb.to_csv(out / 'dropped_bins.tsv', sep='\t', index=False)
    bg = by_gene(ev)
    bg.to_csv(out / 'dropped_genes.tsv', sep='\t', index=False)
    (out / 'dropped_model.txt').write_text('\n'.join(lines) + '\n')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print('\n'.join(lines))
        print(bt.T.to_string())
        print(bb.to_string(index=False))
        print(bg.head(40).to_string(index=False))
    if a.dropped:
        write_dropped(ev, a.statoff, a.dropped)
    if a.records:
        kn = measure(a, [cutoff, statoff_cutoff], ev)
        kn.to_csv(out / 'dropped_known.tsv', sep='\t', index=False)
        print(kn.to_string(index=False))


if __name__ == '__main__':
    main()
