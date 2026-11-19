# -*- coding: utf-8 -*-

""" PR 24: rMATS folders made without the statistics (--statoff) or from one sample
(reports/2026-11-15-read-count-model-audit.md §§15.7, 15.8). With its statistics on, rMATS leaves out the events with
no reads in one sample, no inclusion read or no skipping read; without them it writes every event, with PValue and FDR
NA, and from one sample, the second sample's count fields empty. The folders here are the bundled chr15 folder
rewritten that way, as rMATS 4.4.0 writes them """

import argparse
import csv
import math

import numpy as np
import pandas as pd
import pytest

from jcast import fates, main, model, params
from jcast.junctions import RmatsResults, read_count
from jcast.outcomes import SPLICE_TYPES
from jcast.sequences import Sequence

from tests.conftest import DATA, HUMAN_GENOME, HUMAN_GTF, RMATS
from tests.test_uniprot import UNIPROT, _no_api

COUNTS = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2']

# events added to each file, as copies of its first rows with new IDs: what rMATS's statistics leave out, and what
# they leave out but JCAST keeps
NO_READS = {'IJC_SAMPLE_1': 0, 'SJC_SAMPLE_1': 0, 'IJC_SAMPLE_2': 0, 'SJC_SAMPLE_2': 0}
INCLUSION_ONLY = {'SJC_SAMPLE_1': 0, 'SJC_SAMPLE_2': 0}
SKIPPING_ONLY = {'IJC_SAMPLE_1': 0, 'IJC_SAMPLE_2': 0}
LEFT_OUT = [NO_READS, INCLUSION_ONLY, SKIPPING_ONLY]
SAMPLE_1_ONLY = {'IJC_SAMPLE_2': 0, 'SJC_SAMPLE_2': 0}

OWN = ('logfile.log', 'junctions.tsv')        # the output files that are a run's own


def _zeros(field) -> str:
    """ a count field with each replicate's count 0 """
    return ','.join('0' for _ in str(field).split(','))


def _without_statistics(out, added=(), one_sample=False):
    """
    The bundled rMATS folder as rMATS writes it with --statoff: PValue and FDR NA. From one sample (one_sample), also
    IJC_SAMPLE_2, SJC_SAMPLE_2 and IncLevel2 empty, and IncLevelDifference NA. added: events added to each file, each
    a copy of one of its rows with a new ID and these count fields 0

    :return: out, and the IDs added to each file: {splice type: [the ID of each event of added]}
    """
    out.mkdir()
    ids = {}
    for t in SPLICE_TYPES:
        rows = pd.read_table(RMATS / f'{t}.MATS.JC.txt', dtype=str, keep_default_na=False)
        new = []
        for i, counts in enumerate(added):
            row = rows.iloc[i].copy()
            row['ID'] = row['ID.1'] = str(rows['ID'].astype(int).max() + 1 + i)
            for c in counts:
                row[c] = _zeros(row[c])
            new.append(row)
        ids[t] = [r['ID'] for r in new]
        rows = pd.concat([rows, pd.DataFrame(new)], ignore_index=True)
        rows['PValue'] = rows['FDR'] = 'NA'
        if one_sample:
            rows['IJC_SAMPLE_2'] = rows['SJC_SAMPLE_2'] = rows['IncLevel2'] = ''
            rows['IncLevelDifference'] = 'NA'
        rows.columns = [c.split('.')[0] for c in rows.columns]     # ID.1 back to ID
        rows.to_csv(out / f'{t}.MATS.JC.txt', sep='\t', index=False)
    return out, ids


def _args(rmats, out, cache_dir, **extra):
    args = dict(rmats_folder=str(rmats), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME), out=str(out), read=1,
                model=False, canonical=False, qvalue=[0, 1], splice_type={'RI', 'A5SS'}, uniprot=str(UNIPROT),
                processes=1, cache_dir=str(cache_dir), verbose=False, overwrite=False)
    args.update(extra)
    return argparse.Namespace(**args)


def _run(rmats, out, cache_dir, monkeypatch, **extra):
    """ a run with --uniprot (no network), and its junctions.tsv, the IDs as strings """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    monkeypatch.setattr(params, 'use_gtf_only', False)
    main.run_jcast(_args(rmats, out, cache_dir, **extra))
    return pd.read_table(out / 'junctions.tsv', dtype=str, keep_default_na=False)


# ---- the read count: an empty field counts 0 for its sample only ----

@pytest.mark.parametrize('counts, expected', [
    (('4,6', '3,5', float('nan'), float('nan')), 8),    # one sample: min(ijc 10, sjc 8)
    (('4,6', '3,5', '', ''), 8),
    ((float('nan'), float('nan'), 2, 7), 2),            # the other one
    (('0,1', '0,0', float('nan'), float('nan')), 0),    # no skipping read
])
def test_read_count_of_one_sample(counts, expected):
    """ 2026-11-18 an empty field made both samples' counts 0, so a folder from one sample had no reads """
    assert read_count(*counts) == expected


# ---- the events left out as the folder is read ----

def test_a_folder_with_statistics_loses_nothing(logger):
    """ rMATS's statistics have left these events out already: every bundled folder is read as before """
    for folder in [RMATS] + sorted(p for p in (DATA / 'model').iterdir() if p.is_dir()):
        results = RmatsResults(logger=logger, rmats_dir=folder)
        assert set(results.left_out.values()) <= {0}, folder.name
        for t, table in results.tables.items():
            raw = pd.read_table(folder / f'{t}.MATS.JC.txt', sep='\t', low_memory=False)
            assert len(table) == len(raw)


def test_events_without_a_form_are_left_out(logger, tmp_path):
    """ no reads, no skipping read or no inclusion read over all samples: left out of the tables, the counts and the
    model's array. An event with reads in one sample only is kept """
    folder, ids = _without_statistics(tmp_path / 'statoff', added=LEFT_OUT + [SAMPLE_1_ONLY])
    results = RmatsResults(logger=logger, rmats_dir=folder)
    bundled = RmatsResults(logger=logger, rmats_dir=RMATS)

    assert results.left_out == {t: 3 for t in SPLICE_TYPES}
    for t in SPLICE_TYPES:
        kept = results.tables[t].id.astype(str).tolist()
        assert kept == bundled.tables[t].id.astype(str).tolist() + [ids[t][3]]
        assert results.tables[t].min_jc.min() >= 1
    results.get_junction_count_array()
    assert len(results.counts) == len(results.min_jc_array) == sum(len(t) for t in bundled.tables.values()) + 5


def test_a_form_without_length_is_left_out(logger, tmp_path):
    """ as rMATS leaves it out with its statistics; none of the folders measured has one """
    folder, _ = _without_statistics(tmp_path / 'statoff')
    path = folder / 'SE.MATS.JC.txt'
    rows = pd.read_table(path, dtype=str, keep_default_na=False)
    rows.loc[0, 'IncFormLen'] = '0'
    rows.loc[1, 'SkipFormLen'] = '0'
    rows.columns = [c.split('.')[0] for c in rows.columns]
    rows.to_csv(path, sep='\t', index=False)

    assert RmatsResults(logger=logger, rmats_dir=folder, splice_types=['SE']).left_out == {'SE': 2}


def test_the_model_does_not_see_them(tmp_path):
    """ -m on a --statoff folder: the events left out do not move the cutoff (26, the bundled folder's). Fitted with
    their counts of 0, as before, it was 1 """
    folder, _ = _without_statistics(tmp_path / 'statoff', added=LEFT_OUT * 20)
    counts, _, _, cutoff, no_fit = model.preview(folder)

    assert (len(counts), cutoff, no_fit) == (1105, 26, None)
    with_zeros = np.concatenate([np.array([[c + 1] for c in counts.min_jc]), np.ones((3 * 20 * 5, 1), dtype=int)])
    assert model.fit(with_zeros)[2] < 26


# ---- an NA FDR: kept by -q 0 1 only ----

@pytest.mark.parametrize('fdr, cutoff, qvalue, expected', [
    (math.nan, 1, (0, 1), None),
    (math.nan, 1, [0.0, 1.0], None),            # as the command line gives it
    (math.nan, 1, (0, 0.05), fates.SKIPPED_P_VALUE),
    (math.nan, 1, (0.01, 1), fates.SKIPPED_P_VALUE),
    (math.nan, 9, (0, 1), fates.SKIPPED_LOW_COUNT),
    (0.01, 1, (0, 0.05), None),
    (0.5, 1, (0, 0.05), fates.SKIPPED_P_VALUE),
])
def test_filtered_without_an_fdr(fdr, cutoff, qvalue, expected):
    """ 2026-11-18 a junction without an FDR was always skipped_fdr, so a --statoff folder gave an empty database """
    assert fates.filtered(5, fdr, cutoff, qvalue) == expected


def test_preview_of_a_folder_without_statistics(tmp_path):
    """ the GUI's Model page: -q 0 1 keeps what the bundled folder's -q 0 1 keeps; any other range skips it all """
    folder, _ = _without_statistics(tmp_path / 'statoff')
    counts = model.preview(folder)[0]
    bundled = model.preview(RMATS)[0]

    assert model.filter_counts(counts, 5, (0, 1)).equals(model.filter_counts(bundled, 5, (0, 1)))
    narrow = model.filter_counts(counts, 5, (0, 0.05))
    assert (narrow.kept == 0).all() and (narrow.skipped_fdr == narrow.junctions - narrow.skipped_low_count).all()


# ---- runs ----

@pytest.mark.slow
def test_a_run_without_statistics(tmp_path, cache_dir, monkeypatch):
    """ a --statoff folder gives the bundled folder's run, with FDR NA; the events left out are not in junctions.tsv,
    and the log says how many there were """
    folder, ids = _without_statistics(tmp_path / 'statoff', added=LEFT_OUT)
    statoff = _run(folder, tmp_path / 'out', cache_dir, monkeypatch)
    bundled = _run(RMATS, tmp_path / 'bundled', cache_dir, monkeypatch)

    assert (statoff.fdr == 'NA').all()
    assert statoff.drop(columns='fdr').equals(bundled.drop(columns='fdr'))
    for t in ['RI', 'A5SS']:
        assert not set(ids[t]) & set(statoff.junction_id)
    files = sorted(p.name for p in (tmp_path / 'bundled').iterdir() if p.name not in OWN)
    assert files == sorted(p.name for p in (tmp_path / 'out').iterdir() if p.name not in OWN)
    assert 'psq_T1.fasta' in files
    for name in files:
        assert (tmp_path / 'out' / name).read_bytes() == (tmp_path / 'bundled' / name).read_bytes(), name

    log = (tmp_path / 'out' / 'logfile.log').read_text()
    for t in ['RI', 'A5SS']:
        assert (f'rMATS {t}: 3 events left out, with no inclusion or no skipping read in any sample, as rMATS leaves '
                f'them out when it runs its statistics.') in log
    assert (f'{len(statoff)} junctions have no rMATS FDR (rMATS was run with --statoff, or one sample); -q 0 1 keeps '
            f'them.') in log
    assert 'events left out' not in (tmp_path / 'bundled' / 'logfile.log').read_text()


@pytest.mark.slow
def test_a_narrower_q_skips_them_and_warns(tmp_path, cache_dir, monkeypatch, capsys):
    folder, _ = _without_statistics(tmp_path / 'statoff')
    table = _run(folder, tmp_path / 'out', cache_dir, monkeypatch, qvalue=[0, 0.05])

    assert set(table.fate) == {'skipped_fdr'}
    assert (f'{len(table)} junctions have no rMATS FDR (rMATS was run with --statoff, or one sample); -q 0 0.05 skips '
            f'them (skipped_fdr). -q 0 1 keeps them.') in capsys.readouterr().err


@pytest.mark.slow
def test_a_run_from_one_sample(tmp_path, cache_dir, monkeypatch):
    """ the second sample's fields empty: the read counts are the first sample's, and each junction kept has the fate
    it has in the bundled folder's run. 2026-11-18 every read count was 0 """
    folder, _ = _without_statistics(tmp_path / 'one', one_sample=True)
    one = _run(folder, tmp_path / 'out', cache_dir, monkeypatch)
    bundled = _run(RMATS, tmp_path / 'bundled', cache_dir, monkeypatch).set_index(['junction_type', 'junction_id'])

    rows = pd.concat([pd.read_table(RMATS / f'{t}.MATS.JC.txt', dtype=str) for t in ['RI', 'A5SS']], ignore_index=True)
    sample_1 = {i: read_count(ijc, sjc, '', '') for i, ijc, sjc in zip(rows.ID, rows.IJC_SAMPLE_1, rows.SJC_SAMPLE_1)}
    assert sorted(one.junction_id, key=int) == sorted((i for i, n in sample_1.items() if n), key=int)
    assert [int(n) for n in one.min_jc] == [sample_1[i] for i in one.junction_id]
    assert (one.fdr == 'NA').all() and len(one) > 50
    expected = bundled.loc[list(zip(one.junction_type, one.junction_id))]
    assert one.fate.tolist() == expected.fate.tolist()


@pytest.mark.slow
def test_a_run_from_one_sample_with_m(tmp_path, cache_dir, monkeypatch):
    """ -m fits the first sample's counts, and the log names it """
    folder, _ = _without_statistics(tmp_path / 'one', one_sample=True)
    _run(folder, tmp_path / 'out', cache_dir, monkeypatch, model=True)
    counts = model.preview(folder)[0]

    assert counts.min_jc.min() >= 1
    cutoff = model.fit(np.array([[c + 1] for c in counts.min_jc]))[2]
    log = (tmp_path / 'out' / 'logfile.log').read_text()
    assert f'Translating junctions with a read count of at least {cutoff} (-m)' in log
