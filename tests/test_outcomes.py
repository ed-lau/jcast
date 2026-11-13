# -*- coding: utf-8 -*-

""" The run record (PR 16; reports/2026-11-12-log-audit.md): junctions.tsv, the end-of-run summary, and the log with
and without --verbose """

import argparse
import csv

import pandas as pd
import pytest
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import fates, main, params
from jcast import helpers as h
from jcast.outcomes import COLUMNS
from jcast.sequences import Sequence

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS
from tests.test_junction import _junction
from tests.test_uniprot import UNIPROT, _no_api

API_LOOKUP = Sequence.get_canonical_aa_uniprot


def _args(out, cache_dir, **extra):
    args = dict(rmats_folder=str(RMATS), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME), out=str(out), read=1,
                model=False, canonical=True, mask=None, qvalue=[0, 1], splice_type={'RI', 'A5SS'},
                uniprot=str(UNIPROT), processes=1, cache_dir=str(cache_dir))
    args.update(extra)
    return argparse.Namespace(**args)


def _table(out):
    return pd.read_table(out / 'junctions.tsv', dtype=str, keep_default_na=False)


@pytest.fixture(scope='module')
def runs(tmp_path_factory, cache_dir):
    """ the RI and A5SS junctions of chr15 with --uniprot, with and without --verbose """
    patch = pytest.MonkeyPatch()
    patch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    patch.setattr(params, 'use_gtf_only', False)
    tmp = tmp_path_factory.mktemp('outcomes')
    try:
        main.run_jcast(_args(tmp / 'default', cache_dir))
        main.run_jcast(_args(tmp / 'verbose', cache_dir, verbose=True))
    finally:
        patch.undo()
    return tmp / 'default', tmp / 'verbose'


# ---- fate names ----

@pytest.mark.parametrize('fate, expected', [
    (fates.TIER1.format('0', '0'), ('T1', '')),
    (fates.TIER1_ALT_TERMINUS.format('1', '1', 'T1_alt_cterm'), ('T1_alt_terminus', '')),
    (fates.TIER2.format('0/2', '0/2'), ('T2', '')),
    (fates.TIER4.format(2), ('T4', 'slice 2')),
    (fates.FAIL, ('fail', '')),
    (fates.EXCEPTION.format('KeyError', "'chr99'"), ('exception', "KeyError: 'chr99'")),
    (fates.SKIPPED_LOW_COUNT, ('skipped_low_count', '')),
    (fates.SKIPPED_P_VALUE, ('skipped_fdr', '')),
    (fates.SKIPPED_MASKED, ('skipped_masked', '')),
    (fates.SKIPPED_OUTSIDE_CDS, ('skipped_outside_cds', '')),
])
def test_fate_names(fate, expected):
    assert fates.name(fate) == expected


def test_unknown_fate():
    with pytest.raises(ValueError):
        fates.name('SOMETHING ELSE.')


# ---- junctions.tsv ----

def test_one_row_per_junction_in_file_order(runs, rmats_results):
    table = _table(runs[0])
    assert list(table.columns) == COLUMNS
    expected = [('RI', str(i)) for i in rmats_results.rmats_ri['id']] + \
               [('A5SS', str(i)) for i in rmats_results.rmats_a5ss['id']]
    assert list(zip(table['junction_type'], table['junction_id'])) == expected


def test_rows_agree_with_the_output(runs):
    """ every record written is counted once; the files named exist; the canonicals are uniprot_canonicals.tsv's """
    out = runs[0]
    table = _table(out)
    records = sum(1 for p in out.glob('psq_*.fasta') for _ in SeqIO.parse(p, 'fasta'))
    assert table['new_records'].astype(int).sum() == records
    named = {f for files in table['files'] for f in files.split(';') if f}
    assert {f'psq_{f}.fasta' for f in named} == {p.name for p in out.glob('psq_*.fasta')}

    chosen = {r['gene_id']: (r['source'], r['accession'])
              for r in csv.DictReader(open(out / 'uniprot_canonicals.tsv'), delimiter='\t')}
    assert all(chosen[g] == (s, a) for g, s, a in
               zip(table['gene_id'], table['canonical_source'], table['canonical_accession']))


def test_rows_by_fate(runs):
    table = _table(runs[0])
    skipped = table[table['fate'].str.startswith('skipped_low')]
    translated = table[~table['fate'].str.startswith('skipped')]
    assert (skipped['min_jc'].astype(int) < 1).all() if len(skipped) else True
    assert (skipped[['slice1_phase', 'slice2_phase', 'cds_start', 'cds_end']] == '').all().all()
    assert len(translated) > 0 and (translated['cds_start'] != '').any()
    # a T1 junction's slices are each translated in their annotated frame
    t1 = table[table['fate'] == 'T1']
    assert len(t1) > 0 and t1['slice1_phase'].isin(['0', '1', '2']).all() and t1['slice2_phase'].isin(['0', '1', '2']).all()
    assert set(table.loc[table['fate'] == 'T4', 'fate_detail']) <= {'slice 1', 'slice 2'}
    assert (table.loc[table['fate'] != 'T4', 'fate_detail'] == '').all()
    # -c: every junction's canonical is written, so every row names its gene's canonical source
    assert (table['canonical_source'] != '').all()


def test_table_and_output_do_not_depend_on_verbose(runs):
    default, verbose = runs
    for p in sorted(default.iterdir()):
        if p.name != 'logfile.log':
            assert p.read_bytes() == (verbose / p.name).read_bytes(), p.name


# ---- the log ----

def test_default_log_has_no_junction_detail(runs):
    log = (runs[0] / 'logfile.log').read_text()
    assert '>>>>>> Doing' not in log and 'Retrieved nucleotide' not in log and 'Canonical for' not in log
    assert '~module:' not in log
    assert 'Arguments: Namespace(' in log and 'rMATS RI: ' in log and 'sha256 ' in log and 'GTF: ' in log


def test_verbose_log_has_each_junction(runs):
    log = (runs[1] / 'logfile.log').read_text()
    table = _table(runs[1])
    assert log.count('>>>>>> Doing') == len(table)
    assert log.count('~DEBUG~Retrieved nucleotide') == 2 * (~table['fate'].str.startswith('skipped_low')).sum()
    assert '~module:' not in log


def test_summary_counts_the_table(runs):
    out = runs[0]
    table = _table(out)
    log = (out / 'logfile.log').read_text()
    summary = log[log.index('~INFO~Summary.'):].split('\n')
    header = summary[2].split()
    assert header == ['fate', 'RI', 'A5SS', 'total']
    rows = {line.split()[0]: [int(v.replace(',', '')) for v in line.split()[1:]] for line in summary[3:14]}
    assert list(rows) == [s for _, s in fates.NAMES] + ['total']
    for short, (ri, a5ss, total) in rows.items():
        of = table if short == 'total' else table[table['fate'] == short]
        assert [(of['junction_type'] == 'RI').sum(), (of['junction_type'] == 'A5SS').sum(), len(of)] == \
               [ri, a5ss, total]
    records = summary[14]
    for p in sorted(out.glob('psq_*.fasta')):
        assert f'{p.name} {sum(1 for _ in SeqIO.parse(p, "fasta")):,}' in records
    assert summary[15].split('~')[2].startswith('Time: ')


def test_cli_verbose(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(main, 'run_jcast', lambda args: seen.append(args.verbose))
    for flag in ([], ['--verbose'], ['-v']):
        monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), str(HUMAN_GTF), str(HUMAN_GENOME),
                                         '-o', str(tmp_path / 'out')] + flag)
        main.main()
    assert seen == [False, True, True]


# ---- the parts ----

def test_writer_reports_duplicates(tmp_path):
    record = SeqRecord(Seq('MAAA'), id='a')
    assert h.write_seqrecord_to_fasta(record, str(tmp_path), 'T1') is True
    assert h.write_seqrecord_to_fasta(SeqRecord(Seq('MAAA'), id='b'), str(tmp_path), 'T1') is False


def test_canonical_source_gtf_only(gtf, genome, rmats_results, monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    sequence = Sequence(_junction(rmats_results.rmats_se, 6426))
    sequence.get_canonical_aa(gtf, genome.genome)
    assert sequence.canonical_source == ('gtf', '')


@pytest.mark.parametrize('answer, source', [(b'>sp|P1|A_HUMAN A\nMAAA\n', ('sp by API', 'P1')),
                                            (b'>tr|Q9|B_HUMAN B\nMBBB\n', ('tr by API', 'Q9')),
                                            (b'', ('none', ''))])
def test_canonical_source_api(gtf, genome, rmats_results, monkeypatch, tmp_path, answer, source):
    from jcast import cache
    monkeypatch.setattr(params, 'use_gtf_only', False)
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', API_LOOKUP)
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'api_cache'))

    def get(self, url, **kwargs):
        r = requests.models.Response()
        r.status_code, r._content, r.encoding = 200, answer, 'utf-8'
        return r
    monkeypatch.setattr(requests.Session, 'get', get)

    sequence = Sequence(_junction(rmats_results.rmats_se, 6426))
    sequence.get_canonical_aa(gtf, genome.genome)
    assert sequence.canonical_source == source
