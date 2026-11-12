# -*- coding: utf-8 -*-

""" Real runs (PR 7; reports/2026-11-12-thresholds-runs-canonical-audit.md §3): the output directory, reruns, one
failing junction, the logger across runs, and the caches """

import argparse
import logging
import multiprocessing
import os

import gtfparse
import pandas as pd
import pytest
import requests
from Bio import SeqIO

from jcast import annots, cache, main, params
from jcast.annots import ReadAnnotations
from jcast.sequences import Sequence

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS
from tests.test_junction import _junction
from tests.test_uniprot import UNIPROT, _no_api

FAULT = ('SE', '11950')     # the junction made to raise

API_LOOKUP = Sequence.get_canonical_aa_uniprot


def _args(out, cache_dir, **extra):
    """ the chr15 test set with --uniprot, as the CLI passes it """
    args = dict(rmats_folder=str(RMATS), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME), out=str(out), read=1,
                model=False, canonical=True, mask=None, qvalue=[0, 1], splice_type={'RI'}, uniprot=str(UNIPROT),
                processes=1, cache_dir=str(cache_dir))
    args.update(extra)
    return argparse.Namespace(**args)


def _sequences(out):
    """ psq file -> set of sequences """
    return {p.name: {str(r.seq) for r in SeqIO.parse(p, 'fasta')} for p in sorted(out.glob('psq_*.fasta'))}


# log messages that differ between runs of the same input
VARYING = ('INFO~Time: ', 'INFO~Translating junctions in ')


def _messages(out):
    """ log messages from the read-count cutoff on, without time stamps, the run time or the number of processes """
    lines = (out / 'logfile.log').read_text().split('\n')
    messages = [line.split('~', 1)[1] if line[:4].isdigit() else line for line in lines]
    start = next(i for i, m in enumerate(messages) if 'Translating junctions with' in m)
    return [m for m in messages[start:] if not m.startswith(VARYING)]


def _table(out):
    """ junctions.tsv """
    return pd.read_table(out / 'junctions.tsv', dtype=str, keep_default_na=False)


@pytest.fixture(autouse=True)
def _uniprot_only(monkeypatch):
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    monkeypatch.setattr(params, 'use_gtf_only', False)


# ---- the output directory (MED-4, MED-12) ----

def test_output_directory_is_created(tmp_path, cache_dir):
    out = tmp_path / 'new' / 'dir'
    main.run_jcast(_args(out, cache_dir))
    assert (out / 'logfile.log').exists() and list(out.glob('psq_*.fasta'))


def test_directory_with_output_is_refused(tmp_path, cache_dir):
    out = tmp_path / 'out'
    out.mkdir()
    (out / 'psq_T1.fasta').write_text('>old\nMOLD\n')

    with pytest.raises(FileExistsError, match='--overwrite'):
        main.run_jcast(_args(out, cache_dir))
    assert (out / 'psq_T1.fasta').read_text() == '>old\nMOLD\n'


def test_overwrite_replaces_only_jcast_output(tmp_path, cache_dir):
    """ a rerun with --overwrite gives the output of a fresh run; other files stay """
    fresh, reused = tmp_path / 'fresh', tmp_path / 'reused'
    main.run_jcast(_args(fresh, cache_dir))
    main.run_jcast(_args(reused, cache_dir, read=30))
    (reused / 'notes.txt').write_text('mine')

    main.run_jcast(_args(reused, cache_dir, overwrite=True))

    assert _sequences(reused) == _sequences(fresh)
    assert (reused / 'notes.txt').read_text() == 'mine'
    assert _messages(reused) == _messages(fresh)
    assert (reused / 'logfile.log').read_text().count('~INFO~Arguments: Namespace(') == 1
    assert (reused / 'junctions.tsv').read_bytes() == (fresh / 'junctions.tsv').read_bytes()


def test_cli_refuses_directory_with_output(tmp_path, monkeypatch, capsys):
    (tmp_path / 'logfile.log').write_text('')
    monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), str(HUMAN_GTF), str(HUMAN_GENOME), '-o', str(tmp_path)])

    with pytest.raises(SystemExit):
        main.main()
    assert 'already holds JCAST output (logfile.log); use --overwrite' in capsys.readouterr().err


def test_other_files_do_not_block_a_run(tmp_path):
    (tmp_path / 'notes.txt').write_text('mine')
    assert main.existing_outputs(tmp_path) == []
    assert main.existing_outputs(tmp_path / 'missing') == []


# ---- one failing junction (masked-bases audit §7, item 3) ----

@pytest.mark.slow
@pytest.mark.skipif('fork' not in multiprocessing.get_all_start_methods(), reason='needs fork')
def test_failing_junction_does_not_stop_the_run(tmp_path, monkeypatch, capsys, cache_dir, gtf, genome,
                                               rmats_results):
    translate_one = main._translate_one

    def faulty(junction, *args, **kwargs):
        if (junction.junction_type, junction.name) == FAULT:
            raise RuntimeError('injected fault')
        return translate_one(junction, *args, **kwargs)

    # -r 2: at -r 1 another junction writes the same sequences as the failing one
    normal = tmp_path / 'normal'
    main.run_jcast(_args(normal, cache_dir, splice_type={'SE'}, read=2))
    capsys.readouterr()

    monkeypatch.setattr(main, '_translate_one', faulty)
    outs = {}
    for processes in (1, 3):
        outs[processes] = tmp_path / f'fault_p{processes}'
        main.run_jcast(_args(outs[processes], cache_dir, splice_type={'SE'}, processes=processes, read=2))
        err = capsys.readouterr().err
        assert err.count('Junction SE 11950 of gene') == 1 and 'injected fault' in err
        assert 'Traceback' not in err and '1 junction(s) failed with an exception' in err

    # the same output and log in one process and in three
    assert _sequences(outs[1]) == _sequences(outs[3])
    for name in list(_sequences(outs[1])) + ['junctions.tsv']:
        assert (outs[1] / name).read_bytes() == (outs[3] / name).read_bytes(), name
    assert [m.replace(str(outs[1]), '<out>') for m in _messages(outs[1])] == \
           [m.replace(str(outs[3]), '<out>') for m in _messages(outs[3])]

    # the failure is in junctions.tsv, and its traceback in the log
    table = _table(outs[1])
    assert len(table) == 720
    failed = table[table['fate'] == 'exception']
    assert failed[['junction_type', 'junction_id', 'fate_detail', 'files']].values.tolist() == \
           [['SE', '11950', 'RuntimeError: injected fault', '']]
    log = (outs[1] / 'logfile.log').read_text()
    assert 'Traceback (most recent call last)' in log and '>>>>>> Doing' not in log

    # every other junction's output is unchanged: only records of the failing junction can be missing
    monkeypatch.setattr(main, '_translate_one', translate_one)
    args = _args(tmp_path / 'alone', cache_dir, read=2)
    from jcast.uniprot import UniProtFasta
    fate, records, _ = main._translate_safely(_junction(rmats_results.rmats_se, int(FAULT[1])), gtf, genome, args,
                                              pred_bound=2, uniprot=UniProtFasta(logging.getLogger('jcast-test'),
                                                                                 UNIPROT))
    own = {str(r.seq) for _, _, r in records}
    with_fault, without = _sequences(outs[1]), _sequences(normal)
    assert set(with_fault) <= set(without)
    for name, seqs in without.items():
        assert with_fault.get(name, set()) <= seqs
        assert seqs - with_fault.get(name, set()) <= own, name
    assert sum(len(s) for s in without.values()) > sum(len(s) for s in with_fault.values())


def test_exception_fate(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ the failing junction's fate names the exception; its records are not kept """
    monkeypatch.setattr(main, '_translate_one', lambda *a, **k: (_ for _ in ()).throw(KeyError('chr99')))
    fate, records, found = main._translate_safely(_junction(rmats_results.rmats_se, 6426), gtf, genome,
                                                  _args(tmp_path, tmp_path), pred_bound=1, uniprot=None)
    assert fate.startswith("FAILURE (exception). KeyError: 'chr99'.") and records == [] and found == {}


# ---- the logger across runs (MED-13) ----

def test_two_runs_in_one_process_log_separately(tmp_path, cache_dir):
    first, second = tmp_path / 'first', tmp_path / 'second'
    main.run_jcast(_args(first, cache_dir))
    first_log = (first / 'logfile.log').read_text()
    main.run_jcast(_args(second, cache_dir))

    assert (first / 'logfile.log').read_text() == first_log
    assert _messages(first) == _messages(second)
    assert logging.getLogger('jcast').handlers == []


# ---- the caches (MED-5, MED-6) ----

GTF_TEXT = (
    '15\tensembl\tgene\t100\t400\t.\t+\t.\tgene_id "G1"; gene_name "A"; gene_biotype "protein_coding";\n'
    '15\tensembl\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "A"; '
    'transcript_biotype "protein_coding"; transcript_support_level "NA";\n'
    '15\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1"; '
    'transcript_support_level "NA";\n'
    '15\tensembl\tCDS\t110\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "1"; '
    'transcript_support_level "NA";\n'
    '15\tensembl\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "10"; '
    'transcript_support_level "1";\n'
)


def _gtf_files(directory):
    return sorted(p.name for p in directory.iterdir())


def test_gtf_cache_in_the_cache_directory(tmp_path, monkeypatch):
    """ the cache is not written next to the GTF; the table is the one the old <gtf>.cached gave, without its index """
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'cache'))
    gtf_dir = tmp_path / 'ref'
    gtf_dir.mkdir()
    path = gtf_dir / 'a.gtf'
    path.write_text(GTF_TEXT)

    annot = ReadAnnotations(logging.getLogger('jcast-test'), path).annot

    assert _gtf_files(gtf_dir) == ['a.gtf']
    assert len(list((tmp_path / 'cache' / 'gtf').glob('a.gtf.*.pkl'))) == 1
    gtfparse.read_gtf(path).to_csv(tmp_path / 'old.cached', encoding='utf-8', sep='\t')
    old = pd.read_table(tmp_path / 'old.cached', sep='\t', low_memory=False).drop(columns='Unnamed: 0')
    assert annot.equals(old)
    assert sorted(annot['exon_number'].dropna()) == [1, 1, 10] and annot['transcript_support_level'].isna().sum() == 4


def test_gtf_cache_is_used_then_rebuilt_when_the_gtf_changes(tmp_path, monkeypatch):
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'cache'))
    path = tmp_path / 'a.gtf'
    path.write_text(GTF_TEXT)
    logger = logging.getLogger('jcast-test')
    first = ReadAnnotations(logger, path).annot

    read_gtf = gtfparse.read_gtf
    monkeypatch.setattr(gtfparse, 'read_gtf', lambda p: pytest.fail('the cache was not used'))
    assert ReadAnnotations(logger, path).annot.equals(first)

    monkeypatch.setattr(gtfparse, 'read_gtf', read_gtf)
    path.write_text(GTF_TEXT.replace('gene_name "A"', 'gene_name "B"'))
    assert set(ReadAnnotations(logger, path).annot['gene_name'].dropna()) == {'B'}
    assert len(list((tmp_path / 'cache' / 'gtf').glob('a.gtf.*.pkl'))) == 1


def test_damaged_gtf_cache_is_read_again(tmp_path, monkeypatch, caplog):
    """ a damaged cache is replaced; a leftover partial file of another run is never read """
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'cache'))
    path = tmp_path / 'a.gtf'
    path.write_text(GTF_TEXT)
    logger = logging.getLogger('jcast-test')
    expected = ReadAnnotations(logger, path).annot
    cached = annots._cache_path(str(path), str(tmp_path / 'cache' / 'gtf'))

    with open(cached, 'wb') as fh:
        fh.write(b'not a pickle')
    with open(cached + '.999.tmp', 'wb') as fh:
        fh.write(b'partial')
    assert ReadAnnotations(logger, path).annot.equals(expected)
    assert 'Cannot read the GTF cache' in caplog.text
    assert ReadAnnotations(logger, path).annot.equals(expected)


def test_unusable_cache_directory(tmp_path, monkeypatch, caplog):
    """ a cache directory that cannot be created: the GTF is read without it """
    blocked = tmp_path / 'blocked'
    blocked.write_text('a file, not a directory')
    monkeypatch.setattr(cache, '_directory', str(blocked / 'cache'))
    monkeypatch.setattr(cache, '_unusable', set())
    path = tmp_path / 'a.gtf'
    path.write_text(GTF_TEXT)

    annot = ReadAnnotations(logging.getLogger('jcast-test'), path).annot

    assert set(annot['gene_id']) == {'G1'}
    assert 'Cannot create the cache directory' in caplog.text


def test_api_cache_in_the_cache_directory(tmp_path, monkeypatch, rmats_results):
    """ answers are cached in the cache directory, not under the working directory (tmp_path) """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', API_LOOKUP)
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'jcast_cache'))
    calls = []

    def get(self, url, **kwargs):
        calls.append(url)
        r = requests.models.Response()
        r.status_code, r._content, r.encoding = 200, b'>sp|P1|A_HUMAN A\nMAAA\n', 'utf-8'
        return r

    monkeypatch.setattr(requests.Session, 'get', get)
    junction = _junction(rmats_results.rmats_se, 6426)
    assert str(Sequence(junction).get_canonical_aa_uniprot().seq) == 'MAAA'
    assert str(Sequence(junction).get_canonical_aa_uniprot().seq) == 'MAAA'

    assert len(calls) == 1
    assert (tmp_path / 'jcast_cache' / 'uniprot-cache.db').exists()
    assert not os.path.exists('cache')
