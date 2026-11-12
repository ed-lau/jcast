# -*- coding: utf-8 -*-

""" Translating in several processes, -p (PR 15; reports/2026-11-12-parallel-audit.md) """

import argparse
import logging
import multiprocessing

import pytest

from jcast import main, params
from jcast.sequences import Sequence

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS
from tests.test_uniprot import UNIPROT, _no_api


def _run(tmp_path, name, processes):
    """ run_jcast on the chr15 test set with --uniprot; returns the output directory """
    logging.getLogger('jcast').handlers = []     # get_logger adds handlers on every call (REVIEW.md MED-13)
    gtf_copy = tmp_path / HUMAN_GTF.name
    if not gtf_copy.exists():
        gtf_copy.write_bytes(HUMAN_GTF.read_bytes())
    out = tmp_path / name
    out.mkdir()
    args = argparse.Namespace(rmats_folder=str(RMATS), gtf_file=str(gtf_copy), genome=str(HUMAN_GENOME),
                              out=str(out), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1],
                              splice_type={'MXE', 'SE', 'RI', 'A5SS', 'A3SS'}, uniprot=str(UNIPROT),
                              processes=processes)
    main.run_jcast(args)
    logging.getLogger('jcast').handlers = []
    return out


def _messages(out):
    """ log messages from the first junction on, without timestamps """
    lines = (out / 'logfile.log').read_text().split('\n')
    messages = [line.split('~', 1)[1] if line[:4].isdigit() else line for line in lines]
    first = next(i for i, m in enumerate(messages) if '>>>>>> Doing' in m)
    return messages[first:]


@pytest.mark.slow
@pytest.mark.skipif('fork' not in multiprocessing.get_all_start_methods(), reason='needs fork')
def test_processes_give_the_same_output_and_log(tmp_path, monkeypatch):
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    monkeypatch.setattr(params, 'use_gtf_only', False)

    one = _run(tmp_path, 'p1', 1)
    three = _run(tmp_path, 'p3', 3)

    files = sorted(p.name for p in one.iterdir() if p.name != 'logfile.log')
    assert files == sorted(p.name for p in three.iterdir() if p.name != 'logfile.log')
    for name in files:
        assert (one / name).read_bytes() == (three / name).read_bytes(), name
    assert _messages(one) == _messages(three)
    assert 'Translating junctions in 3 worker processes.' in (three / 'logfile.log').read_text()


def test_more_than_one_process_needs_uniprot():
    logger = logging.getLogger('jcast-test')
    assert main._number_of_processes(argparse.Namespace(processes=4), None, logger) == 1
    assert main._number_of_processes(argparse.Namespace(processes=4), object(), logger) == (
        4 if 'fork' in multiprocessing.get_all_start_methods() else 1)
    assert main._number_of_processes(argparse.Namespace(), None, logger) == 1


@pytest.mark.parametrize('extra, message', [
    (['-p', '2'], '-p/--processes > 1 needs --uniprot'),
    (['-p', '0'], '-p/--processes must be at least 1'),
])
def test_cli_rejects_processes(monkeypatch, capsys, extra, message):
    monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), str(HUMAN_GTF), str(HUMAN_GENOME)] + extra)

    with pytest.raises(SystemExit):
        main.main()
    assert message in capsys.readouterr().err
