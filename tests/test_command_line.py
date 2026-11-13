# -*- coding: utf-8 -*-

""" The command line's checks, and the rMATS files read (PR 12; reports/2026-11-13-command-line-checks.md) """

import csv
import shutil

import pytest

from jcast import main
from jcast.junctions import RmatsResults

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS


def _folder_with(tmp_path, *splice_types):
    """ an rMATS folder holding only these splice types' files """
    folder = tmp_path / 'rmats'
    folder.mkdir()
    for splice_type in splice_types:
        shutil.copy(RMATS / f'{splice_type}.MATS.JC.txt', folder)
    return folder


def _cli(monkeypatch, *argv):
    monkeypatch.setattr('sys.argv', ['jcast'] + [str(a) for a in argv])
    main.main()


def test_only_the_selected_files_are_read(tmp_path, logger):
    """ 2026-11-13 all five were read, so a folder without MXE.MATS.JC.txt failed even with -s SE (MED-8) """
    results = RmatsResults(logger, _folder_with(tmp_path, 'SE'), splice_types={'SE'})

    assert list(results.tables) == ['SE']
    assert len(results.rmats_se) == 720 and results.rmats_mxe is None


def test_the_model_reads_the_other_files_present(logger):
    """ with -m, the files -s does not select are read too, when present: the model fits on every junction """
    for others, junctions in [(False, 720), (True, 1105)]:
        results = RmatsResults(logger, RMATS, splice_types={'SE'}, others_if_present=others)
        results.get_junction_count_array()
        assert len(results.min_jc_array) == junctions


def test_a_run_with_only_the_selected_file(tmp_path, monkeypatch, cache_dir):
    """ 2026-11-13 this failed with FileNotFoundError for MXE.MATS.JC.txt (MED-8) """
    out = tmp_path / 'out'
    _cli(monkeypatch, _folder_with(tmp_path, 'SE'), HUMAN_GTF, HUMAN_GENOME, '-s', 'SE', '-r', 100000, '-o', out,
         '--cache_dir', cache_dir)

    with open(out / 'junctions.tsv') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    assert len(rows) == 720 and {r['junction_type'] for r in rows} == {'SE'}


@pytest.mark.slow
def test_the_model_cutoff_with_fewer_splice_types(tmp_path, monkeypatch, cache_dir):
    """ -m -s SE fits on the five files, as before PR 12, and the log names the files read only for the model; on SE
    alone, the cutoff would be 6 (reports/2026-11-13-command-line-checks.md §2) """
    out = tmp_path / 'out'
    _cli(monkeypatch, RMATS, HUMAN_GTF, HUMAN_GENOME, '-s', 'SE', '-m', '-o', out, '--cache_dir', cache_dir)

    log = (out / 'logfile.log').read_text()
    assert 'modeled read count cutoff 26.' in log
    assert 'rMATS SE: ' in log and 'rMATS MXE (not selected by -s; read for the -m model): ' in log


@pytest.mark.parametrize('splice_types, missing', [(['SE', 'MXE'], 'MXE.MATS.JC.txt'), ([], 'MXE.MATS.JC.txt')])
def test_a_missing_selected_file_is_a_usage_error(tmp_path, monkeypatch, capsys, splice_types, missing):
    """ 2026-11-13 an assert, which python -O skips """
    selected = ['-s'] + splice_types if splice_types else []
    with pytest.raises(SystemExit) as e:
        _cli(monkeypatch, _folder_with(tmp_path, 'SE'), HUMAN_GTF, HUMAN_GENOME, *selected, '-o', tmp_path / 'out')

    assert e.value.code == 2
    err = capsys.readouterr().err
    assert 'rMATS file not found: ' in err and missing in err and 'SE.MATS.JC.txt' not in err


def test_a_missing_rmats_folder_is_a_usage_error(tmp_path, monkeypatch, capsys):
    with pytest.raises(SystemExit) as e:
        _cli(monkeypatch, tmp_path / 'none', HUMAN_GTF, HUMAN_GENOME, '-o', tmp_path / 'out')

    assert e.value.code == 2
    assert 'rMATS folder not found: ' in capsys.readouterr().err


def test_g_or_ln_is_gone(tmp_path, monkeypatch, capsys):
    """ 2026-11-13 --g_or_ln was accepted and read by nothing """
    with pytest.raises(SystemExit) as e:
        _cli(monkeypatch, tmp_path / 'none', HUMAN_GTF, HUMAN_GENOME, '--g_or_ln', 1)

    assert e.value.code == 2
    assert 'unrecognized arguments: --g_or_ln' in capsys.readouterr().err
