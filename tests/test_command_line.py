# -*- coding: utf-8 -*-

""" The command line: its checks and the rMATS files read (PR 12; reports/2026-11-13-command-line-checks.md), and the
command line on Typer (PR 18; reports/2026-11-13-gui-cli-scope.md) """

import argparse
import csv
import json
import shutil
import signal
import subprocess
import sys
import time

import pytest

import jcast.gui
from jcast import main, __version__
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


@pytest.mark.parametrize('option', ['--g_or_ln', '--mask'])
def test_removed_options(tmp_path, monkeypatch, capsys, option):
    """ 2026-11-13 --g_or_ln was accepted and read by nothing (removed in PR 12); --mask was deprecated and ignored
    (removed in PR 18) """
    with pytest.raises(SystemExit) as e:
        _cli(monkeypatch, tmp_path / 'none', HUMAN_GTF, HUMAN_GENOME, option, 1)

    assert e.value.code == 2
    assert f'No such option: {option}' in capsys.readouterr().err


# ---- PR 18: the command line on Typer ----

INPUTS = [RMATS, HUMAN_GTF, HUMAN_GENOME]


def _parsed(monkeypatch, *argv):
    """ the arguments the command line passes to run_jcast """
    calls = []
    monkeypatch.setattr(main, 'run_jcast', lambda args, progress=None: calls.append(args))
    main.main([str(a) for a in argv])
    return calls[0]


def test_run_is_the_default_command(monkeypatch, tmp_path):
    """ jcast RMATS GTF GENOME, as before PR 18, is jcast run RMATS GTF GENOME """
    bare = _parsed(monkeypatch, *INPUTS, '-o', tmp_path, '-m', '-s', 'SE')
    named = _parsed(monkeypatch, 'run', *INPUTS, '-o', tmp_path, '-m', '-s', 'SE')

    assert bare == named
    assert bare.model and bare.splice_type == ['SE']


@pytest.mark.parametrize('argv, splice_types', [
    ([], {'MXE', 'SE', 'RI', 'A5SS', 'A3SS'}),
    (['-s', 'SE', 'MXE'], ['SE', 'MXE']),
    (['-s', 'SE', '-s', 'MXE'], ['SE', 'MXE']),
    (['--splice_type', 'RI', 'A3SS', '-m'], ['RI', 'A3SS']),
])
def test_several_splice_types(monkeypatch, tmp_path, argv, splice_types):
    """ -s SE MXE as before PR 18: Click gives an option a fixed number of values, so main gives each its own -s """
    assert _parsed(monkeypatch, *INPUTS, '-o', tmp_path, *argv).splice_type == splice_types


def test_the_arguments_run_jcast_gets(monkeypatch, tmp_path):
    """ the Namespace argparse gave run_jcast, in its order (the log's Arguments line), less --mask and func """
    args = _parsed(monkeypatch, *INPUTS, '-o', tmp_path, '--overwrite', '--cache_dir', tmp_path, '-r', 3,
                   '-q', 0, 0.05, '-mc', '-v')

    assert list(vars(args)) == ['rmats_folder', 'gtf_file', 'genome', 'out', 'overwrite', 'cache_dir', 'read', 'model',
                                'canonical', 'uniprot', 'processes', 'verbose', 'splice_type', 'qvalue']
    assert (args.out, args.overwrite, args.cache_dir, args.read, args.model, args.canonical, args.verbose) == \
           (str(tmp_path), True, str(tmp_path), 3, True, True, True)
    assert args.qvalue == [0, 0.05] and args.processes == 1 and args.uniprot is None


def test_every_problem_is_reported(tmp_path, capsys):
    """ before PR 18, only the first. -q takes negative numbers as values, as before """
    with pytest.raises(SystemExit) as e:
        main.main([str(RMATS), 'missing.gtf', 'missing.fa', '-q', '-0.1', '-0.2', '-p', '0', '-o', str(tmp_path)])

    assert e.value.code == 2
    err = capsys.readouterr().err
    for message in ['qvalue lower bound must be between 0 and 1', 'qvalue upper bound must be between 0 and 1',
                    'qvalue lower bound must be smaller than upper bound', '-p/--processes must be at least 1',
                    'file not found: missing.gtf', 'file not found: missing.fa']:
        assert message in err


def test_validate_names_the_arguments(tmp_path):
    """ for the GUI, which shows each problem next to its field """
    (tmp_path / 'logfile.log').write_text('')
    args = argparse.Namespace(rmats_folder=str(tmp_path / 'none'), gtf_file=str(HUMAN_GTF), genome='missing.fa',
                              out=str(tmp_path), overwrite=False, uniprot=None, processes=2, splice_type={'SE'},
                              qvalue=[0, 1.5])

    assert [name for name, _ in main.validate(args)] == ['qvalue', 'processes', 'genome', 'rmats_folder', 'out']
    assert main.validate(argparse.Namespace(rmats_folder=str(RMATS), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME),
                                            out=str(tmp_path / 'new'), overwrite=False, uniprot=None, processes=1,
                                            splice_type={'SE'}, qvalue=[0, 1])) == []


def test_gui_command(monkeypatch):
    """ jcast gui, and --gui as before PR 18 (python -m jcast --gui), open the GUI; its exit status is the command's """
    opened = []
    monkeypatch.setattr(jcast.gui, 'launch_gui', lambda: opened.append(True) or 0)
    main.main(['gui'])
    main.main(['--gui'])
    assert opened == [True, True]

    monkeypatch.setattr(jcast.gui, 'launch_gui', lambda: 1)
    with pytest.raises(SystemExit) as e:
        main.main(['gui'])
    assert e.value.code == 1


def test_version(capsys):
    main.main(['--version'])
    assert capsys.readouterr().out == f'jcast {__version__}\n'


def test_help(capsys):
    """ no arguments show the commands; run -h the options, spelled as before PR 18 """
    main.main([])
    commands = capsys.readouterr().out
    assert ' run ' in commands and ' gui ' in commands

    main.main(['run', '-h'])
    options = capsys.readouterr().out
    for option in ['--out', '--overwrite', '--cache_dir', '--read', '--model', '--qvalue', '--splice_type', '--canonical',
                   '--uniprot', '--processes', '--verbose']:
        assert option in options
    assert '--mask' not in options and '--cache-dir' not in options


# ---- PR 19: progress and SIGTERM, for the GUI ----

UNIPROT = RMATS.parent / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'


def test_progress(tmp_path, cache_dir, capsys):
    """ run --progress writes each step as a line of JSON on stdout, the junctions about every 1 % """
    sigterm = signal.getsignal(signal.SIGTERM)
    main.main([str(a) for a in [RMATS, HUMAN_GTF, HUMAN_GENOME, '-s', 'RI', 'SE', '-r', 100000, '-o', tmp_path / 'out',
                                '--cache_dir', cache_dir, '--progress']])
    steps = [json.loads(line) for line in capsys.readouterr().out.splitlines()]

    assert [s['stage'] for s in steps if s['stage'] != 'junctions'] == ['rmats', 'gtf', 'genome', 'outputs', 'done']
    for splice_type, total in [('SE', 720), ('RI', 88)]:
        done = [s['done'] for s in steps if s.get('splice_type') == splice_type]
        assert done[0] == 0 and done[-1] == total and done == sorted(done) and len(done) < 110
    assert signal.getsignal(signal.SIGTERM) == sigterm    # the run's handler is removed after it


@pytest.mark.slow
@pytest.mark.parametrize('processes', [1, 2])
def test_sigterm_stops_the_run(tmp_path, cache_dir, processes):
    """ SIGTERM, the GUI's cancel, stops a run cleanly: exit status 143, the log says so, no worker process is left,
    and junctions.tsv holds the junctions done. The junction it interrupts is not taken for a failed one """
    out = tmp_path / 'out'
    run = subprocess.Popen([sys.executable, '-m', 'jcast', 'run', RMATS, HUMAN_GTF, HUMAN_GENOME, '-c', '--uniprot',
                            UNIPROT, '-p', str(processes), '-o', out, '--cache_dir', cache_dir, '--progress'],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in run.stdout:
        step = json.loads(line)
        if step['stage'] == 'junctions' and step['done'] > 0:
            run.send_signal(signal.SIGTERM)
            break
    _, err = run.communicate(timeout=120)

    assert run.returncode == 143
    assert 'Stopped by SIGTERM: the output is incomplete.' in err
    assert (out / 'logfile.log').read_text().splitlines()[-1].endswith('~WARNING~Stopped by SIGTERM: the output is '
                                                                      'incomplete. Run again with --overwrite to '
                                                                      'replace it.')
    with open(out / 'junctions.tsv') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    assert 0 < len(rows) < 1105 and not [r for r in rows if r['fate'] == 'exception']
    time.sleep(0.5)
    assert subprocess.run(['pgrep', '-f', str(out)], capture_output=True, text=True).stdout == ''
