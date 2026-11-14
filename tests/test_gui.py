# -*- coding: utf-8 -*-

"""
The GUI (PR 13; reports/2026-11-13-gui-cli-scope.md §4), without a display: Qt's offscreen platform. Skipped without
PySide6 or pytest-qt: pip install '.[test-gui]'
"""

import argparse
import math
import os
from types import SimpleNamespace

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
pytest.importorskip('PySide6.QtWidgets')
pytest.importorskip('pytestqt')

from PySide6.QtGui import QDesktopServices                          # noqa: E402
from PySide6.QtWidgets import QMessageBox                           # noqa: E402

from jcast import main, model, outcomes                             # noqa: E402
from jcast.gui import command                                       # noqa: E402
from jcast.gui.main_window import MainWindow                        # noqa: E402
from jcast.gui.results_page import SEARCH_MERGED, SEARCH_PSQ, read_output     # noqa: E402
from jcast.gui.run_page import RunPage                              # noqa: E402
from jcast.outcomes import COLUMNS, SPLICE_TYPES                    # noqa: E402

from tests.conftest import DATA, HUMAN_GENOME, HUMAN_GTF, RMATS     # noqa: E402
from tests.test_command_line import UNIPROT, _folder_with, _parsed  # noqa: E402

STOMACH = DATA / 'model' / 'encode_human_stomach'       # a degenerate fit
EMPTY = DATA / 'model' / 'encode_mouse_fetalheart'      # header-only rMATS files


@pytest.fixture
def window(qtbot):
    w = MainWindow()
    qtbot.addWidget(w)
    w.show()
    yield w
    if w.run.running():     # closing the window would ask whether to stop it
        w.run.process.kill()
        w.run.process.waitForFinished()


def _answer(monkeypatch, answer):
    """ QMessageBox.question gives this answer, and records its questions """
    asked = []
    monkeypatch.setattr(QMessageBox, 'question', lambda *a, **k: asked.append(a[2]) or answer)
    return asked


def _fill(window, out, rmats=RMATS, uniprot=None, processes=1, splice_types=SPLICE_TYPES, **options):
    """ the Setup page's fields, checked at once """
    setup = window.setup
    setup.rmats_folder.setText(str(rmats))
    setup.gtf_file.setText(str(HUMAN_GTF))
    setup.genome.setText(str(HUMAN_GENOME))
    setup.out.setText(str(out))
    setup.uniprot.setText(str(uniprot or ''))
    setup.refresh()
    for t, box in setup.splice_types.items():
        box.setChecked(t in splice_types)
    setup.processes.setValue(processes)
    for name, value in options.items():
        widget = getattr(setup, name)
        if name == 'cache_dir':
            widget.setText(str(value))
        elif name in ('q_lo', 'q_hi'):
            widget.setValue(value)
        else:
            widget.setChecked(value)
    setup.refresh()


def _errors(window) -> dict:
    return {name: label.text() for name, label in window.setup.errors.items() if not label.isHidden()}


def test_four_pages(window):
    assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == ['Setup', 'Model', 'Run', 'Results']


# ---- Setup ----

def test_the_settings_are_a_command_line(window, tmp_path, monkeypatch, cache_dir):
    """ the Setup page's settings are the arguments run_jcast takes; the Run page shows them as a command line, which
    the command line parses into the same arguments """
    _fill(window, tmp_path / 'out', uniprot=UNIPROT, processes=2, splice_types=['RI', 'A5SS'], q_hi=0.05,
          canonical=True, verbose=True, cache_dir=cache_dir)
    args = window.arguments()

    assert args == argparse.Namespace(rmats_folder=str(RMATS), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME),
                                      out=str(tmp_path / 'out'), overwrite=False, cache_dir=str(cache_dir), read=1,
                                      model=False, canonical=True, uniprot=str(UNIPROT), processes=2, verbose=True,
                                      splice_type=['RI', 'A5SS'], qvalue=[0, 0.05])
    assert window.setup.problems == [] and _errors(window) == {}
    assert window.run.command.toPlainText() == command.command_line(args)
    assert _parsed(monkeypatch, 'run', *command.arguments(args)) == args


def test_problems_next_to_their_fields(window, tmp_path):
    """ the command line's checks (main.validate), each under its field; empty fields are shown once a run is
    started """
    _fill(window, tmp_path / 'out', q_lo=0.5, q_hi=0.1)
    window.setup.gtf_file.setText(str(tmp_path / 'missing.gtf'))
    window.setup.genome.setText('')
    window.setup.refresh()

    assert _errors(window) == {'gtf_file': f'file not found: {tmp_path / "missing.gtf"}',
                               'qvalue': 'qvalue lower bound must be smaller than upper bound'}
    assert ('genome', 'required') in window.setup.problems
    assert window.run.command.toPlainText().startswith('The command appears here once')

    window.tabs.setCurrentWidget(window.run)
    window.run.start.click()
    assert window.tabs.currentWidget() is window.setup and window.run.process is None
    assert _errors(window)['genome'] == 'required'


def test_the_splice_types_of_the_folder(window, tmp_path):
    """ only the splice types whose rMATS file is in the folder can be selected: -s, as the command line needs """
    _fill(window, tmp_path / 'out', rmats=_folder_with(tmp_path, 'SE', 'RI'))

    assert [t for t, box in window.setup.splice_types.items() if box.isEnabled()] == ['SE', 'RI']
    assert window.arguments().splice_type == ['SE', 'RI'] and window.setup.problems == []
    assert window.run.command.toPlainText().endswith(' -s SE RI')


def test_processes_need_uniprot(window, tmp_path):
    _fill(window, tmp_path / 'out')
    assert not window.setup.processes.isEnabled() and window.arguments().processes == 1

    _fill(window, tmp_path / 'out', uniprot=UNIPROT, processes=3)
    assert window.setup.processes.isEnabled() and window.arguments().processes == 3


def test_an_output_folder_holding_output(window, tmp_path):
    """ is not a problem: starting the run asks before replacing it """
    (tmp_path / 'logfile.log').write_text('')
    _fill(window, tmp_path)

    assert window.setup.problems == []
    assert _errors(window)['out'].startswith('This folder holds the output of an earlier run: starting the run asks')


# ---- Model ----

def _fitted(window, qtbot):
    window.tabs.setCurrentWidget(window.model)
    qtbot.waitUntil(lambda: window.model._fit is not None or 'could not' in window.model.status.text(), timeout=30000)
    return window.model


def _table(page) -> list:
    return [[page.table.item(i, j).text() for j in range(page.table.columnCount())]
            for i in range(page.table.rowCount())]


def test_the_model_and_its_cutoff(window, qtbot, tmp_path):
    """ the model of the folder, fitted in the background when the page is shown; its cutoff (-m), or another set in
    the number field or on the plot (-r); and the junctions each would keep """
    _fill(window, tmp_path / 'out')
    page = _fitted(window, qtbot)

    assert page.status.text().endswith('its cutoff is 26 reads.') and '1,105 junctions' in page.status.text()
    assert page.use_model.isEnabled() and page.options() == (False, 1)
    page.use_model.setChecked(True)
    assert page.options() == (True, 1) and window.arguments().model
    assert window.run.command.toPlainText().endswith(' -m')
    assert page._cutoff[0].get_xdata()[0] == pytest.approx(math.log(27))

    page.update_table()
    expected = model.filter_counts(model.preview(RMATS)[0], 26, (0, 1))
    assert _table(page) == [[t] + [f'{n:,}' for n in expected.loc[t]] for t in SPLICE_TYPES] + \
           [['total'] + [f'{n:,}' for n in expected.sum()]]

    # a click on the plot: -r there
    page._clicked(SimpleNamespace(name='button_press_event', xdata=math.log(11), button=1))
    assert page.options() == (False, 10) and page.read.value() == 10
    assert page._cutoff[0].get_xdata()[0] == pytest.approx(math.log(11))
    assert window.run.command.toPlainText().endswith(' -r 10')
    # moving the mouse without a button pressed does nothing
    page._clicked(SimpleNamespace(name='motion_notify_event', xdata=math.log(101), button=None))
    assert page.read.value() == 10

    # the table: the splice types selected, and -q
    window.setup.splice_types['SE'].setChecked(False)
    window.setup.q_hi.setValue(0.05)
    window.setup.refresh()
    page.update_table()
    expected = model.filter_counts(model.preview(RMATS)[0], 10, (0, 0.05))
    assert [row[0] for row in _table(page)] == ['MXE', 'RI', 'A5SS', 'A3SS', 'total']
    assert _table(page)[1] == ['RI'] + [f'{n:,}' for n in expected.loc['RI']]

    # the number field
    page.use_model.setChecked(True)
    page.read.setValue(40)
    assert page.options() == (False, 40)


@pytest.mark.parametrize('folder, reason', [(STOMACH, 'the fit was degenerate'),
                                            (EMPTY, 'the fit failed: the rMATS files read hold no junctions')])
def test_a_folder_without_a_fit(window, qtbot, tmp_path, folder, reason):
    """ no -m: why, and the histogram with the -r cutoff, if there are junctions """
    _fill(window, tmp_path / 'out', rmats=folder)
    page = _fitted(window, qtbot)

    assert reason in page.status.text() and 'A run with -m would fall back to -r.' in page.status.text()
    assert not page.use_model.isEnabled() and page.options() == (False, 1)
    assert (page._cutoff is not None) == (folder == STOMACH)
    page.update_table()
    assert len(_table(page)) == (6 if folder == STOMACH else 0)


def test_another_folder_is_fitted_again(window, qtbot, tmp_path):
    """ the model's cutoff chosen for one folder is not kept for another: -r until it is chosen again """
    _fill(window, tmp_path / 'out')
    page = _fitted(window, qtbot)
    page.use_model.setChecked(True)

    window.setup.rmats_folder.setText(str(STOMACH))
    window.setup.refresh()
    assert page.options() == (False, 1) and not page.use_model.isEnabled()
    qtbot.waitUntil(lambda: page._fit is not None and page._fit[0] == str(STOMACH), timeout=30000)
    assert 'degenerate' in page.status.text()

    # a folder that is no longer one, while its model is fitted: the fit is dropped when it ends
    window.setup.rmats_folder.setText(str(RMATS))
    window.setup.refresh()
    assert page.status.text().startswith('Fitting') and not page.busy.isHidden()
    window.setup.rmats_folder.setText('')
    window.setup.refresh()
    assert page.status.text() == 'Choose an rMATS folder on the Setup page.' and page.busy.isHidden()
    qtbot.wait(500)
    assert page._fit is None and page.table.rowCount() == 0


# ---- Run and Results ----

def test_starting_asks_to_replace_earlier_output(window, tmp_path, monkeypatch):
    launched = []
    monkeypatch.setattr(RunPage, 'launch', lambda self, args: launched.append(args))
    (tmp_path / 'logfile.log').write_text('')
    _fill(window, tmp_path)

    asked = _answer(monkeypatch, QMessageBox.No)
    window.run.start.click()
    assert launched == [] and 'holds the output of an earlier JCAST run: logfile.log.' in asked[0]

    _answer(monkeypatch, QMessageBox.Yes)
    window.run.start.click()
    assert launched[0].overwrite and not window.arguments().overwrite


@pytest.mark.slow
def test_a_run(window, qtbot, tmp_path, monkeypatch, cache_dir):
    """ a run from the GUI is the run of its command line: the same output. Its progress, then its results """
    _fill(window, tmp_path / 'gui', uniprot=UNIPROT, processes=2, splice_types=['RI', 'A5SS'], canonical=True,
          cache_dir=cache_dir)
    with qtbot.waitSignal(window.run.finished, timeout=120000) as finished:
        window.run.start.click()
    assert finished.args == [str(tmp_path / 'gui'), 0]
    assert window.run.status.text().startswith('Finished in')
    assert {t: (b.value(), b.maximum()) for t, b in window.run.bars.items()} == {'RI': (88, 88), 'A5SS': (71, 71)}

    qtbot.waitUntil(lambda: window.results.shown is not None, timeout=30000)
    assert window.tabs.currentWidget() is window.results
    shown = window.results.shown
    assert shown['complete'] and shown['cutoff'] == '1 read (-r)' and shown['merged']
    assert shown['fates'] == outcomes.fate_table(outcomes.read_counts(tmp_path / 'gui' / 'junctions.tsv'), None)
    assert window.results.search.text() == SEARCH_MERGED

    # the command line, with another output folder
    args = window.arguments()
    args.out = str(tmp_path / 'cli')
    main.main(['run'] + command.arguments(args))
    names = sorted(os.listdir(tmp_path / 'gui'))
    assert names == sorted(os.listdir(tmp_path / 'cli'))
    for name in names:
        if name != 'logfile.log':
            assert (tmp_path / 'gui' / name).read_bytes() == (tmp_path / 'cli' / name).read_bytes(), name


def _started(window, qtbot, tmp_path, cache_dir):
    """ a run of every chr15 junction, under way """
    _fill(window, tmp_path / 'out', uniprot=UNIPROT, processes=2, cache_dir=cache_dir)
    window.run.start.click()
    qtbot.waitUntil(lambda: any(bar.value() > 0 for bar in window.run.bars.values()), timeout=60000)


@pytest.mark.slow
def test_stop(window, qtbot, tmp_path, monkeypatch, cache_dir):
    """ Stop: the run ends as on SIGTERM (exit status 143), its output marked incomplete, and can be deleted """
    _started(window, qtbot, tmp_path, cache_dir)
    with qtbot.waitSignal(window.run.finished, timeout=30000) as finished:
        window.run.stop_run(ask=False)

    assert finished.args == [str(tmp_path / 'out'), main.STOPPED]
    assert 'the output in' in window.run.status.text() and window.run.status.text().startswith('Stopped after')
    assert 'WARNING -- Stopped by SIGTERM' in window.run.messages.toPlainText()
    qtbot.waitUntil(lambda: window.results.shown is not None, timeout=30000)
    assert window.tabs.currentWidget() is not window.results
    assert not window.results.shown['complete'] and window.results.status.text().startswith('Incomplete')

    assert not window.run.delete.isHidden()
    _answer(monkeypatch, QMessageBox.Yes)
    window.run.delete.click()
    assert main.existing_outputs(str(tmp_path / 'out')) == [] and window.run.delete.isHidden()


@pytest.mark.slow
def test_quitting_stops_the_run(window, qtbot, tmp_path, monkeypatch, cache_dir):
    _started(window, qtbot, tmp_path, cache_dir)
    asked = _answer(monkeypatch, QMessageBox.Yes)
    window.close()

    assert asked == ['A run is in progress. Its output will be incomplete.']
    assert not window.run.running() and window.run.process.exitCode() == main.STOPPED


def _output(folder, merged=False):
    """ the output of a run on the command line: one T1 junction, one that failed, one skipped """
    folder.mkdir()
    rows = [dict(junction_type='SE', junction_id='1', fate='T1', gene_symbol='A', gene_id='G1'),
            dict(junction_type='SE', junction_id='2', fate='exception', fate_detail="KeyError: 'chr99'",
                 gene_symbol='B', gene_id='G2'),
            dict(junction_type='RI', junction_id='3', fate='skipped_low_count', gene_symbol='C', gene_id='G3')]
    (folder / 'junctions.tsv').write_text('\n'.join(['\t'.join(COLUMNS)] + ['\t'.join(r.get(c, '') for c in COLUMNS)
                                                                            for r in rows]) + '\n')
    (folder / 'psq_T1.fasta').write_text('>a\nMAAA\n>b\nMBBB\n')
    stamp = '2026-11-14 10:00:00,000~'
    (folder / 'logfile.log').write_text('\n'.join([
        stamp + 'INFO~JCAST 0.5.0', stamp + 'INFO~Arguments: Namespace()',
        stamp + 'INFO~Translating junctions with a read count of at least 26 (-m): the smaller of …',
        stamp + "ERROR~Junction SE 2 of gene B G2 failed:", 'Traceback (most recent call last):', "KeyError: 'chr99'",
        stamp + 'WARNING~Junction SE 2 of gene B G2 failed: KeyError. The run continues.',
        stamp + 'INFO~Summary.', 'Junctions by fate …', stamp + 'INFO~Time: 1 min 05 s in all, of which 50 s …',
        stamp + 'WARNING~1 junction(s) failed with an exception.']) + '\n')
    if merged:
        (folder / 'jcast_uniprot_merged.fasta').write_text('>sp|P1|A\nMAAA\n')
    return folder


def test_results_of_any_output_folder(window, qtbot, tmp_path, monkeypatch):
    """ the summary of a run made on the command line: from junctions.tsv, by the code that writes the log's
    summary; the junctions that failed; the log's warnings and errors; the file to search """
    folder = _output(tmp_path / 'out')
    found = read_output(str(folder))
    assert (found['complete'], found['version'], found['cutoff'], found['time']) == \
           (True, '0.5.0', '26 reads (-m)', '1 min 05 s')
    assert found['fates'] == outcomes.fate_table(outcomes.read_counts(folder / 'junctions.tsv'), None)
    assert found['failed'] == [['SE', '2', 'B', 'G2', "KeyError: 'chr99'"]]
    assert found['messages'] == ['ERROR -- Junction SE 2 of gene B G2 failed:',
                                 'WARNING -- Junction SE 2 of gene B G2 failed: KeyError. The run continues.',
                                 'WARNING -- 1 junction(s) failed with an exception.']
    assert found['records'] == {'psq_T1.fasta': 2} and not found['merged'] and not found['model']

    page = window.results
    page.load(str(folder))
    qtbot.waitUntil(lambda: page.shown is not None)
    assert page.status.text() == 'Complete; JCAST 0.5.0; read count cutoff 26 reads (-m); 1 min 05 s in all.'
    assert page.search.text() == SEARCH_PSQ
    assert page.fates.rowCount() == 11 and page.fates.item(10, 0).text() == 'total'
    assert [page.fates.item(0, j).text() for j in range(4)] == ['T1', '1', '0', '1']
    assert page.failed.rowCount() == 1 and page.details.tabText(0) == 'Warnings and errors (3)'
    assert {name for name, b in page.buttons.items() if b.isEnabled()} == {'', 'logfile.log', 'junctions.tsv'}

    opened = []
    monkeypatch.setattr(QDesktopServices, 'openUrl', lambda url: opened.append(url.toLocalFile()))
    page.buttons['logfile.log'].click()
    assert opened == [str(folder / 'logfile.log')]

    page.load(str(_output(tmp_path / 'merged', merged=True)))
    qtbot.waitUntil(lambda: page.shown is not None and page.shown['folder'] == str(tmp_path / 'merged'))
    assert page.search.text() == SEARCH_MERGED

    page.load(str(tmp_path))
    qtbot.waitUntil(lambda: page.shown is not None and page.shown['folder'] == str(tmp_path))
    assert page.status.text() == 'No JCAST output in this folder.' and page.fates.rowCount() == 0
