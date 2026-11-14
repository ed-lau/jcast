# -*- coding: utf-8 -*-

"""
The Run page: the run's command line, run as `python -P -m jcast run … --progress` in a child process (a QProcess),
its progress and warnings, and a stop. The child process forks the -p workers itself, runs no Qt, and a crash ends the
run only; Stop sends SIGTERM, which the run handles: the workers ended, the log says so, exit status 143
(reports/2026-11-13-gui-cli-scope.md §3)
"""

import copy
import json
import os
import time

from PySide6.QtCore import QElapsedTimer, QProcess, QProcessEnvironment, Qt, QTimer, Signal
from PySide6.QtGui import QFontDatabase, QGuiApplication
from PySide6.QtWidgets import (QGridLayout, QHBoxLayout, QLabel, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
                               QVBoxLayout, QWidget)

from jcast import main
from jcast.gui.command import command_line, process
from jcast.outcomes import SPLICE_TYPES

# the steps of a run, as run --progress writes them
STAGES = {'rmats': 'Reading the rMATS files', 'gtf': 'Reading the GTF', 'genome': 'Reading the genome',
          'uniprot': 'Reading the UniProt FASTA', 'model': 'Fitting the -m model', 'junctions': 'Translating junctions',
          'outputs': 'Writing the outputs', 'done': 'Done'}

BARS_TIP = 'The junctions done, of each splice type, in file order; with -p, the junctions whose results have been ' \
           'written. A junction skipped by the read count or FDR filters is done at once.'

# after Stop, the time the run has to end its workers and log the stop before it is killed (as on Windows, where
# QProcess.terminate does not stop a console program)
KILL_AFTER_MS = 15000


class RunPage(QWidget):
    """ finished is emitted when a run ends, with its output folder and exit status (-1 if it crashed) """

    finished = Signal(str, int)
    setup_needed = Signal()         # Start, while the Setup page has problems

    def __init__(self, parent: QWidget = None):
        super().__init__(parent)
        self._args, self._problems = None, []   # the settings of the Setup and Model pages
        self.args = None                        # those of the run started last
        self.process = None
        self._stopping = False
        self._started = None                    # when it started (time.time())
        self._stage = ''

        self.command = QPlainTextEdit()
        self.command.setReadOnly(True)
        self.command.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.command.setMaximumHeight(90)
        self.command.setStatusTip('The command line of the run. The GUI runs it as python -P -m jcast run … '
                                  '--progress, with its own Python; the output and the log are those of the command '
                                  'line.')
        self.copy = QPushButton('Copy')
        self.copy.setStatusTip('Copy the command, to run it in a terminal or keep it with the results.')
        self.copy.clicked.connect(lambda: QGuiApplication.clipboard().setText(self.command.toPlainText()))
        self.start = QPushButton('Start')
        self.start.setStatusTip('Start the run, in a separate process. An output folder holding the output of an '
                                'earlier run is replaced only after a question (--overwrite).')
        self.start.clicked.connect(self._start)
        self.stop = QPushButton('Stop')
        self.stop.setStatusTip('Stop the run: it ends its worker processes and its log says that it stopped (exit '
                               'status 143). Its output is then incomplete.')
        self.stop.setEnabled(False)
        self.stop.clicked.connect(lambda: self.stop_run(ask=True))
        self.status = QLabel('Not started.')
        self.status.setWordWrap(True)
        self.delete = QPushButton('Delete its output')
        self.delete.setStatusTip('Delete the files the run wrote in the output folder; the folder and any other file '
                                 'in it are kept.')
        self.delete.hide()
        self.delete.clicked.connect(self._delete_output)
        self.bars, self.counts = {}, {}             # each splice type's progress bar, and its count
        self._bars = QGridLayout()
        self._bars.setColumnStretch(1, 1)
        self.last = QLabel()
        self.last.setTextFormat(Qt.PlainText)
        self.last.setStatusTip(BARS_TIP)
        self.messages = QPlainTextEdit()
        self.messages.setReadOnly(True)
        self.messages.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.messages.setStatusTip('What the run writes on the console: its warnings and errors, as they come. The '
                                   'log (logfile.log) has them too.')

        buttons = QHBoxLayout()
        for button in (self.start, self.stop, self.copy):
            buttons.addWidget(button)
        buttons.addStretch()
        status = QHBoxLayout()
        status.addWidget(self.status, 1)
        status.addWidget(self.delete)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel('The command (the GUI runs it with this Python, as python -m jcast run):'))
        layout.addWidget(self.command)
        layout.addLayout(buttons)
        layout.addLayout(status)
        layout.addLayout(self._bars)
        layout.addWidget(self.last)
        layout.addWidget(QLabel('Warnings and errors, as the run writes them on the console:'))
        layout.addWidget(self.messages, 1)

        self._clock = QElapsedTimer()
        self._tick = QTimer(self)
        self._tick.setInterval(1000)
        self._tick.timeout.connect(self._show_stage)
        self._kill = QTimer(self)
        self._kill.setSingleShot(True)
        self._kill.timeout.connect(self._kill_run)

    def set_arguments(self, args, problems: list) -> None:
        """ The settings of the Setup and Model pages, and the Setup page's problems """
        self._args, self._problems = args, problems
        if not self.running():
            self.command.setPlainText(command_line(args) if not problems else
                                      'The command appears here once the settings on the Setup page are complete.')
            self.copy.setEnabled(not problems)

    def running(self) -> bool:
        return self.process is not None and self.process.state() != QProcess.NotRunning

    # ---- start ----

    def _start(self) -> None:
        """ Start a run with the settings, after asking to replace the output of an earlier run in the folder """
        if self._problems:
            self.setup_needed.emit()
            return
        args = copy.copy(self._args)
        found = main.existing_outputs(args.out)
        if found:
            answer = QMessageBox.question(
                self, 'Replace the earlier output?',
                f'{args.out} holds the output of an earlier JCAST run: {", ".join(found[:4])}'
                f'{", …" if len(found) > 4 else ""}.\n\nReplace it (--overwrite)?')
            if answer != QMessageBox.Yes:
                return
            args.overwrite = True
        self.launch(args)

    def launch(self, args) -> None:
        """ Run `python -P -m jcast run <arguments> --progress` in a child process (command.process) """
        self.args, self._stopping, self._started = args, False, time.time()
        self.command.setPlainText(command_line(args))
        self.messages.clear()
        self.delete.hide()
        self._show_bars([t for t in SPLICE_TYPES if t in args.splice_type])

        program, arguments, pythonpath = process(args)
        environment = QProcessEnvironment.systemEnvironment()
        if pythonpath:
            environment.insert('PYTHONPATH', pythonpath)
        if self.process is not None:
            self.process.deleteLater()
        self.process = QProcess(self)
        self.process.setProgram(program)
        self.process.setArguments(arguments)
        self.process.setProcessEnvironment(environment)
        self.process.readyReadStandardOutput.connect(self._read_progress)
        self.process.readyReadStandardError.connect(self._read_messages)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)
        self._stage = 'Starting'
        self._clock.start()
        self._tick.start()
        self.start.setEnabled(False)
        self.stop.setEnabled(True)
        self.copy.setEnabled(True)
        self.process.start()
        self._show_stage()

    def _show_bars(self, splice_types: list) -> None:
        """ A progress bar for each splice type, and its count beside it: the macOS style draws no text in a bar """
        while self._bars.count():
            self._bars.takeAt(0).widget().deleteLater()
        self.bars, self.counts = {}, {}
        for row, t in enumerate(splice_types):
            self.bars[t], self.counts[t] = QProgressBar(), QLabel('waiting')
            self.bars[t].setTextVisible(False)
            self.bars[t].setValue(0)
            for column, widget in enumerate((QLabel(t), self.bars[t], self.counts[t])):
                widget.setStatusTip(BARS_TIP)
                self._bars.addWidget(widget, row, column)
        self.last.clear()

    # ---- the run's output ----

    def _read_progress(self) -> None:
        """ run --progress: a line of JSON for each step """
        while self.process.canReadLine():
            line = bytes(self.process.readLine()).decode('utf-8', 'replace').rstrip()
            try:
                step = json.loads(line)
            except ValueError:
                self.messages.appendPlainText(line)
                continue
            self._stage = STAGES.get(step['stage'], step['stage'])
            if step['stage'] == 'junctions' and step['splice_type'] in self.bars:
                splice_type, done, total = step['splice_type'], step['done'], step['total']
                self.bars[splice_type].setMaximum(max(total, 1))
                self.bars[splice_type].setValue(done if total else 1)
                self.counts[splice_type].setText(f'{done:,} of {total:,} ({done / total:.0%})' if total else
                                                 'no junctions')
                if 'gene' in step:
                    self.last.setText(f'Last junction done: {splice_type} {step["junction"]}, gene {step["gene"]} '
                                      f'({step["gene_id"]})')
            self._show_stage()

    def _read_messages(self) -> None:
        """ The warnings and errors the run writes on the console """
        text = bytes(self.process.readAllStandardError()).decode('utf-8', 'replace')
        cursor = self.messages.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        cursor.insertText(text)
        self.messages.setTextCursor(cursor)

    def _show_stage(self) -> None:
        if self.running():
            self.status.setText(f'{self._stage}… ({main._duration(self._clock.elapsed() / 1000)})')

    # ---- the end ----

    def _finished(self, code: int, status: QProcess.ExitStatus) -> None:
        self._read_progress()
        self._read_messages()
        self._tick.stop()
        self._kill.stop()
        self.start.setEnabled(True)
        self.stop.setEnabled(False)
        took = main._duration(self._clock.elapsed() / 1000)
        out = self.args.out
        if status == QProcess.NormalExit and code == 0:
            self.status.setText(f'Finished in {took}. The Results page shows its summary.')
        elif self._stopping:
            self.status.setText(f'Stopped after {took}: the output in {out} is incomplete.')
        elif status == QProcess.CrashExit:
            self.status.setText(f'The run ended unexpectedly after {took}: the output in {out} is incomplete.')
        else:
            self.status.setText(f'The run failed after {took} (exit status {code}); see the errors below.')
        self.delete.setVisible(not (status == QProcess.NormalExit and code == 0) and self._wrote())
        self.finished.emit(out, code if status == QProcess.NormalExit else -1)

    def _wrote(self) -> bool:
        """ Whether the run wrote output: not if it ended before its log was opened (the output of an earlier run,
        which a run with --overwrite replaces once started, is still there then) """
        log = os.path.join(self.args.out, 'logfile.log')
        return os.path.isfile(log) and os.path.getmtime(log) >= self._started - 1

    def _error(self, error: QProcess.ProcessError) -> None:
        if error == QProcess.FailedToStart:     # finished is not emitted
            self._tick.stop()
            self.start.setEnabled(True)
            self.stop.setEnabled(False)
            self.status.setText(f'The run could not start: {self.process.errorString()}')

    def stop_run(self, ask: bool) -> None:
        """ Stop the run: SIGTERM, then, if it has not ended after KILL_AFTER_MS, kill it """
        if not self.running():
            return
        if ask and QMessageBox.question(self, 'Stop the run?', 'Its output will be incomplete.') != QMessageBox.Yes:
            return
        self._stopping = True
        self._stage = 'Stopping'
        self._show_stage()
        self.process.terminate()
        self._kill.start(KILL_AFTER_MS)

    def _kill_run(self) -> None:
        if self.running():
            self.process.kill()

    def wait(self) -> None:
        """ Wait for a stopped run to end (quitting) """
        if self.running():
            self.process.waitForFinished(KILL_AFTER_MS + 5000)

    def _delete_output(self) -> None:
        """ Delete the output of a run that did not finish: the files a run writes, not the folder """
        out = self.args.out
        found = main.existing_outputs(out)
        if not found or QMessageBox.question(
                self, 'Delete the output?',
                f'Delete the {len(found)} files the run wrote in {out}?') != QMessageBox.Yes:
            return
        for name in found:
            os.remove(os.path.join(out, name))
        self.delete.hide()
        self.status.setText(self.status.text() + ' Its output was deleted.')
