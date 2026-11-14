# -*- coding: utf-8 -*-

"""
The Results page: the summary of a run, from its output folder, for this run or any earlier one, also one run on the
command line: junctions by fate and splice type (from junctions.tsv, by the code that writes the log's summary),
records per file, the junctions that failed, the warnings, and the file to search (reports/2026-11-13-gui-cli-scope.md
§4)
"""

import csv
import os
import re

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices, QFont
from PySide6.QtWidgets import (QFileDialog, QHBoxLayout, QHeaderView, QLabel, QPlainTextEdit, QPushButton, QSplitter,
                               QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from jcast import main, outcomes
from jcast.gui.tasks import Background
from jcast.uniprot import MERGED_FILE

LOG = 'logfile.log'
MODEL = 'model.pdf'

# lines of logfile.log
STAMP = r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~'
MESSAGE = re.compile(STAMP + r'(WARNING|ERROR)~(.*)$')
VERSION = re.compile(STAMP + r'INFO~JCAST (\S+)$')
CUTOFF = re.compile(STAMP + r'INFO~Translating junctions with a read count of at least (\d+) \((-[mr])\)')
TIME = re.compile(STAMP + r'INFO~Time: (.*) in all')
SUMMARY = re.compile(STAMP + r'INFO~Summary\.$')

FAILED_COLUMNS = ['junction_type', 'junction_id', 'gene_symbol', 'gene_id', 'fate_detail']

SEARCH_MERGED = f'<b>The file to search: {MERGED_FILE}</b>: the UniProt FASTA given with --uniprot, with the new T1 ' \
                f'and T1 alternative-terminus isoforms.'
SEARCH_PSQ = '<b>The file to search:</b> append psq_T1.fasta (and, for annotated alternative termini, ' \
             'psq_T1_alt_cterm.fasta and psq_T1_alt_nterm.fasta) to psq_canonical.fasta (with -c) or to your ' \
             'canonical database, and search the combined file. T2 and T4 sequences are frameshifted or truncated: ' \
             'use them with caution.'


def read_output(folder: str) -> dict:
    """
    What the Results page shows of an output folder. Run in the background: a whole genome's junctions.tsv and FASTA
    files take seconds to read

    :return: complete (the log has the summary), version, cutoff (e.g. '26 reads (-m)'), time, messages (the log's
        warnings and errors), fates (outcomes.fate_table; None without junctions.tsv), failed (the rows of the
        junctions that failed with an exception), records (outcomes.records_per_file), merged and model (whether those
        files are there), files (JCAST output files in the folder)
    """
    found = dict(folder=folder, complete=False, version=None, cutoff=None, time=None, messages=[], fates=None,
                 failed=[])
    log = os.path.join(folder, LOG)
    if os.path.isfile(log):
        with open(log, errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if m := MESSAGE.match(line):
                    found['messages'].append(f'{m.group(1)} -- {m.group(2)}')
                elif m := VERSION.match(line):
                    found['version'] = m.group(1)
                elif m := CUTOFF.match(line):
                    found['cutoff'] = f'{m.group(1)} read{"" if m.group(1) == "1" else "s"} ({m.group(2)})'
                elif m := TIME.match(line):
                    found['time'] = m.group(1)
                elif SUMMARY.match(line):
                    found['complete'] = True

    table = os.path.join(folder, outcomes.JUNCTIONS_FILE)
    if os.path.isfile(table):
        found['fates'] = outcomes.fate_table(outcomes.read_counts(table), None)
        with open(table, newline='') as fh:
            found['failed'] = [[row[c] for c in FAILED_COLUMNS] for row in csv.DictReader(fh, delimiter='\t')
                               if row['fate'] == 'exception']
    found['records'] = outcomes.records_per_file(folder)
    found['merged'] = os.path.isfile(os.path.join(folder, MERGED_FILE))
    found['model'] = os.path.isfile(os.path.join(folder, MODEL))
    found['files'] = main.existing_outputs(folder)
    return found


class ResultsPage(QWidget):

    def __init__(self, parent: QWidget = None):
        super().__init__(parent)
        self.folder = None
        self.shown = None           # read_output of the folder shown
        self._task = None

        self.path = QLabel('No output folder open.')
        self.path.setTextInteractionFlags(Qt.TextSelectableByMouse)
        choose = QPushButton('Open an output folder…')
        choose.clicked.connect(self._choose)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.search = QLabel()
        self.search.setWordWrap(True)

        self.fates = self._table()
        self.records = self._table(['File', 'Records'])
        self.records.itemDoubleClicked.connect(lambda item: self._open(self.records.item(item.row(), 0).text()))
        self.failed = self._table(['Splice type', 'rMATS ID', 'Gene', 'Gene ID', 'Exception'])
        self.messages = QPlainTextEdit()
        self.messages.setReadOnly(True)
        self.details = QTabWidget()
        self.details.addTab(self.messages, 'Warnings and errors')
        self.details.addTab(self.failed, 'Junctions that failed')

        self.buttons = {}
        row = QHBoxLayout()
        for name, label in (('', 'Open the folder'), (LOG, 'Open the log'), (outcomes.JUNCTIONS_FILE,
                            'Open junctions.tsv'), (MODEL, 'Open model.pdf')):
            button = QPushButton(label)
            button.clicked.connect(lambda checked=False, name=name: self._open(name))
            self.buttons[name] = button
            row.addWidget(button)
        row.addStretch()

        left = QWidget()
        column = QVBoxLayout(left)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(QLabel('Junctions by fate and splice type (each junction is a row of junctions.tsv):'))
        column.addWidget(self.fates)
        right = QWidget()
        column = QVBoxLayout(right)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(QLabel('Records per file (double-click to open):'))
        column.addWidget(self.records)
        tables = QSplitter()
        tables.addWidget(left)
        tables.addWidget(right)
        tables.setStretchFactor(0, 3)
        tables.setStretchFactor(1, 2)

        top = QHBoxLayout()
        top.addWidget(self.path, 1)
        top.addWidget(choose)
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.status)
        layout.addWidget(self.search)
        layout.addLayout(row)
        layout.addWidget(tables, 3)
        layout.addWidget(self.details, 1)
        self._show(None)

    @staticmethod
    def _table(headers: list = ()) -> QTableWidget:
        table = QTableWidget(0, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.verticalHeader().hide()
        table.verticalHeader().setDefaultSectionSize(table.fontMetrics().height() + 8)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        table.horizontalHeader().setStretchLastSection(True)
        return table

    def load(self, folder: str) -> None:
        """ Show the output folder's summary, read in the background """
        self.folder = folder
        self.path.setText(folder)
        self.status.setText('Reading the output…')
        self._task = Background(self, lambda: read_output(folder), lambda r: self._read(folder, r),
                                lambda e: self._read(folder, e))

    def _choose(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, 'A JCAST output folder', self.folder or '')
        if folder:
            self.load(folder)

    def _read(self, folder: str, result) -> None:
        if folder != self.folder:   # another folder was opened meanwhile
            return
        self._task = None
        if isinstance(result, Exception):
            self._show(None)
            self.path.setText(folder)
            self.status.setText(f'The output could not be read: {type(result).__name__}: {result}')
            return
        self._show(result)

    def _show(self, found) -> None:
        self.shown = found
        for name, button in self.buttons.items():
            button.setEnabled(found is not None and (name == '' or name in found['files']))
        if found is None:
            for table in (self.fates, self.records, self.failed):
                table.setRowCount(0)
            self.messages.clear()
            self.search.clear()
            return

        if not found['files']:
            state = 'No JCAST output in this folder'
        elif found['complete']:
            state = 'Complete'
        else:
            state = 'Incomplete: the run did not finish (its log has no summary)'
        settings = [f'JCAST {found["version"]}' if found['version'] else '',
                    f'read count cutoff {found["cutoff"]}' if found['cutoff'] else '',
                    f'{found["time"]} in all' if found['time'] else '']
        self.status.setText('; '.join([state] + [s for s in settings if s]) + '.')
        self.search.setText((SEARCH_MERGED if found['merged'] else SEARCH_PSQ) if found['files'] else '')

        fates = found['fates'] or [['fate']]
        self.fates.setColumnCount(len(fates[0]))
        self.fates.setHorizontalHeaderLabels(fates[0])
        self.fates.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.fates.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self._fill(self.fates, fates[1:], bold_last=True)
        self._fill(self.records, [[name, n] for name, n in found['records'].items()])
        self._fill(self.failed, found['failed'])
        self.messages.setPlainText('\n'.join(found['messages']) or 'None.')
        self.details.setTabText(0, f'Warnings and errors ({len(found["messages"])})')
        self.details.setTabText(1, f'Junctions that failed ({len(found["failed"])})')

    @staticmethod
    def _fill(table: QTableWidget, rows: list, bold_last: bool = False) -> None:
        table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                number = isinstance(value, int)
                item = QTableWidgetItem(f'{value:,}' if number else str(value))
                item.setTextAlignment((Qt.AlignRight if number else Qt.AlignLeft) | Qt.AlignVCenter)
                if bold_last and i == len(rows) - 1:
                    font = QFont(item.font())
                    font.setBold(True)
                    item.setFont(font)
                table.setItem(i, j, item)

    def _open(self, name: str) -> None:
        """ Open the folder, or a file in it, in the application the system chooses """
        QDesktopServices.openUrl(QUrl.fromLocalFile(os.path.join(self.folder, name) if name else self.folder))
