# -*- coding: utf-8 -*-

"""
The Setup page: the run's inputs, output and options, with the command line's checks (main.validate) shown next to the
fields. The read count cutoff is set on the Model page (reports/2026-11-13-gui-cli-scope.md §4). The main window keeps
the page's settings between sessions (PR 26)
"""

import argparse
import os

from PySide6.QtCore import QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
                               QLineEdit, QMessageBox, QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget)

from jcast import main
from jcast.gui.command import option_tip
from jcast.junctions import rmats_path
from jcast.outcomes import SPLICE_TYPES

# the fields a run needs
REQUIRED = ['rmats_folder', 'gtf_file', 'genome', 'out']

# the page's settings (SetupPage.values), as a new page has them; the main window keeps them between sessions
PATHS = ['rmats_folder', 'gtf_file', 'genome', 'out', 'cache_dir', 'uniprot']
DEFAULTS = {**{name: '' for name in PATHS}, 'splice_types': ','.join(SPLICE_TYPES), 'q_lo': 0.0, 'q_hi': 1.0,
            'canonical': False, 'processes': 1, 'verbose': False}

FASTA = 'FASTA (*.fa *.fa.gz *.fasta *.fasta.gz *.fna *.fna.gz);;All files (*)'
GTF = 'GTF (*.gtf);;All files (*)'

ERROR_STYLE = 'color: #d32f2f'


class PathField(QWidget):
    """ A path, typed or chosen with Browse… """

    changed = Signal()

    def __init__(self, folder: bool, filter_: str = '', placeholder: str = '', parent: QWidget = None):
        super().__init__(parent)
        self.folder, self.filter = folder, filter_
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.textChanged.connect(self.changed)
        browse = QPushButton('Browse…')
        browse.clicked.connect(self._browse)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.edit, 1)
        layout.addWidget(browse)

    def text(self) -> str:
        return self.edit.text().strip()

    def setText(self, text: str) -> None:
        self.edit.setText(text)

    def _browse(self) -> None:
        if self.folder:
            path = QFileDialog.getExistingDirectory(self, '', self.text())
        else:
            path, _ = QFileDialog.getOpenFileName(self, '', os.path.dirname(self.text()), self.filter)
        if path:
            self.edit.setText(path)


class SetupPage(QWidget):
    """ The settings of a run but the read count cutoff; changed is emitted when they change, once checked """

    changed = Signal()

    def __init__(self, parent: QWidget = None):
        super().__init__(parent)
        self.problems = []              # the problems found by the last check, as (argument, message)
        self._show_required = False     # 'required' is shown once a run was started without the field
        self.errors = {}                # argument -> the label showing its problems

        self.rmats_folder = PathField(folder=True)
        self.gtf_file = PathField(folder=False, filter_=GTF)
        self.genome = PathField(folder=False, filter_=FASTA)
        self.out = PathField(folder=True)
        self.cache_dir = PathField(folder=True, placeholder='the user cache directory')
        self.uniprot = PathField(folder=False, filter_=FASTA, placeholder='none: canonicals from the UniProt API')

        self.splice_types = {t: QCheckBox(t) for t in SPLICE_TYPES}
        for box in self.splice_types.values():
            box.setChecked(True)
        self.q_lo, self.q_hi = QDoubleSpinBox(), QDoubleSpinBox()
        for spin, value in ((self.q_lo, 0), (self.q_hi, 1)):
            spin.setStatusTip(option_tip('qvalue'))
            spin.setRange(0, 1)
            spin.setDecimals(4)
            spin.setSingleStep(0.01)
            spin.setValue(value)
        self.canonical = QCheckBox("Write each gene's canonical protein to psq_canonical.fasta (-c)")
        self.processes = QSpinBox()
        self.processes.setRange(1, os.cpu_count() or 1)
        self.processes.setMaximumWidth(80)
        self.verbose = QCheckBox("Log each junction's steps and sequences (-v)")

        inputs = QFormLayout()
        self._row(inputs, 'rMATS folder', self.rmats_folder, 'rmats_folder')
        self._row(inputs, 'GTF file', self.gtf_file, 'gtf_file')
        self._row(inputs, 'Genome FASTA', self.genome, 'genome')
        types = QHBoxLayout()
        for box in self.splice_types.values():
            types.addWidget(box)
        types.addStretch()
        self._row(inputs, 'Splice types (-s)', types, 'splice_type')

        output = QFormLayout()
        self._row(output, 'Output folder (-o)', self.out, 'out')
        self._row(output, 'Cache folder (--cache_dir)', self.cache_dir, 'cache_dir')

        filters = QFormLayout()
        q = QHBoxLayout()
        q.addWidget(self.q_lo)
        q.addWidget(QLabel('to'))
        q.addWidget(self.q_hi)
        q.addStretch()
        self._row(filters, 'FDR range (-q)', q, 'qvalue')
        filters.addRow('', QLabel('The read count cutoff (-r or -m) is set on the Model page.'))

        canonicals = QFormLayout()
        canonicals.addRow('', self.canonical)
        self.canonical.setStatusTip(option_tip('canonical'))
        self._row(canonicals, 'UniProt FASTA (--uniprot)', self.uniprot, 'uniprot')

        the_run = QFormLayout()
        self._row(the_run, 'Processes (-p)', self.processes, 'processes')
        the_run.addRow('', self.verbose)
        self.verbose.setStatusTip(option_tip('verbose'))

        self.restore = QPushButton('Restore the defaults')
        self.restore.setStatusTip('Clear every path and set every option to its default. The settings of this page '
                                  'are kept between sessions; the read count cutoff is not.')
        self.restore.clicked.connect(self._restore_defaults)
        restore = QHBoxLayout()
        restore.addStretch()
        restore.addWidget(self.restore)

        page = QWidget()
        column = QVBoxLayout(page)
        for title, form in (('Inputs', inputs), ('Output', output), ('Filters', filters), ('Canonicals', canonicals),
                            ('The run', the_run)):
            box = QGroupBox(title)
            box.setLayout(form)
            column.addWidget(box)
        column.addLayout(restore)
        column.addStretch()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(page)
        QVBoxLayout(self).addWidget(scroll)

        # each change is checked once the typing stops
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self.refresh)
        for field in (self.rmats_folder, self.gtf_file, self.genome, self.out, self.cache_dir, self.uniprot):
            field.changed.connect(self._timer.start)
        for box in self.splice_types.values():
            box.toggled.connect(self._timer.start)
        for spin in (self.q_lo, self.q_hi, self.processes):
            spin.valueChanged.connect(self._timer.start)
        for box in (self.canonical, self.verbose):
            box.toggled.connect(self._timer.start)
        self.refresh()

    def _row(self, form: QFormLayout, label: str, field, name: str) -> None:
        """ A field, with its problems under it; the command line's help is its status tip, and its label's """
        error = QLabel()
        error.setStyleSheet(ERROR_STYLE)
        error.setWordWrap(True)
        error.hide()
        self.errors[name] = error
        cell = QVBoxLayout()
        cell.setSpacing(2)
        if isinstance(field, QWidget):
            field.setStatusTip(option_tip(name))
            cell.addWidget(field)
        else:
            cell.addLayout(field)
        cell.addWidget(error)
        title = QLabel(label)
        title.setStatusTip(option_tip(name))
        form.addRow(title, cell)

    def arguments(self) -> argparse.Namespace:
        """ The run's arguments, as run_jcast and main.validate take them, with the default cutoff: -r 1 """
        types = [t for t, box in self.splice_types.items() if box.isChecked() and box.isEnabled()]
        uniprot = self.uniprot.text()
        return argparse.Namespace(rmats_folder=self.rmats_folder.text(), gtf_file=self.gtf_file.text(),
                                  genome=self.genome.text(), out=self.out.text(), overwrite=False,
                                  cache_dir=self.cache_dir.text() or None, read=1, model=False,
                                  canonical=self.canonical.isChecked(), uniprot=uniprot or None,
                                  processes=self.processes.value() if uniprot else 1, verbose=self.verbose.isChecked(),
                                  splice_type=set(SPLICE_TYPES) if types == SPLICE_TYPES else types,
                                  qvalue=[self.q_lo.value(), self.q_hi.value()],
                                  )

    def values(self) -> dict:
        """ The page's settings, as DEFAULTS names them: what the main window keeps between sessions """
        return {**{name: getattr(self, name).text() for name in PATHS},
                'splice_types': ','.join(t for t, box in self.splice_types.items() if box.isChecked()),
                'q_lo': self.q_lo.value(), 'q_hi': self.q_hi.value(), 'canonical': self.canonical.isChecked(),
                'processes': self.processes.value(), 'verbose': self.verbose.isChecked()}

    def set_values(self, values: dict) -> None:
        """ Set the page's settings, as values gives them, and check them now """
        for name in PATHS:
            getattr(self, name).setText(values[name])
        chosen = values['splice_types'].split(',')
        for t, box in self.splice_types.items():
            box.setChecked(t in chosen)
        self.q_lo.setValue(values['q_lo'])
        self.q_hi.setValue(values['q_hi'])
        self.canonical.setChecked(values['canonical'])
        self.processes.setValue(values['processes'])      # after --uniprot, which -p needs
        self.verbose.setChecked(values['verbose'])
        self.refresh()

    def _restore_defaults(self) -> None:
        if QMessageBox.question(self, 'Restore the defaults?', 'Clear every path and set every option to its '
                                'default?') != QMessageBox.Yes:
            return
        self._show_required = False
        self.set_values(DEFAULTS)

    def show_required(self) -> None:
        """ Show the fields a run needs that are empty: a run was started without them """
        self._show_required = True
        self.refresh()

    def refresh(self) -> None:
        """ Check the settings now, show the problems next to their fields, and emit changed """
        self._timer.stop()
        self._enable()
        self.problems = self._check()
        notes = {}
        found = main.existing_outputs(self.out.text()) if self.out.text() else []
        if found:
            notes['out'] = 'This folder holds the output of an earlier run: starting the run asks before replacing ' \
                           'it (--overwrite).'
        for name, label in self.errors.items():
            messages = [m for n, m in self.problems if n == name and (m != 'required' or self._show_required)]
            label.setText('\n'.join(messages) or notes.get(name, ''))
            label.setStyleSheet(ERROR_STYLE if messages else '')
            label.setVisible(bool(messages or name in notes))
        self.changed.emit()

    def _enable(self) -> None:
        """ The splice types whose rMATS file is in the folder can be selected; -p needs --uniprot """
        folder = self.rmats_folder.text()
        for t, box in self.splice_types.items():
            present = not os.path.isdir(folder) or os.path.isfile(rmats_path(folder, t))
            box.setEnabled(present)
            box.setStatusTip(option_tip('splice_type') if present else option_tip(
                'splice_type', f'{os.path.basename(rmats_path(folder, t))} is not in the rMATS folder: {t} cannot be '
                               f'selected'))
        self.processes.setEnabled(bool(self.uniprot.text()))
        if not self.uniprot.text():
            self.processes.setValue(1)

    def _check(self) -> list:
        """
        What keeps the run from starting: the fields left empty, no splice type, and the command line's checks. An
        output folder holding JCAST output is not one: the Run page asks before replacing it
        """
        args = self.arguments()
        problems = [(name, 'required') for name in REQUIRED if not getattr(args, name)]
        if not args.splice_type:
            problems.append(('splice_type', 'select at least one splice type'))
        empty = {name for name, _ in problems}
        return problems + [(name, message) for name, message in main.validate(args)
                           if name not in empty and name != 'out']
