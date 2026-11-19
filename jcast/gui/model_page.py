# -*- coding: utf-8 -*-

"""
The Model page: the -m model of the rMATS folder before the run (model.preview), and the read count cutoff: the
model's (-m), or another (-r), set in a number field or by clicking on the plot. A table shows how many junctions of
each splice type the run would keep (reports/2026-11-13-gui-cli-scope.md §4). Two more charts, beside the model, show
what a cutoff keeps (PR 27; reports/2026-11-15-read-count-model-audit.md §16.2): the junctions a run would translate at
each cutoff, annotated and novel, and by read count the events whose alternative region keeps the frame
"""

import math

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QButtonGroup, QHBoxLayout, QHeaderView, QLabel, QProgressBar, QRadioButton, QSpinBox,
                               QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from jcast import model
from jcast.gui.command import option_tip
from jcast.gui.tasks import Background
from jcast.outcomes import SPLICE_TYPES

HEADERS = ['Splice type', 'Junctions', 'Read count below the cutoff', 'FDR outside -q', 'Kept']


class ModelPage(QWidget):
    """ changed is emitted when the cutoff's option (options) changes """

    changed = Signal()

    def __init__(self, parent: QWidget = None):
        super().__init__(parent)
        self._folder = None             # the Setup page's rMATS folder; None until it is one
        self._splice_types, self._qvalue = list(SPLICE_TYPES), (0, 1)
        self._task = None               # the fit running, and its folder
        self._fit = None                # the fit shown: its folder, then what model.preview returns
        self._cutoff = None             # the cutoff's line and label on the plot
        self._chart_cutoffs = {}        # chart -> the cutoff's line and label on it
        self._charts_for = None         # the splice types and -q the charts were drawn for

        self.status = QLabel('Choose an rMATS folder on the Setup page.')
        self.status.setWordWrap(True)
        self.busy = QProgressBar()
        self.busy.setRange(0, 0)
        self.busy.setMaximumWidth(160)
        self.busy.hide()

        self.figure = Figure(figsize=(8, 3.2), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(260)
        self.canvas.setStatusTip('Click or drag on the plot to set a read count cutoff of your own (-r). The plot: the '
                                 'read counts of every junction in the rMATS folder, the -m model fitted to them '
                                 '(trace and expressed components, and their sum), and the cutoff.')
        self.canvas.mpl_connect('button_press_event', self._clicked)
        self.canvas.mpl_connect('motion_notify_event', self._clicked)

        self.kept_figure = Figure(figsize=(8, 3.2), constrained_layout=True)
        self.kept_canvas = FigureCanvasQTAgg(self.kept_figure)
        self.kept_canvas.setMinimumHeight(260)
        self.kept_canvas.setStatusTip(
            'Click or drag on the plot to set a read count cutoff of your own (-r). The plot: how many junctions of '
            'the splice types selected on the Setup page a run would translate at each cutoff, with -q; for the events '
            'rMATS calls annotated and those it calls novel (its fromGTF.novelJunction and fromGTF.novelSpliceSite '
            'files), if they are in the folder.')
        self.kept_canvas.mpl_connect('button_press_event', self._clicked_kept)
        self.kept_canvas.mpl_connect('motion_notify_event', self._clicked_kept)
        self.frame_figure = Figure(figsize=(8, 3.2), constrained_layout=True)
        self.frame_canvas = FigureCanvasQTAgg(self.frame_figure)
        self.frame_canvas.setMinimumHeight(260)
        self.frame_canvas.setStatusTip(
            'By read count, the share of events whose alternative region (the skipped exon or the retained intron; for '
            'the other splice types, the difference between the two forms) is a multiple of 3 nt, so that both forms '
            'keep the reading frame. A random length does in one event of three. Events of the splice types selected '
            'on the Setup page; annotated and novel, as rMATS calls them.')
        self.charts = QTabWidget()
        for canvas, title in ((self.canvas, 'Read counts and the model'),
                              (self.kept_canvas, 'Junctions translated by cutoff'),
                              (self.frame_canvas, 'Frame by read count')):
            self.charts.addTab(canvas, title)
        self._charts = {'kept': (self.kept_figure, self.kept_canvas), 'frame': (self.frame_figure, self.frame_canvas)}

        self.use_read = QRadioButton('Read count cutoff (-r)')
        self.use_read.setStatusTip(option_tip('read'))
        self.read = QSpinBox()
        self.read.setRange(0, 10 ** 9)
        self.read.setValue(1)
        self.read.setMinimumWidth(110)
        self.read.setStatusTip(option_tip('read'))
        self.use_model = QRadioButton("The model's cutoff (-m)")
        self.use_model.setStatusTip(option_tip('model'))
        self.use_model.setEnabled(False)
        options = QButtonGroup(self)
        options.addButton(self.use_read)
        options.addButton(self.use_model)
        self.use_read.setChecked(True)
        self.use_read.toggled.connect(self._option_changed)
        self.read.valueChanged.connect(self._read_changed)

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(self.table.fontMetrics().height() + 8)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setStatusTip('The junctions of each splice type selected on the Setup page: those whose read count '
                                'is below the cutoff, those whose FDR is outside -q, and those kept, which the run '
                                'translates.')
        self.table_note = QLabel('Junctions of the splice types selected on the Setup page. Those kept are translated; '
                                 'the others are skipped (skipped_low_count and skipped_fdr in junctions.tsv).')
        self.table_note.setWordWrap(True)

        top = QHBoxLayout()
        top.addWidget(self.status, 1)
        top.addWidget(self.busy)
        cutoff = QHBoxLayout()
        cutoff.addWidget(self.use_read)
        cutoff.addWidget(self.read)
        cutoff.addWidget(QLabel('reads'))
        cutoff.addSpacing(24)
        cutoff.addWidget(self.use_model)
        cutoff.addStretch()
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.charts, 3)
        layout.addLayout(cutoff)
        layout.addWidget(QLabel('Click or drag on the model, or on the junctions translated, to set a cutoff of your '
                                'own.'))
        layout.addWidget(self.table, 2)
        layout.addWidget(self.table_note)

        # the table, once the cutoff stops moving: on a whole genome it takes a moment
        self._table_timer = QTimer(self)
        self._table_timer.setSingleShot(True)
        self._table_timer.setInterval(150)
        self._table_timer.timeout.connect(self.update_table)

    # ---- the settings, and the fit ----

    def set_settings(self, folder, splice_types: list, qvalue: tuple) -> None:
        """ The Setup page's rMATS folder (None if it is not one), splice types and -q """
        if (splice_types, qvalue) != (self._splice_types, self._qvalue):
            self._splice_types, self._qvalue = splice_types, qvalue
            self._table_timer.start()
        if folder != self._folder:
            self._folder = folder
            self._clear()
            if self.isVisible():
                self.fit()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self.fit()

    def fit(self) -> None:
        """ Fit the model of the folder in the background, unless it is shown or being fitted """
        folder = self._folder
        if folder is None or (self._fit is not None and self._fit[0] == folder) or \
                (self._task is not None and self._task[1] == folder):
            return
        self.status.setText(f'Fitting the -m model of every rMATS file in {folder}…')
        self.busy.show()
        task = Background(self, lambda: model.preview(folder), lambda r: self._fitted(folder, r),
                          lambda e: self._failed(folder, e))
        self._task = (task, folder)

    def _fitted(self, folder, result: tuple) -> None:
        if folder != self._folder:      # the folder changed meanwhile
            return
        self._task = None
        self.busy.hide()
        self._fit = (folder, *result)
        counts, ln_sjc, params, cutoff, no_fit = result
        files = sorted(set(counts.junction_type), key=SPLICE_TYPES.index)
        read = f'{len(counts):,} junctions of every rMATS file in the folder ({", ".join(files) or "none"})'
        if cutoff is not None:
            self.status.setText(f'The -m model of {read}: its cutoff is {cutoff} reads.')
            self.use_model.setText(f"The model's cutoff (-m): {cutoff} reads")
            self.use_model.setEnabled(True)
        else:
            self.status.setText(f'No -m model of {read}: the fit {no_fit}. A run with -m would fall back to -r.')
        self._draw()
        self.update_table()

    def _failed(self, folder, error: Exception) -> None:
        if folder != self._folder:
            return
        self._task = None
        self.busy.hide()
        self.status.setText(f'The model could not be read from {folder}: {error}')

    def _clear(self) -> None:
        """ Another folder: the fit shown is not its model, and a fit still running is dropped when it ends """
        self._fit, self._cutoff, self._task = None, None, None
        self.busy.hide()
        if self.use_model.isChecked():
            self.use_read.setChecked(True)
        self.use_model.setText("The model's cutoff (-m)")
        self.use_model.setEnabled(False)
        self.status.setText('Choose an rMATS folder on the Setup page.' if self._folder is None else '')
        self._draw()
        self.update_table()

    # ---- the cutoff ----

    def options(self) -> tuple:
        """ -m and -r for the run: (True, 1) with the model's cutoff, (False, the cutoff) with another """
        return (True, 1) if self.use_model.isChecked() else (False, self.read.value())

    def cutoff(self) -> int:
        """ The cutoff shown: the model's with -m """
        return self._fit[4] if self.use_model.isChecked() else self.read.value()

    def _option_changed(self) -> None:
        self._move()
        self.changed.emit()

    def _read_changed(self) -> None:
        if self.use_read.isChecked():
            self._move()
            self.changed.emit()
        else:
            self.use_read.setChecked(True)

    def _clicked(self, event) -> None:
        """ A click or a drag on the plot sets the cutoff there (-r) """
        if self._cutoff is None or event.xdata is None or (event.name == 'motion_notify_event' and not event.button):
            return
        self.read.setValue(max(0, round(math.exp(event.xdata) - 1)))
        self.use_read.setChecked(True)

    def _clicked_kept(self, event) -> None:
        """ A click or a drag on the junctions translated sets the cutoff there (-r) """
        if 'kept' not in self._chart_cutoffs or event.xdata is None or \
                (event.name == 'motion_notify_event' and not event.button):
            return
        self.read.setValue(max(1, round(event.xdata)))
        self.use_read.setChecked(True)

    def _move(self) -> None:
        if self._cutoff is not None:
            model.move_cutoff(*self._cutoff, self.cutoff())
            self.canvas.draw_idle()
        for chart, (line, label) in self._chart_cutoffs.items():
            model.move_chart_cutoff(line, label, self.cutoff())     # its numbers come with the table
            self._charts[chart][1].draw_idle()
        self._table_timer.start()

    # ---- the plot and the table ----

    def _draw(self) -> None:
        self.figure.clear()
        self._cutoff = None
        if self._fit is not None:
            _, counts, ln_sjc, params, _, _ = self._fit
            if len(ln_sjc):
                self._cutoff = model.draw_model(self.figure.subplots(), ln_sjc, params, self.cutoff())
            else:
                self.figure.text(0.5, 0.5, 'The rMATS files hold no junctions.', ha='center')
        self.canvas.draw_idle()
        self._draw_charts()

    def _draw_charts(self) -> None:
        """ The junctions translated at each cutoff, and the frame by read count, for the splice types and -q """
        self._chart_cutoffs = {}
        self._charts_for = (list(self._splice_types), tuple(self._qvalue))
        for figure, _ in self._charts.values():
            figure.clear()
        if self._fit is not None:
            counts = self._fit[1]
            if len(counts):
                kept = model.kept_at_cutoffs(counts, self._qvalue, self._splice_types)
                self._chart_cutoffs['kept'] = model.draw_kept(self.kept_figure.subplots(), kept, self.cutoff())
                frame = model.frame_by_count(counts, self._splice_types)
                self._chart_cutoffs['frame'] = model.draw_frame(self.frame_figure.subplots(), frame, self.cutoff())
            else:
                for figure, _ in self._charts.values():
                    figure.text(0.5, 0.5, 'The rMATS files hold no junctions.', ha='center')
        for _, canvas in self._charts.values():
            canvas.draw_idle()

    def update_table(self) -> None:
        """ The junctions of each selected splice type, and how many the cutoff and -q keep (model.filter_counts) """
        self._table_timer.stop()
        rows = []
        if self._fit is not None:
            counts = self._fit[1]
            table = model.filter_counts(counts, self.cutoff(), self._qvalue)
            table = table.loc[[t for t in table.index if t in self._splice_types]]
            rows = [[t] + list(r) for t, r in zip(table.index, table.itertuples(index=False))]
            if len(rows) > 1:
                rows.append(['total'] + list(table.sum()))
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                item = QTableWidgetItem(value if j == 0 else f'{int(value):,}')
                item.setTextAlignment((Qt.AlignLeft if j == 0 else Qt.AlignRight) | Qt.AlignVCenter)
                self.table.setItem(i, j, item)

        # the charts, for other splice types or -q; and the junctions translated at the cutoff, on its label
        if self._charts_for != (list(self._splice_types), tuple(self._qvalue)):
            self._draw_charts()
        if 'kept' in self._chart_cutoffs:
            at = model.kept_at_cutoffs(self._fit[1], self._qvalue, self._splice_types, cutoffs=[self.cutoff()])
            text = ', '.join(f'{n:,} {name}' if name != 'all' else f'{n:,} junctions' for name, n in at.iloc[0].items())
            model.move_chart_cutoff(*self._chart_cutoffs['kept'], self.cutoff(), text)
            self.kept_canvas.draw_idle()
