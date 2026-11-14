# -*- coding: utf-8 -*-

"""
The main window: four pages, as tabs, since a rerun with another cutoff means going back: Setup, Model, Run and
Results (reports/2026-11-13-gui-cli-scope.md §4); and under them, the help of what the mouse is over
"""

from PySide6.QtCore import QEvent, Qt
from PySide6.QtWidgets import QFrame, QLabel, QMainWindow, QMessageBox, QTabWidget, QVBoxLayout, QWidget

from jcast import __version__
from jcast.gui.app import icon
from jcast.gui.model_page import ModelPage
from jcast.gui.results_page import ResultsPage
from jcast.gui.run_page import RunPage
from jcast.gui.setup_page import SetupPage
from jcast.outcomes import SPLICE_TYPES


HELP = 'Point at a setting, a button or a table to see what it does.'


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"JCAST {__version__}")
        self.setWindowIcon(icon())
        self.resize(980, 800)

        self.setup = SetupPage()
        self.model = ModelPage()
        self.run = RunPage()
        self.results = ResultsPage()
        self.tabs = QTabWidget()
        for page, title in ((self.setup, 'Setup'), (self.model, 'Model'), (self.run, 'Run'),
                            (self.results, 'Results')):
            self.tabs.addTab(page, title)

        # the help: the status tip of what the mouse is over (the command line's help for a setting), in a panel of
        # its own rather than the status bar, which shows one line
        self.help = QLabel(HELP)
        self.help.setTextFormat(Qt.RichText)
        self.help.setWordWrap(True)
        self.help.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.help.setFrameShape(QFrame.StyledPanel)
        self.help.setMargin(6)
        self.help.setMinimumHeight(4 * self.help.fontMetrics().lineSpacing())
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.addWidget(self.tabs, 1)
        layout.addWidget(self.help)
        self.setCentralWidget(central)

        self.setup.changed.connect(self._settings_changed)
        self.model.changed.connect(self._settings_changed)
        self.run.setup_needed.connect(self._show_setup)
        self.run.finished.connect(self._run_finished)
        self._settings_changed()

    def event(self, event) -> bool:
        """ A status tip, sent when the mouse enters or leaves a widget that has one: shown in the help panel """
        if event.type() == QEvent.StatusTip:
            self.help.setText(event.tip() or HELP)
            return True
        return super().event(event)

    def arguments(self):
        """ The run's arguments: the Setup page's, with the Model page's cutoff """
        args = self.setup.arguments()
        args.model, args.read = self.model.options()
        return args

    def _settings_changed(self) -> None:
        args, problems = self.arguments(), self.setup.problems
        folder = args.rmats_folder if args.rmats_folder and 'rmats_folder' not in dict(problems) else None
        self.model.set_settings(folder, [t for t in SPLICE_TYPES if t in args.splice_type], tuple(args.qvalue))
        self.run.set_arguments(args, problems)

    def _show_setup(self) -> None:
        self.setup.show_required()
        self.tabs.setCurrentWidget(self.setup)
        self.statusBar().showMessage('Complete the settings marked on the Setup page to start a run.', 10000)

    def _run_finished(self, out: str, code: int) -> None:
        self.results.load(out)
        if code == 0:
            self.tabs.setCurrentWidget(self.results)

    def closeEvent(self, event) -> None:
        """ Quitting during a run stops it, as Stop does """
        if self.run.running():
            if QMessageBox.question(self, 'Stop the run and quit?',
                                    'A run is in progress. Its output will be incomplete.') != QMessageBox.Yes:
                event.ignore()
                return
            self.run.stop_run(ask=False)
            self.run.wait()
        event.accept()
