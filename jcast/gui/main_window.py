# -*- coding: utf-8 -*-

"""
The main window: four pages, as tabs, since a rerun with another cutoff means going back: Setup, Model, Run and
Results (reports/2026-11-13-gui-cli-scope.md §4)
"""

from PySide6.QtWidgets import QMainWindow, QMessageBox, QTabWidget

from jcast import __version__
from jcast.gui.model_page import ModelPage
from jcast.gui.results_page import ResultsPage
from jcast.gui.run_page import RunPage
from jcast.gui.setup_page import SetupPage
from jcast.outcomes import SPLICE_TYPES


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"JCAST {__version__}")
        self.resize(980, 760)

        self.setup = SetupPage()
        self.model = ModelPage()
        self.run = RunPage()
        self.results = ResultsPage()
        self.tabs = QTabWidget()
        for page, title in ((self.setup, 'Setup'), (self.model, 'Model'), (self.run, 'Run'),
                            (self.results, 'Results')):
            self.tabs.addTab(page, title)
        self.setCentralWidget(self.tabs)

        self.setup.changed.connect(self._settings_changed)
        self.model.changed.connect(self._settings_changed)
        self.run.setup_needed.connect(self._show_setup)
        self.run.finished.connect(self._run_finished)
        self._settings_changed()

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
