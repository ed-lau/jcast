# -*- coding: utf-8 -*-

"""Main window for the JCAST GUI.

Skeleton: holds a QStackedWidget that future PRs will populate with the
wizard (step 1-2), model view (step 3), and run panel (step 4). The status
bar shows the current step and any in-flight job info.
"""

from PySide6.QtWidgets import QMainWindow, QStackedWidget

from jcast import __version__
from jcast.gui.wizard import WizardPage
from jcast.gui.model_view import ModelView
from jcast.gui.run_panel import RunPanel


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"JCAST {__version__}")
        self.resize(880, 640)

        self._stack = QStackedWidget(self)
        self._stack.addWidget(WizardPage(self))
        self._stack.addWidget(ModelView(self))
        self._stack.addWidget(RunPanel(self))
        self.setCentralWidget(self._stack)

        self.statusBar().showMessage("Ready")
