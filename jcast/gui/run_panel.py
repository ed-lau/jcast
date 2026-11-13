# -*- coding: utf-8 -*-

"""Run panel.

Skeleton placeholder. The real implementation:

* runs ``jcast.main.run_jcast`` in a QThread worker so the GUI stays
  responsive,
* tails ``out_dir/logfile.log`` into a QPlainTextEdit,
* shows per-splice-type progress bars,
* increments per-tier counters (T1 / T1 ALT / T2 / T4 / canonical / orphan) by
  parsing the ``fates`` strings written to the log,
* offers "Open output folder" on completion.
"""

from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget


class RunPanel(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Run progress — placeholder."))
        # TODO(v0.4): QThread worker; log tail; tier counters.
