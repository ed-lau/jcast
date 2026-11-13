# -*- coding: utf-8 -*-

"""Input wizard.

Skeleton placeholder. The real implementation collects:

  Step 1 — paths
    * rMATS folder
    * GTF file
    * Genome FASTA (optionally .gz)
    * Output directory
    * Splice-type checkboxes (MXE / SE / RI / A5SS / A3SS)

  Step 2 — filters
    * q-value low/high
    * read-cutoff source: manual integer vs. mixture-model (-m)
    * --canonical toggle

and produces a RunConfig that the run panel hands to ``jcast.main.run_jcast``
in a QThread worker.
"""

from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget


class WizardPage(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Inputs & filters — placeholder."))
        # TODO(v0.4): replace with the wizard described in REVIEW.md §8.
