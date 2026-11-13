# -*- coding: utf-8 -*-

"""Mixture-model preview widget.

Skeleton placeholder. The real implementation embeds
``jcast.model.plot_general_mixture_model`` (or its scikit-learn replacement
in v0.4.0) into a ``FigureCanvasQTAgg``, with the decision boundary drawn as
a draggable vertical line so the user can override the modeled cutoff before
launching the run.
"""

from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget


class ModelView(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Model preview — placeholder."))
        # TODO(v0.4): embed FigureCanvasQTAgg; draggable decision boundary.
