# -*- coding: utf-8 -*-

"""The optional PySide6 GUI for JCAST.

Opened by ``jcast gui``, ``python -m jcast --gui`` or the ``jcast-gui``
console script. The Qt dependency is optional; install with
``pip install 'jcast[gui]'``.

Four pages (reports/2026-11-13-gui-cli-scope.md §4): Setup, the run's
settings; Model, the -m model before the run, where the read count cutoff
is set; Run, which runs ``python -m jcast run …`` in a child process; and
Results, the summary of any run's output folder.
"""

from jcast.gui.app import launch_gui

__all__ = ["launch_gui"]
