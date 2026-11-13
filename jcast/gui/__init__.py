# -*- coding: utf-8 -*-

"""Optional PySide6 GUI for JCAST.

This subpackage is reachable via ``python -m jcast --gui`` or the
``jcast-gui`` console script. The Qt dependency is optional; install with
``pip install 'jcast[gui]'``.

The current build is a skeleton: it opens an empty main window and reserves
the widget structure. See REVIEW.md (§8) for the full design and the v0.4.0
implementation plan.
"""

from jcast.gui.app import launch_gui

__all__ = ["launch_gui"]
