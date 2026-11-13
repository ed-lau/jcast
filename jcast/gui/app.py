# -*- coding: utf-8 -*-

"""GUI bootstrap. Creates the QApplication and shows the MainWindow."""

import sys


_MISSING_PYSIDE_MSG = (
    "JCAST's graphical interface requires PySide6, which is an optional "
    "dependency. Install it with:\n\n"
    "    pip install 'jcast[gui]'\n"
)


def launch_gui() -> int:
    """Launch the JCAST GUI.

    Returns the Qt exit code. Exits with code 1 (and prints an install hint)
    if PySide6 is not importable.
    """
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print(_MISSING_PYSIDE_MSG, file=sys.stderr)
        return 1

    from jcast.gui.main_window import MainWindow

    app = QApplication.instance() or QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(launch_gui())
