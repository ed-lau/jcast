# -*- coding: utf-8 -*-

""" A function run in a background thread, with its result delivered in the GUI thread """

import threading

from PySide6.QtCore import QObject, Signal, Slot


class Background(QObject):
    """
    Runs fn() in a thread and calls done(result), or failed(exception), in the GUI thread: the model preview takes up
    to 20 s on a whole genome (reports/2026-11-13-gui-cli-scope.md §6). A daemon thread, so that quitting does not wait
    for it. The caller drops a result it no longer wants (a newer task replaced this one)
    """

    _ended = Signal(object, object)

    def __init__(self, parent: QObject, fn, done, failed):
        super().__init__(parent)
        self._done, self._failed = done, failed
        self.running = True
        # self lives in the GUI thread, so _ended, emitted in the background thread, is delivered there
        self._ended.connect(self._end)
        threading.Thread(target=self._work, args=(fn,), daemon=True).start()

    def _work(self, fn) -> None:
        try:
            result, error = fn(), None
        except Exception as e:
            result, error = None, e
        try:
            self._ended.emit(result, error)
        except RuntimeError:    # the window closed meanwhile, and this object with it
            pass

    @Slot(object, object)
    def _end(self, result, error) -> None:
        self.running = False
        if error is None:
            self._done(result)
        else:
            self._failed(error)
        self.deleteLater()
