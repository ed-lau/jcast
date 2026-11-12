# -*- coding: utf-8 -*-

""" Where JCAST keeps its caches: the GTF read into a table, and the UniProt API's answers.

2026-11-12 the caches were written next to the GTF and into the working directory; they are now in the user's cache
directory, or --cache_dir. See reports/2026-11-12-thresholds-runs-canonical-audit.md §3.5
"""

import logging
import os

import platformdirs

# --cache_dir; None for the user's cache directory (~/Library/Caches/jcast, ~/.cache/jcast, ...)
_directory = None

# directories that could not be created, warned about once
_unusable = set()


def set_directory(path) -> None:
    """ Use this directory for the caches; None for the user's cache directory """
    global _directory
    _directory = None if path is None else os.fspath(path)


def directory(sub: str = '') -> str:
    """
    The cache directory (or a subdirectory of it), created if needed.

    :param sub: subdirectory
    :return: its path, or None if it cannot be created; the caller then runs without the cache
    """
    path = os.path.join(_directory or platformdirs.user_cache_dir('jcast'), sub)
    if path in _unusable:
        return None
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        _unusable.add(path)
        logging.getLogger('jcast').warning(f'Cannot create the cache directory {path} ({e}); running without it.')
        return None
    return path
