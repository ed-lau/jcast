# -*- coding: utf-8 -*-

""" Logger. """

import logging
import os

def get_logger(name: str,
               out_path: str,
               verbose: bool = False,
               ) -> logging.Logger:
    """
    Get logger: <out_path>/logfile.log, and warnings and errors on the console. A record logged with
    extra={'jcast_console': False} goes to the log file only. Each junction's steps and sequences are logged at DEBUG
    level, which goes to the log file only if verbose.

    2026-11-12 every junction was logged with its sequences, about 1.7 kB per junction, and each line ended with
    ~module:<module>~function:<module>; junctions.tsv now records each junction (reports/2026-11-12-log-audit.md).

    2026-11-12 handlers from an earlier run in the same process are closed and removed first; they used to stay, so
    a second run logged every line twice, once into the first run's log (REVIEW.md MED-13). close_logger closes this
    run's handlers when it ends.
    """

    file_formatter = logging.Formatter('%(asctime)s~%(levelname)s~%(message)s')
    console_formatter = logging.Formatter('%(levelname)s -- %(message)s')

    file_handler = logging.FileHandler(os.path.join(out_path, "logfile.log"), mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(console_formatter)
    console_handler.addFilter(lambda record: getattr(record, 'jcast_console', True))

    logger = logging.getLogger(name)
    close_logger(logger)
    logger.propagate = False
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def close_logger(logger: logging.Logger) -> None:
    """ Close and remove the logger's handlers """
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
