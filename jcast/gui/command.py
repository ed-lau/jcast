# -*- coding: utf-8 -*-

"""
The GUI's run as a command line: the arguments of `jcast run` for the settings of the Setup and Model pages. The GUI
shows them, and runs them as `python -P -m jcast run … --progress` in a child process (process;
reports/2026-11-13-gui-cli-scope.md §3). No Qt here, so that the tests of the command line need no PySide6.
"""

import functools
import os
import shlex
import site
import sys

from jcast.outcomes import SPLICE_TYPES


def arguments(args) -> list:
    """
    The arguments of `jcast run` for a run's arguments, leaving out those at their default: the command line parses
    them into the same arguments again

    :param args: the run's arguments, as run_jcast and main.validate take them (argparse.Namespace)
    :return: the words after `jcast run`
    """
    words = [args.rmats_folder, args.gtf_file, args.genome, '-o', args.out]
    if args.overwrite:
        words.append('--overwrite')
    if args.cache_dir:
        words += ['--cache_dir', args.cache_dir]
    if args.model:
        words.append('-m')
    if args.read != 1:
        words += ['-r', str(args.read)]
    if list(args.qvalue) != [0, 1]:
        words += ['-q'] + [_number(q) for q in args.qvalue]
    if set(args.splice_type) != set(SPLICE_TYPES):
        words += ['-s'] + [t for t in SPLICE_TYPES if t in args.splice_type]
    if args.canonical:
        words.append('-c')
    if args.uniprot:
        words += ['--uniprot', args.uniprot]
    if args.processes != 1:
        words += ['-p', str(args.processes)]
    if args.verbose:
        words.append('-v')
    return words


def command_line(args) -> str:
    """ The command the GUI shows: jcast run and the arguments, quoted for a POSIX shell """
    return shlex.join(['jcast', 'run'] + arguments(args))


def process(args) -> tuple:
    """
    The child process of a run: `python -P -m jcast run <arguments> --progress`, with this Python. The child runs the
    jcast the GUI runs, whatever its working directory: -P keeps that off sys.path (python -m puts it first), and this
    jcast's folder is put first on PYTHONPATH unless it is a site-packages folder (a source checkout; an editable
    install)

    :return: the program, its arguments, and the PYTHONPATH
    """
    import jcast

    root = os.path.dirname(os.path.dirname(os.path.abspath(jcast.__file__)))
    paths = [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]
    if root not in site.getsitepackages() + [site.getusersitepackages()] and root not in paths:
        paths.insert(0, root)
    return sys.executable, ['-P', '-m', 'jcast', 'run'] + arguments(args) + ['--progress'], os.pathsep.join(paths)


def _number(q: float) -> str:
    """ 0.05, not 0.05000000000000000277; repr if the short form would not read back as the same number """
    short = f'{q:g}'
    return short if float(short) == q else repr(float(q))


@functools.lru_cache
def option_help(name: str) -> str:
    """ The help of an argument or option of `jcast run`, as jcast run -h shows it: the GUI's tooltips """
    import typer.main
    from jcast import main

    for param in typer.main.get_command(main.app).commands['run'].params:
        if param.name == name:
            return param.help
    raise KeyError(name)
