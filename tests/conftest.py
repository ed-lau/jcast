# -*- coding: utf-8 -*-

""" Shared fixtures. Test data is loaded once per session; caches are kept in a temporary directory. """

import logging
from pathlib import Path

import pytest

from jcast import cache
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.junctions import RmatsResults

DATA = Path(__file__).parent / 'data'
HUMAN_GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
RAT_GTF = DATA / 'genome' / 'Rattus_norvegicus.Rnor_6.0.100.chromosome.12.gtf'
HUMAN_GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
RMATS = DATA / 'rmats'


@pytest.fixture(scope='session')
def logger():
    return logging.getLogger('jcast-test')


@pytest.fixture(scope='session')
def cache_dir(tmp_path_factory):
    """ The GTF and UniProt API caches, shared by the session's tests, never the user's cache directory """
    path = tmp_path_factory.mktemp('cache')
    cache.set_directory(path)
    return path


@pytest.fixture(scope='session')
def gtf(logger, cache_dir):
    return ReadAnnotations(logger=logger, path=HUMAN_GTF)


@pytest.fixture(scope='session')
def rat_gtf(logger, cache_dir):
    return ReadAnnotations(logger=logger, path=RAT_GTF)


@pytest.fixture(scope='session')
def genome(logger):
    return ReadGenome(logger=logger, path=HUMAN_GENOME)


@pytest.fixture(scope='session')
def rmats_results(logger):
    return RmatsResults(logger=logger, rmats_dir=RMATS)


@pytest.fixture(autouse=True)
def _wide_console(monkeypatch):
    """ Typer prints a usage error in a box as wide as the console, 80 columns when it is not a terminal; a wide one
    keeps each message on one line. And it takes the console for a terminal, with colours and styles in the text, when
    GITHUB_ACTIONS (on CI), FORCE_COLOR or PY_COLORS is set; the tests' console is not one """
    monkeypatch.setenv('COLUMNS', '1000')
    monkeypatch.setattr('typer.rich_utils.FORCE_TERMINAL', False)


@pytest.fixture(autouse=True)
def _run_in_tmp_dir(tmp_path, monkeypatch, cache_dir):
    """ Run each test in a temp dir, with the caches in the session's cache directory (run_jcast resets it from
    args.cache_dir, so it is set again for every test) """
    monkeypatch.chdir(tmp_path)
    cache.set_directory(cache_dir)
