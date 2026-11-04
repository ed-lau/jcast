# -*- coding: utf-8 -*-

""" Shared fixtures. Test data is loaded once per session; caches are kept out of the repo. """

import logging
import shutil
from pathlib import Path

import pytest

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


def _read_gtf_copy(logger, tmp_path_factory, gtf_path):
    """ Read a copy of the GTF so the .cached file is written to a temp dir, not next to the test data. """
    tmp_gtf = tmp_path_factory.mktemp('gtf') / gtf_path.name
    shutil.copy(gtf_path, tmp_gtf)
    return ReadAnnotations(logger=logger, path=tmp_gtf)


@pytest.fixture(scope='session')
def gtf(logger, tmp_path_factory):
    return _read_gtf_copy(logger, tmp_path_factory, HUMAN_GTF)


@pytest.fixture(scope='session')
def rat_gtf(logger, tmp_path_factory):
    return _read_gtf_copy(logger, tmp_path_factory, RAT_GTF)


@pytest.fixture(scope='session')
def genome(logger):
    return ReadGenome(logger=logger, path=HUMAN_GENOME)


@pytest.fixture(scope='session')
def rmats_results(logger):
    return RmatsResults(logger=logger, rmats_dir=RMATS)


@pytest.fixture(autouse=True)
def _run_in_tmp_dir(tmp_path, monkeypatch):
    """ The UniProt cache is written relative to the CWD (REVIEW.md MED-5); keep it out of the repo. """
    monkeypatch.chdir(tmp_path)
