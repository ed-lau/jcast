# -*- coding: utf-8 -*-

""" Canonicals under -c and UniProt API failures (PR 8; reports/2026-11-12-thresholds-runs-canonical-audit.md §4) """

import argparse
import logging

import pytest
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import fates, main, params
from jcast.sequences import Sequence, UniProtAPIError
from jcast.uniprot import UniProtFasta

from tests.test_junction import _junction
from tests.test_uniprot import UNIPROT, _no_api

API_LOOKUP = Sequence.get_canonical_aa_uniprot


def _args(out, **extra):
    a = dict(out=str(out), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
    a.update(extra)
    return argparse.Namespace(**a)


def _canonicals(out):
    path = out / 'psq_canonical.fasta'
    return [str(r.seq) for r in SeqIO.parse(path, 'fasta')] if path.exists() else []


@pytest.fixture
def gtf_only(monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)


@pytest.fixture(scope='module')
def uniprot():
    return UniProtFasta(logger=logging.getLogger('jcast-test'), path=UNIPROT)


# ---- -c on every path (MED-14) ----

def test_canonical_written_for_a_t1_junction(gtf, genome, rmats_results, tmp_path, gtf_only):
    """ TSPAN3's only junction is T1 and neither slice is its canonical; -c used to write nothing for it """
    junction = _junction(rmats_results.rmats_mxe, 441)

    fate = main._translate_one(junction, gtf, genome, _args(tmp_path), pred_bound=1)

    sequence = Sequence(_junction(rmats_results.rmats_mxe, 441))
    sequence.get_canonical_aa(gtf, genome.genome)
    assert fate.startswith('SUCCESS TIER 1.')
    assert _canonicals(tmp_path) == [str(sequence.canonical_aa.seq)]


def test_no_canonical_without_c(gtf, genome, rmats_results, tmp_path, gtf_only):
    main._translate_one(_junction(rmats_results.rmats_mxe, 441), gtf, genome, _args(tmp_path, canonical=False),
                        pred_bound=1)
    assert _canonicals(tmp_path) == []


@pytest.mark.parametrize('use_uniprot', [False, True])
def test_canonical_written_for_a_masked_junction(gtf, genome, rmats_results, tmp_path, monkeypatch, uniprot,
                                                 use_uniprot):
    """ the masked-base check came before the canonical lookup, so -c wrote nothing for a masked junction """
    monkeypatch.setattr(params, 'use_gtf_only', not use_uniprot)
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    junction = _junction(rmats_results.rmats_se, 6076)        # PIAS1's only chr15 junction
    seq = genome.genome['15']
    masked = argparse.Namespace(genome={'15': seq[:junction.alt1_es - 1] + 'N' * (junction.alt1_ee - junction.alt1_es + 1)
                                        + seq[junction.alt1_ee:]})

    fate = main._translate_one(junction, gtf, masked, _args(tmp_path), pred_bound=1,
                               uniprot=uniprot if use_uniprot else None)

    assert fate == fates.SKIPPED_MASKED
    assert len(_canonicals(tmp_path)) == 1


def test_failing_junction_writes_no_canonical(gtf, genome, rmats_results, tmp_path, gtf_only, monkeypatch):
    """ a junction that raises writes nothing, its canonical included (PR 7) """
    monkeypatch.setattr(Sequence, 'translate', lambda self, log=True: 1 / 0)
    fate, records = main._translate_safely(_junction(rmats_results.rmats_mxe, 441), gtf, genome, _args(tmp_path),
                                           pred_bound=1, uniprot=None)
    assert fate.startswith('FAILURE (exception). ZeroDivisionError') and records == []


# ---- write_canonical's API fallback (MED-9) ----

def test_empty_canonical_makes_no_api_request(rmats_results, tmp_path, monkeypatch):
    """ an empty canonical is not looked up again; GTF-only runs make no API request """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    sequence = Sequence(_junction(rmats_results.rmats_se, 6426))

    sequence.write_canonical(out_dir=str(tmp_path))

    assert _canonicals(tmp_path) == []


# ---- UniProt API failures (MED-10) ----

def _response(status, text=''):
    r = requests.models.Response()
    r.status_code, r._content, r.encoding, r.url = status, text.encode(), 'utf-8', 'stub'
    return r


@pytest.fixture
def api(monkeypatch, tmp_path):
    """ the real API lookup, with Session.get replaced; returns the list of calls (url, kwargs) """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', API_LOOKUP)
    from jcast import cache
    monkeypatch.setattr(cache, '_directory', str(tmp_path / 'api_cache'))
    calls = []

    def install(answer):
        def get(self, url, **kwargs):
            calls.append((url, kwargs, self.adapters['https://'].max_retries))
            if isinstance(answer, Exception):
                raise answer
            return answer
        monkeypatch.setattr(requests.Session, 'get', get)
    return calls, install


def test_api_entry_and_no_entry(rmats_results, api):
    calls, install = api
    install(_response(200, '>sp|P1|A_HUMAN A\nMAAA\n'))
    assert str(Sequence(_junction(rmats_results.rmats_se, 6426)).get_canonical_aa_uniprot().seq) == 'MAAA'
    install(_response(200, ''))
    assert len(Sequence(_junction(rmats_results.rmats_se, 6076)).get_canonical_aa_uniprot()) == 0
    url, kwargs, retries = calls[-1]
    assert kwargs['timeout'] == params.uniprot_timeout
    assert 429 in retries.status_forcelist and retries.total == params.uniprot_max_retries


@pytest.mark.parametrize('answer', [_response(429), _response(403), _response(500),
                                    requests.exceptions.ConnectionError('refused'),
                                    requests.exceptions.ReadTimeout('timed out')])
def test_api_failure_raises(rmats_results, api, answer):
    """ a refused or failed request is an error, never "no canonical" (it used to orphan every record) """
    calls, install = api
    install(answer)
    with pytest.raises(UniProtAPIError, match='ENSG'):
        Sequence(_junction(rmats_results.rmats_se, 6426)).get_canonical_aa_uniprot()


def test_api_failure_fails_the_junction(gtf, genome, rmats_results, tmp_path, api, monkeypatch):
    """ in a run, the junction gets the exception fate and writes nothing, not orphans """
    monkeypatch.setattr(params, 'use_gtf_only', False)
    calls, install = api
    install(_response(429))
    fate, records = main._translate_safely(_junction(rmats_results.rmats_mxe, 441), gtf, genome, _args(tmp_path),
                                           pred_bound=1, uniprot=None)
    assert fate.startswith('FAILURE (exception). UniProtAPIError') and records == []


# ---- --uniprot: gene names only when unique in the GTF (mapping audit §7) ----

def test_gene_name_used_only_when_unique(uniprot):
    record = SeqRecord(Seq('MGTF'), id='gtf')

    chosen = uniprot._choose('G1', 'TSPAN3', lambda: set(), lambda: record, unique_name=True)
    assert chosen[1] == 'sp by gene name'

    chosen = uniprot._choose('G2', 'TSPAN3', lambda: set(), lambda: record, unique_name=False)
    assert chosen[1:3] == ('gtf', '')


def test_gene_name_counts(gtf):
    from jcast.annots import gene_name_is_unique
    assert gene_name_is_unique(gtf, 'TSPAN3')
    assert gene_name_is_unique(gtf, 'NOT_A_GENE_NAME')
    names = gtf.annot.loc[gtf.annot['feature'] == 'gene', 'gene_name'].value_counts()
    duplicated = names[names > 1].index
    assert all(not gene_name_is_unique(gtf, n) for n in duplicated)
