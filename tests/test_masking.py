# -*- coding: utf-8 -*-

""" Soft-masked and ambiguous genome bases (PR 3; reports/2026-11-11-masked-bases-audit.md) """

import argparse
import types

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from jcast import fates, main, params
from jcast import helpers as h

from tests.test_junction import _junction


@pytest.mark.parametrize('nt, strand, terminate, expected', [
    ('ATGAAATAG', '+', True, 'MK'),
    ('atgaaatag', '+', True, 'MK'),
    ('ctatttcat', '-', True, 'MK'),
    ('ATGNNNTAG', '+', True, ''),
    ('ATGRAATAG', '+', False, ''),
    ('CTARTTCAT', '-', True, ''),
])
def test_make_pep_never_raises(nt, strand, terminate, expected):
    """ lowercase is translated as uppercase; a codon with an ambiguous base makes the sequence untranslatable """

    assert h.make_pep(Seq(nt), strand, 0, terminate=terminate) == expected


def test_get_complementary_iupac():
    """ the fallback for non-ACGT bases was unreachable (REVIEW.md HIGH-3) """

    assert h.get_complementary('ATGCAA') == 'TTGCAT'
    assert h.get_complementary('acgRN') == 'NYcgt'


def _run(junction, gtf, genome_index, out_dir, monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    out_dir.mkdir()
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
    genome = types.SimpleNamespace(genome=genome_index)
    return main._translate_one(junction, gtf, genome, args, pred_bound=args.read)


def _fastas(out_dir):
    return {p.name: [(r.id, str(r.seq)) for r in SeqIO.parse(p, 'fasta')] for p in sorted(out_dir.glob('*.fasta'))}


@pytest.fixture(scope='module')
def softmasked(genome):
    """ the chr15 test genome entirely in lowercase """
    return {k: v.lower() for k, v in genome.genome.items()}


@pytest.mark.parametrize('jxn_type, rmats_id', [
    ('mxe', 430),   # first junction to crash on a soft-masked genome (audit §4.1)
    ('se', 5739),   # MFGE8, -: T1
    ('se', 28030),  # FAM96A, -: T1 ALT TERMINUS
    ('se', 6592),   # EMC4, +: T1 ALT TERMINUS
])
def test_softmasked_genome_gives_unmasked_output(jxn_type, rmats_id, gtf, genome, softmasked, rmats_results,
                                                 tmp_path, monkeypatch):
    """ a soft-masked genome raised KeyError; it must give exactly the unmasked fates and sequences """

    rma = getattr(rmats_results, f'rmats_{jxn_type}')
    fate = _run(_junction(rma, rmats_id), gtf, genome.genome, tmp_path / 'unmasked', monkeypatch)
    fate_sm = _run(_junction(rma, rmats_id), gtf, softmasked, tmp_path / 'softmasked', monkeypatch)

    assert fate_sm == fate
    assert _fastas(tmp_path / 'softmasked') == _fastas(tmp_path / 'unmasked')
    assert _fastas(tmp_path / 'unmasked')


@pytest.mark.parametrize('base', ['N', 'n', 'R', 'y'])
def test_slice_with_ambiguous_base_is_skipped(base, gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ any base outside ACGT in a slice skips the junction, as N did """

    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8, -: T1 when unmasked
    pos = (junction.alt1_es + junction.alt1_ee) // 2  # inside the included exon (1-based)
    record = genome.genome['15']
    mutated = dict(genome.genome)
    mutated['15'] = record[:pos - 1] + base + record[pos:]

    assert _run(junction, gtf, mutated, tmp_path / 'out', monkeypatch) == fates.SKIPPED_MASKED
