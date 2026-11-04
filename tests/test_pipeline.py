# -*- coding: utf-8 -*-

""" End-to-end golden test: run every chr15 junction through the tiering logic and pin the outcome.

Canonical sequences come from the GTF only and UniProt is stubbed out, so this test needs no network.
If an intentional change alters the output (e.g. REVIEW.md §5.1, the tx0/tx1 fix), update the expected
values below in the same commit and record the difference in CHANGELOG.md.
"""

import argparse
import collections
import hashlib

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import main, params
from jcast.junctions import Junction
from jcast.sequences import Sequence

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS

EXPECTED_FATES = {
    'FAILURE': 360,
    'PARTIAL TIER 4': 236,
    'SKIPPED': 297,
    'SUCCESS TIER 1': 185,
    'SUCCESS TIER 2': 25,
    'SUCCESS TIER 3': 2,
}

# file name -> (number of records, sha256 of the sorted sequences)
EXPECTED_FASTA = {
    'psq_T1.fasta': (152, '155d0579ce053722'),
    'psq_T1_orphan.fasta': (48, 'd706539ee9d4af17'),
    'psq_T2.fasta': (20, '6d3a178e69825516'),
    'psq_T2_orphan.fasta': (13, '67b969f4bb023a98'),
    'psq_T3_orphan.fasta': (4, '95f9daa6904c380d'),
    'psq_T4.fasta': (193, '40bceb3dd60b89fa'),
    'psq_T4_orphan.fasta': (72, '22c1f861b24e8af1'),
    'psq_canonical.fasta': (268, '6b891792cb75c69e'),
}


def _digest(path):
    seqs = sorted(str(r.seq) for r in SeqIO.parse(path, 'fasta'))
    return len(seqs), hashlib.sha256('\n'.join(seqs).encode()).hexdigest()[:16]


@pytest.fixture
def pipeline_output(gtf, genome, rmats_results, tmp_path, monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot',
                        lambda self, reviewed='true': SeqRecord(Seq('')))

    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=True,
                              mask='NnXx', qvalue=[0, 1])
    fates = collections.Counter()
    for rma in rmats_results.rmats_list:
        for i in range(len(rma)):
            junction = Junction(**rma.iloc[i].to_dict())
            fate = main._translate_one(junction, gtf, genome, args, pred_bound=args.read)
            fates[fate.split('.')[0]] += 1

    return fates, {p.name: _digest(p) for p in sorted(tmp_path.glob('psq_*.fasta'))}


@pytest.mark.slow
def test_chr15_golden(pipeline_output):
    fates, fasta = pipeline_output

    assert dict(fates) == EXPECTED_FATES
    assert fasta == EXPECTED_FASTA


def test_cli_rejects_missing_files(monkeypatch, capsys):
    monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), 'missing.gtf', str(HUMAN_GENOME)])

    with pytest.raises(SystemExit):
        main.main()
    assert 'file not found: missing.gtf' in capsys.readouterr().err


def test_read_gzipped_genome(genome):
    """ REVIEW.md HIGH-8: a .fa.gz genome used to crash with UnicodeDecodeError. """
    assert len(genome.genome['15']) > 100_000_000


def test_cli_passes_paths(monkeypatch, tmp_path):
    """ GTF and genome reach run_jcast as plain paths, not open file handles. """
    calls = []
    monkeypatch.setattr(main, 'run_jcast', lambda args: calls.append(args))
    monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), str(HUMAN_GTF), str(HUMAN_GENOME),
                                     '-o', str(tmp_path)])

    main.main()

    assert calls[0].genome == str(HUMAN_GENOME)
    assert isinstance(calls[0].gtf_file, str)
