# -*- coding: utf-8 -*-

""" End-to-end golden test: run every chr15 junction through the tiering logic and pin the outcome.

Canonical sequences come from the GTF only and UniProt is stubbed out, so this test needs no network.
If an intentional change alters the output, update the expected values below in the same commit and record
the difference in CHANGELOG.md. History:
    PR 4 (tx0/tx1 from stop codons) changed 48 fates; reports/2026-11-03-tx1-translation-bounds-audit.md
    PR 4a (trim every exon to the CDS, skip events outside it) changed 103; reports/2026-11-10-utr-exon-trimming-audit.md
    PR 4a2 (RI slice 1 is the retained exon alone) changed 5; reports/2026-11-10-ri-slice-audit.md
    PR 4b candidate (per-slice frame, fallback kept) changed 73; reports/2026-11-10-phase-audit.md
    PR 4b final, 2679edb (three-frame fallback and Tier 3 removed) changed 15 more; same report, §9
    PR 4b2, a0f3338 (frameshift test from coding lengths) changed 17; reports/2026-11-10-frameshift-audit.md
    PR 4c, e11cc86 (alternative C- and N-termini) changed 26; reports/2026-11-10-alt-termini-audit.md §10
    PR 5, 2142711 (joining to the canonical: complete ends, unique anchors, anchor length 7) changed 1 fate and moved
        records between joined and orphan files; reports/2026-11-11-stitch-anchor-audit.md §9
    PR 6 (keep a junction if its read count >= -r; -r 1 now keeps every junction) translates the 297 junctions with 1
        read; reports/2026-11-12-thresholds-runs-canonical-audit.md §9
"""

import argparse
import collections
import hashlib

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import fates, main, params
from jcast.junctions import Junction
from jcast.sequences import Sequence

from tests.conftest import HUMAN_GENOME, HUMAN_GTF, RMATS

EXPECTED_FATES = {
    'FAILURE': 207,
    'PARTIAL TIER 4': 403,
    'SKIPPED: Splice event lies outside the annotated coding region': 119,
    'SUCCESS TIER 1': 301,
    'SUCCESS TIER 1 ALT TERMINUS': 26,
    'SUCCESS TIER 2': 49,
}

# file name -> (number of records, sha256 of the sorted sequences)
EXPECTED_FASTA = {
    'psq_T1.fasta': (244, 'a1a57325f866167f'),
    'psq_T1_alt_cterm.fasta': (18, 'a03a3ea8d813b497'),
    'psq_T1_alt_cterm_orphan.fasta': (2, 'af4d8be95cb5dd05'),
    'psq_T1_alt_nterm.fasta': (6, '8ec341ee7237f78b'),
    'psq_T1_orphan.fasta': (38, '87b70b2314626c70'),
    'psq_T2.fasta': (49, 'a29140a2da79e90e'),
    'psq_T2_orphan.fasta': (6, '105037b254205e97'),
    'psq_T4.fasta': (350, '37e01ff952fcd5eb'),
    'psq_T4_orphan.fasta': (82, '0d2597931bd80267'),
    'psq_canonical.fasta': (266, '2d5de4ca94513af5'),
}


def _label(fate):
    """ fate class, keeping the reason for skips """
    head = fate.split('.')
    return f'{head[0]}:{head[1]}' if head[0] == 'SKIPPED' else head[0]


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
            fates[_label(fate)] += 1

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


def _junction_with_count(rmats_results, count):
    rma = rmats_results.rmats_se
    return next(j for j in (Junction(**rma.iloc[i].to_dict()) for i in range(len(rma))) if j.min_jc == count)


@pytest.mark.parametrize('read, skipped', [(2, False), (3, True)])
def test_read_cutoff_keeps_a_count_equal_to_it(gtf, genome, rmats_results, tmp_path, monkeypatch, read, skipped):
    """ -r N keeps a junction with N reads; it used to need N + 1 (REVIEW.md MED-11) """
    monkeypatch.setattr(params, 'use_gtf_only', True)
    junction = _junction_with_count(rmats_results, 2)
    args = argparse.Namespace(out=str(tmp_path), read=read, model=False, canonical=False, mask=None, qvalue=[0, 1])

    fate = main._translate_one(junction, gtf, genome, args, pred_bound=read)

    assert (fate == fates.SKIPPED_LOW_COUNT) == skipped


def test_filters_come_first(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a skipped junction is not trimmed or translated; with -c its gene's canonical is still written """
    monkeypatch.setattr(params, 'use_gtf_only', True)
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', lambda self, reviewed='true': SeqRecord(Seq('')))
    junction = _junction_with_count(rmats_results, 2)
    args = argparse.Namespace(out=str(tmp_path), read=3, model=False, canonical=True, mask=None, qvalue=[0, 1])

    assert main._translate_one(junction, gtf, genome, args, pred_bound=3) == fates.SKIPPED_LOW_COUNT
    assert junction.tx0 is None and junction.slice_phases == (None, None)
    assert len(list(SeqIO.parse(tmp_path / 'psq_canonical.fasta', 'fasta'))) == 1

    args.qvalue = [0, 0.0001]
    junction = _junction_with_count(rmats_results, 2)
    assert main._translate_one(junction, gtf, genome, args, pred_bound=1) == fates.SKIPPED_P_VALUE
    assert junction.tx0 is None
