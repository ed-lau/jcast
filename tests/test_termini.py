# -*- coding: utf-8 -*-

""" Alternative C- and N-termini (PR 4c; reports/2026-11-10-alt-termini-audit.md) """

from Bio import SeqIO

from jcast import fates

from tests.test_junction import _annotated_protein, _junction, _translate


def _records(out_dir, name):
    path = out_dir / f'psq_{name}.fasta'
    return {str(r.seq) for r in SeqIO.parse(path, 'fasta')} if path.exists() else set()


def test_two_alternative_c_termini(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ both slices end at different annotated stop codons: FAILURE before, both isoforms reconstructed now """

    junction = _junction(rmats_results.rmats_se, 28030)  # FAM96A, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 1 ALT TERMINUS')
    written = _records(tmp_path, 'T1_alt_cterm')
    assert {_annotated_protein(gtf, genome, t) for t in ['FAM96A-206', 'FAM96A-202']} <= written


def test_alternative_c_terminus_and_canonical(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ one slice is the canonical, the other ends at an alternative annotated stop codon (T4 before) """

    junction = _junction(rmats_results.rmats_a3ss, 6113)  # NMB, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 1 ALT TERMINUS')
    assert _annotated_protein(gtf, genome, 'NMB-201') in _records(tmp_path, 'T1_alt_cterm')
    assert _annotated_protein(gtf, genome, 'NMB-202') in _records(tmp_path, 'canonical')


def test_terminal_annotated_stop_is_an_alternative_c_terminus(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a slice ending exactly at an annotated stop that is not the canonical's is joined N-terminally only """

    junction = _junction(rmats_results.rmats_se, 6592)  # EMC4, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 1 ALT TERMINUS')
    assert _annotated_protein(gtf, genome, 'EMC4-201') in _records(tmp_path, 'T1_alt_cterm')


def test_alternative_n_terminus(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a slice lacking the gene's start codon is translated from its own annotated start, joined C-terminally """

    junction = _junction(rmats_results.rmats_se, 5744)  # DET1, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 1 ALT TERMINUS')
    assert _annotated_protein(gtf, genome, 'DET1-201') in _records(tmp_path, 'T1_alt_nterm')


def test_alternative_n_terminus_replaces_recreated_canonical(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ ARPP19 SE 7770 was T1 with its start-less slice recreating the canonical; now it is ARPP19-202 """

    junction = _junction(rmats_results.rmats_se, 7770)  # ARPP19, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 1 ALT TERMINUS')
    assert _annotated_protein(gtf, genome, 'ARPP19-202') in _records(tmp_path, 'T1_alt_nterm')
    assert _annotated_protein(gtf, genome, 'ARPP19-201') in _records(tmp_path, 'T1')


def test_nmd_only_stop_is_still_a_ptc(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a stop codon annotated only on nonsense-mediated decay transcripts is not an alternative C-terminus """

    junction = _junction(rmats_results.rmats_se, 5703)  # SEMA4B, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('PARTIAL TIER 4')
    assert not list(tmp_path.glob('psq_T1_alt*.fasta'))


def test_anchored_frameshift_stays_tier_2(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ slices ending at the same stop in different frames are a frameshift even if one is anchored """

    junction = _junction(rmats_results.rmats_mxe, 3368)  # APH1B, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('SUCCESS TIER 2')


def test_internal_event_ending_at_the_canonical_stop_is_unchanged(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ both slices end at the canonical's stop codon: plain T1, as before """

    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate == fates.TIER1.format(junction.phase_label, junction.phase_label)
    assert _annotated_protein(gtf, genome, 'MFGE8-202') in _records(tmp_path, 'T1')
