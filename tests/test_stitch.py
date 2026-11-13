# -*- coding: utf-8 -*-

""" Joining translated slices to the canonical (PR 5; reports/2026-11-11-stitch-anchor-audit.md) """

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import fates, params
from jcast.sequences import Sequence

from tests.test_junction import _junction, _translate


def _records(out_dir, name, tag):
    """ sequences in psq_<name>.fasta whose header carries tag, e.g. '|SE2|5917|' """
    path = out_dir / f'psq_{name}.fasta'
    return [str(r.seq) for r in SeqIO.parse(path, 'fasta') if tag in r.id] if path.exists() else []


def _canonical(junction, gtf, genome, monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    sequence = Sequence(junction)
    return str(sequence.get_canonical_aa_gtf(gtf, genome.genome).seq)


def test_peptide_ending_at_a_stop_is_not_joined_c_terminally(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a T4 peptide ending at its premature stop was joined at both ends, giving back the canonical """

    junction = _junction(rmats_results.rmats_se, 19972)  # UBE2Q2, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)
    canonical = _canonical(_junction(rmats_results.rmats_se, 19972), gtf, genome, monkeypatch)

    assert fate.startswith('PARTIAL TIER 4')
    [product] = _records(tmp_path, 'T4', '|SE1|19972|')
    assert product != canonical
    assert len(product) == 197 and canonical.startswith(product[:149])


def test_peptide_starting_at_the_gene_start_is_not_joined_n_terminally(gtf, genome, rmats_results, tmp_path,
                                                                        monkeypatch):
    """ COMMD4 exon 2 skipped: MSSVKL... also occurs at canonical Met24, which prepended canonical residues 1-24 """

    junction = _junction(rmats_results.rmats_se, 5917)  # COMMD4, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('PARTIAL TIER 4')
    [product] = _records(tmp_path, 'T4', '|SE2|5917|')
    assert product.startswith('MSSVKL') and len(product) == 175


def test_repeated_anchor_is_not_joined(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ CD276 (duplicated Ig domains): the start anchor occurs twice; the first occurrence gave a wrong product """

    junction = _junction(rmats_results.rmats_a3ss, 1316)  # CD276, +
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate == fates.TIER1.format(junction.phase_label, junction.phase_label)
    assert not _records(tmp_path, 'T1', '|A3SS1|1316|')
    assert _records(tmp_path, 'T1_orphan', '|A3SS1|1316|')


def test_join_rule_on_synthetic_canonical():
    """ an anchor is used only if it occurs once, and the start anchor must not be after the end anchor """

    sequence = Sequence.__new__(Sequence)
    sequence.canonical_aa = SeqRecord(Seq('MDEFGHIKLAAAAAAAPEPTIDEWCCCCCCCRSTVWYAAAAAAAQ'))
    length = params.aa_stitch_length

    # unique anchors: joined at both ends
    assert str(sequence._join('PEPTIDENEWCCCCCCC', 'both').seq) == 'MDEFGHIKLAAAAAAAPEPTIDENEWCCCCCCCRSTVWYAAAAAAAQ'
    # the start anchor (AAAAAAA) occurs twice: not joined
    assert sequence._join('A' * length + 'NEW' + 'RSTVWY', 'both') is None
    # the start anchor is after the end anchor: not joined at both ends, but N-terminally if allowed
    assert sequence._join('RSTVWYA' + 'NEW' + 'PEPTIDE', 'both') is None
    assert str(sequence._join('RSTVWYA' + 'NEW' + 'PEPTIDE', 'both_or_n').seq) == \
        'MDEFGHIKLAAAAAAAPEPTIDEWCCCCCCC' + 'RSTVWYA' + 'NEW' + 'PEPTIDE'


def test_short_peptide_is_joined(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a 7-residue peptide could not be joined with the former 10-residue anchor; joined, it is the canonical """

    junction = _junction(rmats_results.rmats_se, 11950)  # GTF2A2, -
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)
    canonical = _canonical(_junction(rmats_results.rmats_se, 11950), gtf, genome, monkeypatch)

    assert fate.startswith('SUCCESS TIER 2')
    assert not _records(tmp_path, 'T2_orphan', '|SE2|11950|')
    assert _records(tmp_path, 'canonical', '') == [canonical]
