# -*- coding: utf-8 -*-

""" Tests for splice junction related functions """

import argparse
import types

import pytest

from jcast import fates, main, params
from jcast.junctions import Junction
from jcast.sequences import Sequence


@pytest.mark.network
def test_pkm_manual(gtf, genome, rmats_results):
    """ test for translating PKM, a known AS protein """

    rma = rmats_results.rmats_mxe
    # TODO: test negative strand and trimming
    i = 71  # PKM

    junction = Junction(**rma.iloc[i].to_dict())

    junction.trim_cds(gtf)
    junction.get_translated_phase(gtf)

    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)
    sequence.get_canonical_aa(gtf, genome.genome)
    sequence.translate(use_phase=True)

    for i in [1, 2]:
        sequence.stitch_to_canonical_aa(slice_to_stitch=i)

    # Compare with Uniprot retrieved sequences
    retrieved_pkm2 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                      'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                      'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                      'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                      'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                      'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                      'LSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEAS'
                      'FKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                      'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

    retrieved_pkm1 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                      'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                      'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                      'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                      'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                      'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                      'LSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEAS'
                      'YKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                      'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

    assert sequence.slice1_stitched.seq == retrieved_pkm2
    assert sequence.slice2_stitched.seq == retrieved_pkm1


def test_trim_start(gtf, rmats_results):
    """ test for trimming junctions """

    rma = rmats_results.rmats_a5ss
    # TODO: test negative strand and trimming
    i = 6  # SEMA4B A5SS 742

    junction = Junction(**rma.iloc[i].to_dict())
    junction.trim_cds(gtf)

    assert junction.tx0 == 90201579


def _junction(rma, rmats_id):
    return Junction(**rma[rma['id'] == rmats_id].iloc[0].to_dict())


def test_c_terminal_bound_plus_strand(gtf, rmats_results):
    """ tx1 is the end of the stop codon on the + strand (REVIEW.md §5.1) """

    junction = _junction(rmats_results.rmats_se, 6426)  # CASC4
    junction.trim_cds(gtf)

    assert (junction.tx0, junction.tx1) == (44289030, 44413406)
    assert junction.down_ee == 44413406


def test_c_terminal_bound_minus_strand(gtf, rmats_results):
    """ tx0 is the start of the stop codon on the - strand (REVIEW.md §5.1) """

    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8
    junction.trim_cds(gtf)

    assert (junction.tx0, junction.tx1) == (88899395, 88913319)
    assert junction.anc_es == 88899395


def test_c_terminal_event_translates(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ an event spanning the stop-codon exon is no longer read as a premature stop """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=False,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8

    fate = main._translate_one(junction, gtf, genome, args, pred_bound=args.read)

    assert fate.startswith('SUCCESS TIER 1')


def test_bounds_without_tsl_column(gtf, rmats_results):
    """ a GTF without transcript_support_level still yields bounds (REVIEW.md MED-7) """

    no_tsl = types.SimpleNamespace(annot=gtf.annot.drop(columns='transcript_support_level'))
    junction = _junction(rmats_results.rmats_se, 6426)  # CASC4
    junction.trim_cds(no_tsl)

    assert (junction.tx0, junction.tx1) == (44289030, 44413406)


def _exon(junction, exon):
    return getattr(junction, f'{exon}_es'), getattr(junction, f'{exon}_ee')


def test_trim_drops_utr_exon(gtf, rmats_results):
    """ exons lying entirely outside the coding region are removed (REVIEW.md PR 4a) """

    junction = _junction(rmats_results.rmats_se, 6443)  # TMEM62, +: anchor is 5' UTR
    junction.trim_cds(gtf)
    assert _exon(junction, 'anc') == (-1, -1)
    assert _exon(junction, 'alt1') == (43134257, 43134368)

    junction = _junction(rmats_results.rmats_se, 6654)  # UBE3A, -: downstream exon is 5' UTR
    junction.trim_cds(gtf)
    assert _exon(junction, 'down') == (-1, -1)
    assert _exon(junction, 'anc') == (25375465, 25375763)


def test_trim_minus_strand_alternative_sites(gtf, rmats_results):
    """ A5SS/A3SS on the - strand are trimmed at the bound their exons actually contain """

    junction = _junction(rmats_results.rmats_a5ss, 3645)  # WDR73: start codon in the alternative exons
    junction.trim_cds(gtf)
    assert _exon(junction, 'alt1') == (84654227, 84654274)
    assert _exon(junction, 'alt2') == (84654234, 84654274)

    junction = _junction(rmats_results.rmats_a3ss, 1238)  # HDDC3: stop codon in the alternative exons
    junction.trim_cds(gtf)
    assert _exon(junction, 'alt1') == (90931275, 90931556)
    assert _exon(junction, 'alt2') == (90931275, 90931405)


def test_trim_retained_intron_region(gtf, rmats_results):
    """ the retained region of an RI event is trimmed along with its flanking exons """

    junction = _junction(rmats_results.rmats_ri, 3355)  # B2M, +: stop codon in the downstream exon
    junction.trim_cds(gtf)
    assert _exon(junction, 'down') == (44716329, 44716342)
    assert _exon(junction, 'alt1') == (44715423, 44716342)


def test_event_outside_cds_is_skipped(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a splice event lying entirely in the UTR is skipped; the canonical is still written with -c """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=True,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_mxe, 433)  # ZNF592: both alternative exons are 5' UTR

    fate = main._translate_one(junction, gtf, genome, args, pred_bound=args.read)

    assert fate == fates.SKIPPED_OUTSIDE_CDS
    assert [p.name for p in tmp_path.glob('psq_*.fasta')] == ['psq_canonical.fasta']


def test_retained_intron_slices(genome, rmats_results):
    """ RI slice 1 is the retained exon alone; slice 2 splices the flanking exons (REVIEW.md HIGH-11) """

    junction = _junction(rmats_results.rmats_ri, 916)  # UNC45A, +
    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)

    chrom = genome.genome['15'].seq
    retained = chrom[junction.alt1_es - 1:junction.alt1_ee]
    spliced = chrom[junction.anc_es - 1:junction.anc_ee] + chrom[junction.down_es - 1:junction.down_ee]

    assert (len(sequence.slice1_nt), len(sequence.slice2_nt)) == (515, 37 + 324)
    assert sequence.slice1_nt.seq == retained
    assert sequence.slice2_nt.seq == spliced


def test_retained_intron_in_utr_is_skipped(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ an intron retained upstream of the start codon is outside the CDS; it used to yield a doubled slice """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=False,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_ri, 943)  # FAH, +

    assert main._translate_one(junction, gtf, genome, args, pred_bound=args.read) == fates.SKIPPED_OUTSIDE_CDS
