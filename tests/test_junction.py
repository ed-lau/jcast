# -*- coding: utf-8 -*-

""" Tests for splice junction related functions """

import argparse
import types

import pytest

from jcast import main, params
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
