# -*- coding: utf-8 -*-

""" Tests for splice junction related functions """

import pytest

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
