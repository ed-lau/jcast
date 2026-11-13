# -*- coding: utf-8 -*-

""" Tests for gtf related functionalities """


def test_human_tss(gtf):
    """ tests finding the translation start of a human gene """

    gtf0 = gtf.annot.query('gene_id == "ENSG00000104055"')
    assert len(gtf0) > 0

    # Get the coding sequence start codon locs, sort by transcript_support_level
    gtf0_start = gtf0.query('feature == "start_codon"').sort_values(['transcript_support_level']).sort_values(
        ['ccds_id']).loc[:, 'start']
    assert len(gtf0_start) > 0

    # Translation start site
    tx0 = gtf0_start.iloc[0]
    assert tx0 > 0


def test_rat_tss(rat_gtf):
    """ tests finding the translation start of a rat gene """

    gtf0 = rat_gtf.annot.query('gene_id == "ENSRNOG00000055204"')
    assert len(gtf0) > 0

    gtf0_start = gtf0.query('feature == "start_codon"').loc[:, 'start']
    assert len(gtf0_start) > 0

    # Translation start site
    tx0 = gtf0_start.iloc[0]
    assert tx0 > 0
