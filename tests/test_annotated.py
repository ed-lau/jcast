# -*- coding: utf-8 -*-

""" Tests for genome related functionalities """

import pytest

from jcast.junctions import Junction
from jcast.sequences import Sequence


@pytest.mark.network
def test_retrieve_gtf_canonical_transcripts(gtf, genome, rmats_results):
    """ test for DPP8 between GTF and Uniprot canonical sequences """

    rma = rmats_results.rmats_se
    # TODO: test negative strand and trimming
    i = 157  # DPP8

    junction = Junction(**rma.iloc[i].to_dict())
    junction.trim_cds(gtf)
    junction.get_translated_phase(gtf)
    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)

    # Get the GTF aa and Uniprot AA
    gtf_aa = sequence.get_canonical_aa_gtf(gtf,
                                           genome_index=genome.genome)

    uniprot_aa = sequence.get_canonical_aa_uniprot()

    assert gtf_aa.seq == uniprot_aa.seq


@pytest.mark.network
def test_unannotated_transcripts(gtf, genome, rmats_results):
    """ test for GOLGA2P10 which has no annotated translation product on gtf or uniprot """

    rma = rmats_results.rmats_mxe
    # TODO: test negative strand and trimming
    i = 8  # GOLPA2P10

    junction = Junction(**rma.iloc[i].to_dict())
    junction.trim_cds(gtf)
    junction.get_translated_phase(gtf)
    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)
    sequence.get_canonical_aa(gtf=gtf, genome_index=genome.genome)

    # return a 0 length sequence since no canonical annotation
    assert len(sequence.canonical_aa) == 0

    # since there is no canonical, there should be no stitching
    sequence.stitch_to_canonical_aa(slice_to_stitch=1)
    assert sequence.slice1_stitched is None
