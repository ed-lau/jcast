# -*- coding: utf-8 -*-

""" Tests """


import unittest
import os


from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence

class GenomeTest(unittest.TestCase):
    """
    Test cases involving reading genome sequence and references
    """

    def setUp(self):
        """

        :return:
        """

        global gtf, rmats_results, genome, test_data_loc

        test_data_loc = os.path.join('tests', 'data')
        genome_loc = os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
        gtf_loc = os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf')
        rmats_loc = os.path.join(test_data_loc, 'rmats')

        genome = ReadGenome(genome_loc)

        gtf = ReadAnnotations(gtf_loc)
        gtf.read_gtf()

        rmats_results = RmatsResults(rmats_loc)

        pass

    def tearDown(self):

        # Remove the cached gtf file
        os.remove(os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf.cached'))

        pass


    def test_retrieve_gtf_canonical_transcripts(self):
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

        self.assertEqual(gtf_aa.seq, uniprot_aa.seq)


    def test_unannotated_transcripts(self):
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
        self.assertEqual(len(sequence.canonical_aa), 0)

        # since there is no canonical, there should be no stitching
        sequence.stitch_to_canonical_aa(slice_to_stitch=1)
        self.assertIsNone(sequence.slice1_stitched)


