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

    def test_pkm_manual(self):

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

        self.assertEqual(sequence.slice1_stitched.seq, retrieved_pkm2)
        self.assertEqual(sequence.slice2_stitched.seq, retrieved_pkm1)

    def test_trim_start(self):
        """ test for trimming """

        rma = rmats_results.rmats_a5ss
        # TODO: test negative strand and trimming
        i = 6  # SEMA4B A5SS 742

        junction = Junction(**rma.iloc[i].to_dict())
        junction.trim_cds(gtf)

        self.assertEqual(junction.tx0, 90201579)
        pass
