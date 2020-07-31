# -*- coding: utf-8 -*-

""" Tests """


import unittest
import jcast.helpers as h

class GenomeTest(unittest.TestCase):
    """
    Test cases involving reading genome sequence and references
    """

    def setUp(self):
        """

        :return:
        """

        pass

    def tearDown(self):

        pass

    def test_translate_cdna(self):

        # Insulin 4 transcript protein coding region and protein sequences downloaded
        # http://uswest.ensembl.org/Homo_sapiens/Transcript/Sequence_Protein?db=core;g=ENSG00000120211;r=9:5231419-5235304;t=ENST00000239316
        ins4_cdna = ('ATGGCCAGCCTGTTCCGGTCCTATCTGCCAGCAATCTGGCTGCTGCTGAGCCAACTCCTT'
                     'AGAGAAAGCCTAGCAGCAGAGCTGAGGGGATGTGGTCCCCGATTTGGAAAACACTTGCTG'
                     'TCATATTGCCCCATGCCTGAGAAGACATTCACCACCACCCCAGGAGGGTGGCTGCTGGAA'
                     'TCTGGACGTCCCAAAGAAATGGTGTCAACCTCCAACAACAAAGATGGACAAGCCTTAGGT'
                     'ACGACATCAGAATTCATTCCTAATTTGTCACCAGAGCTGAAGAAACCACTGTCTGAAGGG'
                     'CAGCCATCATTGAAGAAAATAATACTTTCCCGCAAAAAGAGAAGTGGACGTCACAGATTT'
                     'GATCCATTCTGTTGTGAAGTAATTTGTGACGATGGAACTTCAGTTAAATTATGTACATAG')

        ins4_protein = ('MASLFRSYLPAIWLLLSQLLRESLAAELRGCGPRFGKHLLSYCPMPEKTFTTTPGGWLLE'
                        'SGRPKEMVSTSNNKDGQALGTTSEFIPNLSPELKKPLSEGQPSLKKIILSRKKRSGRHRF'
                        'DPFCCEVICDDGTSVKLCT')

        self.assertEqual(h.make_pep(ins4_cdna, strand='+', phase=0, terminate=True),
                         ins4_protein)

        # TODO: translate all phases/strands

        pass