# -*- coding: utf-8 -*-

""" Tests """


import unittest
import os
from jcast.annots import ReadAnnotations


class AnnotationTest(unittest.TestCase):
    """
    Test cases involving reading genome sequence and references
    """

    def setUp(self):
        """

        :return:
        """

        global test_data_loc, human_gtf_loc, rat_gtf_loc


        test_data_loc = os.path.join('tests', 'data')
        human_gtf_loc = os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf')
        rat_gtf_loc = os.path.join(test_data_loc, 'genome', 'Rattus_norvegicus.Rnor_6.0.100.chromosome.12.gtf')

        # Remove the cached gtf file
        try:
            os.remove(human_gtf_loc + '.cached')
        except FileNotFoundError:
            pass

        try:
            os.remove(rat_gtf_loc + '.cached')
        except FileNotFoundError:
            pass

    def tearDown(self):

        # Remove the cached gtf file
        try:
            os.remove(human_gtf_loc + '.cached')
        except FileNotFoundError:
            pass

        try:
            os.remove(rat_gtf_loc + '.cached')
        except FileNotFoundError:
            pass

        pass

    def test_human_tss(self):
        """ finds the translation start of a human gene """

        human_gtf = ReadAnnotations(human_gtf_loc)
        human_gtf.read_gtf()

        gtf0 = human_gtf.annot.query('gene_id == "ENSG00000104055"')
        self.assertGreater(len(gtf0), 0)

        # Get the coding sequence start codon locs, sort by transcript_support_level
        gtf0_start = gtf0.query('feature == "start_codon"').sort_values(['transcript_support_level']).sort_values(
            ['ccds_id']).loc[:, 'start']
        self.assertGreater(len(gtf0_start), 0)

        # Translation start site
        tx0 = gtf0_start.iloc[0]
        self.assertGreater(tx0, 0)

        # Get anchor frame
        gtf0 = human_gtf.annot.query('gene_id == "ENSG00000259490"').query('start == "19987656"'). \
            query('end == "19987968"').query('feature == "CDS"')

        pass

    def test_rat_tss(self):
        """ finds the translation start of a rat gene """

        rat_gtf = ReadAnnotations(rat_gtf_loc)
        rat_gtf.read_gtf()

        gtf0 = rat_gtf.annot.query('gene_id == "ENSRNOG00000055204"')
        self.assertGreater(len(gtf0), 0)

        gtf0_start = gtf0.query('feature == "start_codon"').loc[:, 'start']
        self.assertGreater(len(gtf0_start), 0)

        # Translation start site
        tx0 = gtf0_start.iloc[0]
        self.assertGreater(tx0, 0)

        pass


