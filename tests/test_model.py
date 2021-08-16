# -*- coding: utf-8 -*-

""" Tests for read count model related functionalities """

import unittest
import os

from jcast import model
from jcast.junctions import RmatsResults


class  ModelTest(unittest.TestCase):
    """
    Test cases involving reading genome sequence and references
    """

    def setUp(self):
        # Create out folder if not exists
        if not os.path.exists('tests/out'):
            os.makedirs('tests/out')

        self.write_dir = 'tests/out'
        pass

    def tearDown(self):
        pass

    def test_gaussian_encode_human_heart(self):
        """ test for selecting read count cutoff in ENCODE human heart data """

        rmats_results = RmatsResults(rmats_dir='tests/data/model/encode_human_heart')
        rmats_results.get_junction_count_array()

        pt, gmm, min_count = model.gaussian_mixture(sum_sjc_array=rmats_results.sum_sjc_array)
        model.plot_model(sum_sjc_array=rmats_results.sum_sjc_array,
                         pt=pt, gmm=gmm, min_count=min_count,
                         write_dir=self.write_dir, filename='encode_human_heart')

        self.assertGreater(gmm.means_[1], gmm.means_[0])

    def test_gaussian_encode_mouse_heart(self):
        """ test for selecting read count cutoff in ENCODE mouse heart data """

        rmats_results = RmatsResults(rmats_dir='tests/data/model/encode_mouse_heart')
        rmats_results.get_junction_count_array()

        pt, gmm, min_count = model.gaussian_mixture(sum_sjc_array=rmats_results.sum_sjc_array)
        model.plot_model(sum_sjc_array=rmats_results.sum_sjc_array,
                         pt=pt, gmm=gmm, min_count=min_count,
                         write_dir=self.write_dir, filename='encode_mouse_heart')

        self.assertGreater(gmm.means_[1], gmm.means_[0])

    def test_gaussian_ipsc_human_d14(self):
        """ test for selecting read count cutoff in in house human iPSC cardiomyocyte data """

        rmats_results = RmatsResults(rmats_dir='tests/data/model/ipsc_human_d14')
        rmats_results.get_junction_count_array()

        pt, gmm, min_count = model.gaussian_mixture(sum_sjc_array=rmats_results.sum_sjc_array)
        model.plot_model(sum_sjc_array=rmats_results.sum_sjc_array,
                         pt=pt, gmm=gmm, min_count=min_count,
                         write_dir=self.write_dir, filename='ipsc_human_d14')

        self.assertGreater(gmm.means_[1], gmm.means_[0])