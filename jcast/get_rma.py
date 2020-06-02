# -*- coding: utf-8 -*-

""" Methods that concern RMATS results - reading the files and mapping the columns to exons. """


import logging
import os.path
import pandas as pd

class RmatsResults(object):
    """
    Container to hold the rMATS output folder and create individual objects based on the five splice types.

    Note that each splice event will be defined by the anchor exon (anc), potential alternative exons (alt1, alt2), and
    the downstream exon (down). For some splice type, some of these will not be present, for example, Skipped Exons
    either contain the alt1 exon or no alt1 exon in the two alternative forms, and no alt2 is present.

    """

    def __init__(self, rmats_dir):

        self.dir = rmats_dir
        self.rmats_mxe = self._read_rmats_mxe()
        self.rmats_se = self._read_rmats_se()
        self.rmats_ri = self._read_rmats_ri()
        self.rmats_a5ss = self._read_rmats_a5ss()
        self.rmats_a3ss = self._read_rmats_a3ss()

        self.logger = logging.getLogger('jcast.input')

    def _read_rmats_mxe(self):
        """
        Read input data frame for Mutually Exclusive Exons (MXE)
        The slices should be anc-alt1-down and anc-alt2-down

        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'MXE.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                      'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'anc_es',
                      'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                      'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                      'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'MXE'

        return df

    def _read_rmats_se(self):
        """
        Read input data frame for Skipped Exons (SE)
        The slices should be anc-down and anc-alt1-down

        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'SE.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'SE'
        df['alt2_es'] = -1
        df['alt2_ee'] = -1

        return df

    def _read_rmats_ri(self):
        """
        Read input data frame for Retained Introns (RI)
        The slices should be anc-down and anc-alt1-down

        :return:
        """
        df = pd.read_table(os.path.join(self.dir, 'RI.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'RI'
        df['alt2_es'] = -1
        df['alt2_ee'] = -1

        return df

    def _read_rmats_a5ss(self):
        """
        Read input data frame for Alternative 5' Splice Sites (A5SS)
        Note this splice type is without the 'downstream' exon, but the anchor (flanking) is downstream.
        The slices should be alt1-anc and alt2-anc
        Note if the strand is +, alt1_es and alt2_es should be identical and before anchor in genomic position.
        Note if the strand is -, alt1_ee and alt2_ee should be the same and after anchor in genomic position.
        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'A5SS.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'A5SS'
        df['down_es'] = -1
        df['down_ee'] = -1

        return df

    def _read_rmats_a3ss(self):
        """
        Read input data frame for Alternative 3' Splice Sites (A3SS)
        Note this splice type is without the downstream exon, the anchor is the upstream.
        The slices are anc-alt1 and anc-alt2.
        Note if the strand is +, alt1_ee and alt2_ee are identical and after anchor in genomic position.
        If the stand is -, alt1_es and alt2_es are the same, and they are both before anchor in genomic position.
        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'A3SS.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'A3SS'
        df['down_es'] = -1
        df['down_ee'] = -1

        return df