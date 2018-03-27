#
# 	Classes that concern RMATS results - reading the files and mapping the columns to exons
#   (anc, alt1, alt2, down)
#

class RmatsResults(object):
    """
    Container to hold the rMATS output folder and create individual objects based on the five splice types.

    Note that each splice event will be defined by the anchor exon (anc), potential alternative exons (alt1, alt2), and
    the downstream exon (down). For some splice type, some of these will not be present, for example, Skipped Exons
    either contain the alt1 exon or no alt1 exon in the two alternative forms, and no alt2 is present.

    """

    def __init__(self, dir):

        import os.path
        import pandas as pd

        # Code for Mutually Exclusive Exons (MXE)
        # The slices should be anc-alt1-down and anc-alt2-down
        self.rmats_mxe = pd.read_table(os.path.join(dir, 'MXE.MATS.JC.txt'), sep='\t')
        self.rmats_mxe.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                  'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'anc_es',
                                  'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                  'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                  'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_mxe['jxn_type'] = 'MXE'

        # Code for Skipped Exons (SE)
        # The slices should be anc-down and anc-alt1-down
        self.rmats_se = pd.read_table(os.path.join(dir, 'SE.MATS.JC.txt'), sep='\t')
        self.rmats_se.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_se['jxn_type'] = 'SE'
        self.rmats_se['alt2_es'] = -1
        self.rmats_se['alt2_ee'] = -1

        # Code for Retained Introns (RI)
        # The slices should be anc-down and anc-alt1-down
        self.rmats_ri = pd.read_table(os.path.join(dir, 'RI.MATS.JC.txt'), sep='\t')
        self.rmats_ri.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_ri['jxn_type'] = 'RI'
        self.rmats_ri['alt2_es'] = -1
        self.rmats_ri['alt2_ee'] = -1

        # Code for Alternative 5' Splice Sites (A5SS)
        # Note this splice type is without the 'downstream' exon, but the anchor (flanking) is downstream.
        # The slices should be alt1-anc and alt2-anc
        # Note if the strand is +, alt1_es and alt2_es should be identical and before anchor in genomic position.
        # Note if the strand is -, alt1_ee and alt2_ee should be the same and after anchor in genomic position.
        self.rmats_a5ss = pd.read_table(os.path.join(dir, 'A5SS.MATS.JC.txt'), sep='\t')
        self.rmats_a5ss.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_a5ss['jxn_type'] = 'A5SS'
        self.rmats_a5ss['down_es'] = -1
        self.rmats_a5ss['down_ee'] = -1

        # Code for Alternative 3' Splice Sites (A3SS)
        # Note this splice type is without the downstream exon, the anchor is the upstream.
        # The slices are anc-alt1 and anc-alt2.
        # Note if the strand is +, alt1_ee and alt2_ee are identical and after anchor in genomic position.
        # If the stand is -, alt1_es and alt2_es are the same, and they are both before anchor in genomic position.
        self.rmats_a3ss = pd.read_table(os.path.join(dir, 'A3SS.MATS.JC.txt'), sep='\t')
        self.rmats_a3ss.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_a3ss['jxn_type'] = 'A3SS'
        self.rmats_a3ss['down_es'] = -1
        self.rmats_a3ss['down_ee'] = -1
