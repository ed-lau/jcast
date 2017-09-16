#
# 	Classes that concern RMATS results - reading the files and mapping the columns to exons
#   (anc, alt1, alt2, down)
#

class RmatsResults(object):

    def __init__(self, dir):

        import os.path
        import pandas as pd

        self.rmats_mxe = pd.read_table(os.path.join(dir, 'MXE.MATS.JC.txt'), sep='\t')
        self.rmats_mxe.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                  'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'anc_es',
                                  'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                  'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                  'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_mxe['jxn_type'] = 'MXE'

        self.rmats_se = pd.read_table(os.path.join(dir, 'SE.MATS.JC.txt'), sep='\t')
        self.rmats_se.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_se['jxn_type'] = 'SE'
        self.rmats_se['alt2_es'] = -1
        self.rmats_se['alt2_ee'] = -1

        self.rmats_ri = pd.read_table(os.path.join(dir, 'RI.MATS.JC.txt'), sep='\t')
        self.rmats_ri.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_ri['jxn_type'] = 'RI'
        self.rmats_ri['alt2_es'] = -1
        self.rmats_ri['alt2_ee'] = -1

        self.rmats_a5ss = pd.read_table(os.path.join(dir, 'A5SS.MATS.JC.txt'), sep='\t')
        self.rmats_a5ss.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'anc_es',
                                   'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_a5ss['jxn_type'] = 'A5SS'
        self.rmats_a5ss['alt2_es'] = -1
        self.rmats_a5ss['alt2_ee'] = -1

        self.rmats_a3ss = pd.read_table(os.path.join(dir, 'A3SS.MATS.JC.txt'), sep='\t')
        self.rmats_a3ss.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        self.rmats_a3ss['jxn_type'] = 'A3SS'
        self.rmats_a3ss['down_es'] = -1
        self.rmats_a3ss['down_ee'] = -1
