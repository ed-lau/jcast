# -*- coding: utf-8 -*-

""" Methods that concern splice junctions - getting their coordinates, transcription starts/ends, phases. """


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

    def __init__(self,
                 rmats_dir,
                 ):

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


class Junction(object):
    """
    Splice junctions (one row in an rMTAS output file), containing gene ID, and exon coordinates. There are five types
    of exon splice junctions, as from the five rMATS results files.

    """

    def __init__(self, **kwargs):
        self.name = str(kwargs['id'])
        self.gene_id = kwargs['gene_id']
        self.strand = kwargs['strand']
        self.chr = kwargs['chr']
        self.anc_es = kwargs['anc_es']+1
        self.anc_ee = kwargs['anc_ee']
        self.alt1_es = kwargs['alt1_es']+1
        self.alt1_ee = kwargs['alt1_ee']
        self.alt2_es = kwargs['alt2_es']+1
        self.alt2_ee = kwargs['alt2_ee']
        self.down_es = kwargs['down_es']+1
        self.down_ee = kwargs['down_ee']
        self.junction_type = kwargs['jxn_type']
        self.gene_symbol = kwargs['gene_symbol']
        self.tx1 = -1
        self.tx0 = -1
        self.phase = -1
        self.min_read_count = 0
        self.fdr = kwargs['fdr']
        self.sjc_s1 = kwargs['sjc_s1']
        self.sjc_s2 = kwargs['sjc_s2']

        self.logger = logging.getLogger('jcast.junction')


    def __repr__(self):
        """ repr """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    def __str__(self):
        """ str """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    def set_min_read_count(self, count):
        """
        A setter to mark the minimum read count in the junction

        :return:
        """
        self.min_read_count = count

        return True

    def get_translated_region(self, gtf):
        """
        Read the genome annotation (.gtf) file, find the coding sequences (CDS) that share the gene name and coordinates
        of the anchor exon (anc) supplied, then find out where the annotated translation start and end sites are. If the
        splice junction sequences extend BEYOND the translation starts and ends, trim them to avoid running into stop
        codons.

        :param gtf: Genome annotation
        :return: True
        """
        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')

        #
        # 	Get the translation start and end positions
        #
        try:
            gtf0_start = gtf0.query('feature == "start_codon"').sort_values(['transcript_support_level']).sort_values(['ccds_id']).loc[:, 'start']

        # Bypass transcript support level in non-human/mouse Ensembl gtfs
        except KeyError:
            gtf0_start = gtf0.query('feature == "start_codon"').loc[:, 'start']

        if len(gtf0_start) > 0:
            self.tx0 = gtf0_start.iloc[0]
        else:
            self.tx0 = -1

        try:
            gtf0_end = gtf0.query('feature == "stop_codon"').sort_values(['transcript_support_level']).loc[:, 'start']

        # Bypass transcript support level in non-human/mouse Ensembl gtfs
        except KeyError:
            gtf0_end = gtf0.query('feature == "stop_codon"').loc[:, 'start']

        if len(gtf0_end) > 0:
            self.tx1 = gtf0_end.iloc[0]
        else:
            self.tx1 = -1

        #
        # By rMATS convention, if strand is -ve
        # then the upstream is near the end of tx
        #
        if self.strand == '-' and self.tx1 > 0 and self.tx0 > 0:
            self.tx0, self.tx1 = self.tx1, self.tx0
            self.tx1 += 2

        elif self.strand == '+':
            self.tx1 += 2


        return True

    def get_translated_phase(self, gtf):
        """
        Get the annotated translation phase from the GTF file
        :param gtf:  genome annotation
        :return:
        """

        # Select the anchor exon from CDS and get the frame
        # Note 2018-03-24 this is probably not quite right. I think you should find out whether the anchor
        # is really the one that determines the phase here.

        # 2018-03-24: First define the exon we are looking for.
        if self.junction_type in ['MXE', 'SE', 'RI']:
            if self.strand == '+':
                ph0, ph1 = self.anc_es, self.anc_ee
            elif self.strand == '-':
                ph0, ph1 = self.down_es, self.down_ee

        elif self.junction_type == 'A5SS':
            if self.strand == '+':
                ph0, ph1 = self.alt1_ee, self.alt1_es  # Might have to search also for alt2
            elif self.strand == '-':
                ph0, ph1 = self.anc_es, self.anc_ee

        elif self.junction_type == 'A3SS':
            if self.strand == '+':
                ph0, ph1 = self.anc_ee, self.anc_es  # Might have to search also for alt2
            elif self.strand == '-':
                ph0, ph1 = self.alt1_es, self.alt2_ee

        self.logger.info('Anchor start {0} end {1}'.format(ph0,
                                                           ph1))

        # Get the frame of that coding exon from GTF.
        gtf0 = gtf.annot.query('gene_id == @self.gene_id').query('start == @ph0').\
            query('end == @ph1').query('feature == "CDS"')

        # If phases retrieved, get the first value
        if len(gtf0) > 0:
            self.phase = int([x for x in gtf0.loc[:, 'frame'].iloc if x != '.'][0])

        else:
            self.phase = -1

        return True

    def trim(self):
        """
        :return: True

        Trims the junction based on transcription start and end:

        """
        if self.junction_type == 'MXE':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.anc_es = -1
                    self.anc_ee = -1
                    self.alt1_es = self.tx0

                if self.alt2_ee > self.tx0 > self.alt2_es:
                    self.anc_es = -1
                    self.anc_ee = -1
                    self.alt2_es = self.tx0
            except:
                self.logger.info('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.down_es = -1
                    self.down_ee = -1
                    self.alt1_ee = self.tx1

                if self.alt2_ee > self.tx1 > self.alt2_es:
                    self.down_es = -1
                    self.down_ee = -1
                    self.alt2_ee = self.tx1
            except:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'SE':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.anc_es = -1
                    self.anc_ee = -1
                    self.alt1_es = self.tx0

            except:
                self.logger.info('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.down_es = -1
                    self.down_ee = -1
                    self.alt1_ee = self.tx1

            except:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'RI':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except:
                self.logger.info('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

            except:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'A5SS':
            try:
                if self.alt2_ee > self.tx0 > self.alt2_es:
                    self.alt2_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.alt1_es = self.tx0

            except:
                self.logger.info('Trimming start failed.')

            try:
                if self.anc_ee > self.tx1 > self.anc_es:
                    self.anc_ee = self.tx1

            except:
                self.logger.info('Trimming end failed.')


        elif self.junction_type == 'A3SS':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except:
                self.logger.info('Trimming start failed.')

            try:
                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.alt1_ee = self.tx1

                if self.alt2_ee > self.tx1 > self.alt2_es:
                    self.alt2_ee = self.tx1

            except:
                self.logger.info('Trimming end failed.')

        return True
