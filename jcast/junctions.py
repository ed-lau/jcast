# -*- coding: utf-8 -*-

""" Methods that concern splice junctions - getting their coordinates, transcription starts/ends, phases. """

import logging
import os.path
import pandas as pd

from jcast import params

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
        The slices should be alt1-anc and alt2-anc for
        2020-07-235 A5SS may be treated differently here: I think for (-) strand, the slices are
        anc-alt1 and anc-alt2
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
        2020-07-235 A3SS may be treated differently here: I think for (-) strand, the slices are
        alt1-anc and alt2-anc
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
        self.fdr = kwargs['fdr']
        self.sjc_s1 = kwargs['sjc_s1']
        self.sjc_s2 = kwargs['sjc_s2']

        self.tx1 = None
        self.tx0 = None
        self.phase = None
        self.num_start_codons = 0

        self.logger = logging.getLogger('jcast.junction')


    def __repr__(self):
        """ repr """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    def __str__(self):
        """ str """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    @property
    def sum_sjc(self):
        """
        Returns the sum of all SJCs in the junction. if rMATS was run with one technical replicate,
        the count field is an int, otherwise it is a list. Currently this takes the skipped junction count (SJC)
        as filtering criterion because the majority of translatable events are probably SE (skipped exon).
        Essentially this filters out alternative junctions that are very rarely skipped (high inclusion
        level of the exons) that are not likely to be translatable.

        :return:
        """
        try:
            sum_count_sample1 = int(sum([int(x) for x in (str(self.sjc_s1).split(sep=','))]))
            sum_count_sample2 = int(sum([int(x) for x in (str(self.sjc_s2).split(sep=','))]))

        except ValueError:
            sum_count_sample1 = 0
            sum_count_sample2 = 0

        return sum_count_sample1 + sum_count_sample2

    def _get_translated_region(self,
                               gtf,
                               startsite_index: int = 0,
                               ):
        """
        Read the genome annotation (.gtf) file, find the coding sequences (CDS) that share the gene name and coordinates
        of the anchor exon (anc) supplied, then find out where the annotated translation start and end sites are. If the
        splice junction sequences extend BEYOND the translation starts and ends, trim them to avoid running into stop
        codons.

        :param gtf: genome annotation
        :param startsite_index: the index of which start site to retrieve; default is 0 (most upstream TSS).
        :return: True
        """
        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')

        #
        # 	Get the translation start and end positions
        #
        tsl = params.tsl_threshold

        # 2020-07-25 now getting the start codons of all protein coding transcripts at TSL threshold
        gtf0_start = gtf0.query('feature == "start_codon" & '
                                'transcript_biotype == "protein_coding" & '
                                        'transcript_support_level <= @tsl').loc[:, 'start'].drop_duplicates()

        # Number of start sites:
        self.num_start_codons = len(gtf0_start)

        # If there are retrievable start site:
        if self.num_start_codons > 0:

            # Get the start site for the longest transcript (lowest coordinates if strand is +)
            gtf0_start = sorted(gtf0_start, reverse=self.strand == '-')

            self.tx0 = gtf0_start[startsite_index]

        # 2020-07-25 now getting end codons of all protein coding transcripts at TSL threshold
        try:
            gtf0_end = gtf0.query('feature == "start_codon" & '
                                  'transcript_biotype == "protein_coding" & '
                                  'transcript_support_level <= @tsl').loc[:, 'start'].drop_duplicates()

        except KeyError:
            gtf0_end = None

        # Get the longest transcripts (highest coordinates if strand is +)
        if len(gtf0_end) > 0:
            self.tx1 = sorted(gtf0_end, reverse=self.strand == '+')[0]
        else:
            self.tx1 = -1

        #
        # by rMATS convention, if strand is -ve
        # then the upstream is near the end of tx
        # shift by 2 to get to the end of the end codon.
        #
        if self.strand == '-' and self.tx1 > 0 and self.tx0 > 0:
            self.tx0, self.tx1 = self.tx1, self.tx0
            self.tx1 += 2

        elif self.strand == '+':
            self.tx1 += 2

        self.logger.debug('Chosen start codon is {0}; end codon is {1}; tsl is {2}.'.format(self.tx0,
                                                                                            self.tx1,
                                                                                            tsl,
                                                                                            )
                          )
        return True

    def get_translated_phase(self, gtf):
        """
        Get the annotated translation phase from the GTF file
        :param gtf:  genome annotation
        :return:
        """

        """
        Get translation phase from GTF file.
        If there is no phase found in the GTF, use phase -1 for now.
        # TODO: look more closely into GTF file, or try translating from all frames
        """

        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')

        # Select the anchor exon from CDS and get the frame
        # Note 2018-03-24 this is probably not quite right. I think you should find out whether the anchor
        # is really the one that determines the phase here.

        ph0, ph1 = None, None

        # 2018-03-24: First define the exon we are looking for.
        if self.junction_type in ['MXE', 'SE', 'RI']:
            if self.strand == '+':
                ph0, ph1 = self.anc_es, self.anc_ee
            elif self.strand == '-':
                ph0, ph1 = self.down_es, self.down_ee

        elif self.junction_type == 'A5SS':
            if self.strand == '+':
                ph0, ph1 = self.alt1_es, self.alt1_ee  # Might have to search also for alt2
            elif self.strand == '-':
                #
                # 2020-07-25 changed to using the alt1 to retrieve phase because it is also upstream
                # ph0, ph1 = self.anc_es, self.anc_ee #
                ph0, ph1 = self.alt1_es, self.alt1_ee

        elif self.junction_type == 'A3SS':
            if self.strand == '+':
                ph0, ph1 = self.anc_es, self.anc_ee
            elif self.strand == '-':
                #
                # 2020-07-25 changed to using the anc to retrieve phase because it is also upstream
                # ph0, ph1 = self.alt1_es, self.alt1_ee
                ph0, ph1 = self.anc_es, self.anc_ee

        self.logger.info('Anchor exon start {0} Anchor exon end {1}'.format(ph0,
                                                                            ph1))

        # Get the frame of that coding exon from GTF.
        coding_exon = gtf0.query('start == @ph0').\
            query('end == @ph1').query('feature == "CDS" & transcript_biotype == "protein_coding"')
            # 2020-07-25 added protein coding filter in case a nonsense-mediated decay CDS comes first


        # If phases retrieved, get the first value
        if len(coding_exon) > 0:
            self.phase = int([x for x in coding_exon.loc[:, 'frame'].iloc if x != '.'][0])
            # TODO: find the canonical transcript

        else:
            self.phase = None     # dummy value to trigger trying different phases for longest translation
            # TODO: handle phase retrieval failure in a tidier manner


        self.logger.info('Transcription start: {0} Transcript end: {1}'.format(self.tx0,
                                                                               self.tx1))
        self.logger.info('Retrieved phase: {0}'.format(self.phase))

    def trim_cds(self, gtf):
        """
        Wrapper to trimming the exons by CDS; first retrieve the coordinates of the translation starts and ends,
        then
        """

        if self.tx0 is None:
            self._get_translated_region(gtf=gtf,
                                        startsite_index=0)

        self._trim()

        # 2020-07-25 Check if the newly trimmed exon is actually part of a cds
        # TODO: only doing that for MXE and SE for now and only getting second start site for now
        cds = gtf.annot.query('gene_id == @self.gene_id & '
                              'feature == "CDS" & '
                              'transcript_biotype == "protein_coding" &'
                              'start == @self.anc_es & '
                              'end == @self.anc_ee')

        # 2020-07-25 If the newly trimmed anchor is not part of a CDS exon, this may suggest there is a second start site
        if len(cds) == 0 and self.num_start_codons > 1:
            self._get_translated_region(gtf=gtf,
                                        startsite_index=1)
            self._trim()


    def _trim(self):
        """
        Core trim function that shaves off the UTRs
        :return: True

        Trims the junction based on transcription start and end:

        """

        #
        # Subset the gtf file by the current gene_id
        #

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

            except TypeError or ValueError:
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

            except TypeError or ValueError:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'SE':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.anc_es = -1
                    self.anc_ee = -1
                    self.alt1_es = self.tx0

            except TypeError or ValueError:
                self.logger.info('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.down_es = -1
                    self.down_ee = -1
                    self.alt1_ee = self.tx1

            except TypeError or ValueError:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'RI':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except TypeError or ValueError:
                self.logger.info('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

            except TypeError or ValueError:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'A5SS':
            try:
                if self.alt2_ee > self.tx0 > self.alt2_es:
                    self.alt2_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.alt1_es = self.tx0

            except TypeError or ValueError:
                self.logger.info('Trimming start failed.')

            try:
                if self.anc_ee > self.tx1 > self.anc_es:
                    self.anc_ee = self.tx1

            except TypeError or ValueError:
                self.logger.info('Trimming end failed.')

        elif self.junction_type == 'A3SS':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except TypeError or ValueError:
                self.logger.info('Trimming start failed.')

            try:
                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.alt1_ee = self.tx1

                if self.alt2_ee > self.tx1 > self.alt2_es:
                    self.alt2_ee = self.tx1

            except TypeError or ValueError:
                self.logger.info('Trimming end failed.')

        return True
