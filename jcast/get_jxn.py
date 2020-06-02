# -*- coding: utf-8 -*-

""" Classes that concern junctions - getting their coordinates, transcription starts/ends, phases, and annotations. """


import os.path
import pandas as pd
import gtfparse as gtp

# import logging

class Annotation(object):
    """
    Class holds the name and location of the GTF file, plus a pandas dataframe
    """

    def __init__(self, path):
        """

        :param path: Path of the GTF file
        """
        self.path = os.path.join(path)
        self.annot = 0

    def read_gtf(self):
        """
        Read gtf file based on the location and name supplied
        :return:
        """

        ind = self.path[:-4] + 'indexed.txt'

        if os.path.isfile(ind) is False:
            self.annot = gtp.read_gtf(self.path)
            self.annot.to_csv(ind, encoding='utf-8', sep='\t')

        self.annot = pd.read_table(ind, sep='\t', low_memory=False)

        return True

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
        self.junction_type = kwargs['junction_type']
        #self.species = kwargs['species']
        self.gene_symbol = kwargs['gene_symbol']
        self.tx1 = -1
        self.tx0 = -1
        self.phase = -1
        self.min_read_count = 0

    def __str__(self):
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

        :param gtf: The annotation file
        :return: True
        """
        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')



        #
        # 	Get the translation start and end positions
        #
        gtf0_start = gtf0.query('feature == "start_codon"').sort_values(['transcript_support_level']).sort_values(['ccds_id']).loc[:, 'start']
        if len(gtf0_start) > 0:
            self.tx0 = gtf0_start.iloc[0]
        else:
            self.tx0 = -1

        gtf0_end = gtf0.query('feature == "stop_codon"').sort_values(['transcript_support_level']).loc[:, 'start']
        if len(gtf0_end) > 0:
            self.tx1 = gtf0_end.iloc[0]
        else:
            self.tx1 = -1

        #
        # By rMATS convention, if strand is -ve
        # then the upstream is near the end of tx
        #
        if self.strand == '-' and self.tx1 > 0 and self.tx0 > 0:
            self.tx1 += 2
            self.tx0 += 2
            temp = self.tx1
            self.tx1 = self.tx0
            self.tx0 = temp


        return True

    def get_translated_phase(self, gtf):
        """
        Get the annotated translation phase from the GTF file
        :param gtf:
        :return:
        """

        # Select the anchor exon from CDS and get the frame
        # Note 2018-03-24 this is probably not quite right. I think you should find out whether the anchor
        # is really the one that determines the phase here.

        # 2018-03-24: First define the exon we are looking for.
        if self.junction_type in ['MXE', 'SE', 'RI']:
            if self.strand == '+':
                ph0 = self.anc_es
                ph1 = self.anc_ee
            if self.strand == '-':
                ph0 = self.down_es
                ph1 = self.down_ee

        elif self.junction_type == 'A5SS':
            if self.strand == '+':
                ph0 = self.alt1_ee
                ph1 = self.alt1_es # Might have to search also for alt2
            if self.strand == '-':
                ph0 = self.anc_es
                ph1 = self.anc_ee

        elif self.junction_type == 'A3SS':
            if self.strand == '+':
                ph0 = self.anc_ee
                ph1 = self.anc_es  # Might have to search also for alt2
            if self.strand == '-':
                ph0 = self.alt1_es
                ph1 = self.alt2_ee

        # Get the frame of that coding exon from GTF.
        gtf0 = gtf.annot.query('gene_id == @self.gene_id').query('start == @ph0').\
            query('end == @ph1').query('feature == "CDS"')

        if len(gtf0) > 0:

            self.phase = gtf0.loc[:, 'frame'].iloc[0]
            # print(self.phase)
            if self.phase != '.':
                self.phase = int(self.phase)
                #print('Retrieved phase: ' + str(self.phase))

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
                print('Trimming start failed.')

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
                print('Trimming end failed.')

        elif self.junction_type == 'SE':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.anc_es = -1
                    self.anc_ee = -1
                    self.alt1_es = self.tx0

            except:
                print('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.down_es = -1
                    self.down_ee = -1
                    self.alt1_ee = self.tx1

            except:
                print('Trimming end failed.')

        elif self.junction_type == 'RI':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except:
                print('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

            except:
                print('Trimming end failed.')


        elif self.junction_type == 'A5SS':
            try:
                if self.alt2_ee > self.tx0 > self.alt2_es:
                    self.alt2_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.alt1_es = self.tx0

            except:
                print('Trimming start failed.')

            try:
                if self.anc_ee > self.tx1 > self.anc_es:
                    self.anc_ee = self.tx1

            except:
                print('Trimming end failed.')


        elif self.junction_type == 'A3SS':
            try:
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

            except:
                print('Trimming start failed.')

            try:
                if self.alt1_ee > self.tx1 > self.alt1_es:
                    self.alt1_ee = self.tx1

                if self.alt2_ee > self.tx1 > self.alt2_es:
                    self.alt2_ee = self.tx1

            except:
                print('Trimming end failed.')

        return True

    # def write_fate(self, fate, output):
    #     """
    #     Write out the outcome of the attempt to translate each junction into a report file
    #
    #     :param fate:    int         Code for message to be writtebn
    #     :param output:  string      Output directory
    #     :return:
    #     """
    #
    #     import os.path
    #
    #     os.makedirs('out', exist_ok=True)
    #     o = os.path.join('out', output + '_' + 'fate' + '.txt')
    #     print(o)
    #
    #     # Set the stored junction fate as the message
    #     self.fate = fate
    #
    #     assert type(self.fate) is int, 'Junction fate code error.'
    #
    #     if self.fate == -2:
    #         msg = 'DELETED. Junction read counts too low.'
    #
    #     elif self.fate == -1:
    #         msg = 'DELETED. Junction inconsistent across replicates.'
    #
    #     elif self.fate == 0:
    #         msg = ''
    #
    #     elif self.fate == 1:
    #         msg = "SUCCESS 1. Retrieved phase: " + str(
    #                     self.phase) + " Used phase: " + str(self.translated_phase) + ". No Frameshift."
    #
    #     elif self.fate == 2:
    #         msg = "SUCCESS 2. Retrieved phase: " + str(
    #                     self.phase) + " Used phase: " + str(self.translated_phase) + ". Frameshift."
    #
    #     elif self.fate == 3:
    #         msg = "SUCCESS 3. The GTF frame appears to be wrong. Retrieved phase: " + str(
    #         self.phase) + " Used phase: " + str(self.translated_phase)
    #
    #     elif self.fate == 4:
    #         msg = "WARNING 4. Slice 2 hit a stop codon. Used longest phase."
    #
    #     elif self.fate == 5:
    #         msg = "WARNING 5. Slice 1 hit a stop codon. Used longest phase."
    #
    #     elif self.fate == 6:
    #         msg = 'FAILURE. No translation was done. At least one PTC at each frame.'
    #
    #     elif self.fate == 7:
    #         msg = 'SUCCESS. Six frame translation done.'
    #
    #     else:
    #         raise AssertionError
    #
    #     f = open(o, 'a')
    #     f.write(self.junction_type + '\t' + self.name + '\t' + self.gene_symbol + '\t' + msg + '\n')
    #     f.close()
    #
    #     return True

#
#   For doctest
#
