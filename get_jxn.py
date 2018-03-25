#
# 	Classes that concern junctions - getting their coordinates, transcription starts/ends, phases, and annotations.
#

import os.path
import pandas as pd
import gtfparse as gp

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
            self.annot = gp.read_gtf_as_dataframe(self.path)
            self.annot.to_csv(ind, encoding='utf-8', sep='\t')

        self.annot = pd.read_table(ind, sep='\t')

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
        self.species = kwargs['species']
        self.gene_symbol = kwargs['gene_symbol']
        self.tx1 = -1
        self.tx0 = -1
        self.phase = -1

    def __str__(self):
        return 'Splice junction object: ' + self.gene_id + ' ' + self.gene_symbol + ' ' + self.name

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

        print('Anchor exon start: ' + str(self.anc_es) + ' Anchor exon end: ' + str(self.anc_ee))

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

        print('Transcription start: ' + str(self.tx0) + ' Transcription end:' + str(self.tx1))

        return True

    def get_translated_phase(self, gtf):
        """
        Get the annotated translation phase from the GTF file
        :param gtf:
        :return:
        """

        # Select the anchor exon from CDS and get the frame
        gtf0 = gtf.annot.query('gene_id == @self.gene_id').query('start == @self.anc_es').query('end == @self.anc_ee').query('feature == "CDS"')

        if len(gtf0) > 0:

            self.phase = gtf0.loc[:, 'frame'].iloc[0]
            print(self.phase)
            if self.phase != '.':
                self.phase = int(self.phase)
                print('Retrieved phase: ' + str(self.phase))

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
                if self.anc_ee > self.tx0 > self.anc_es:
                    self.anc_es = self.tx0

                if self.alt1_ee > self.tx0 > self.alt1_es:
                    self.alt1_es = self.tx0

            except:
                print('Trimming start failed.')

            try:
                if self.down_ee > self.tx1 > self.down_es:
                    self.down_ee = self.tx1

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

#
#   For doctest
#
if __name__ == '__main__':
    import doctest
    doctest.testmod()
