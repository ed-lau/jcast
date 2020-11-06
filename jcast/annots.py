# -*- coding: utf-8 -*-

""" Genome annotations. """

import os.path
import logging
import gzip

import gtfparse
import pandas as pd
from Bio import SeqIO
#from Bio.Alphabet import generic_dna


class ReadAnnotations(object):
    """
    Class holds the name and location of the GTF file, plus a pandas dataframe
    """

    def __init__(self, path):
        """

        :param path: Path of the GTF file
        """
        self.path = os.path.join(path)
        self.annot = 0

        self.logger = logging.getLogger('jcast.gtf')

    def read_gtf(self):
        """
        Read gtf file based on the location and name supplied
        :return:
        """

        self.logger.info("Reading GTF file. This could take a minute.")

        # Writes cached gtf if it doesn't exist
        cached_gtf = self.path + '.cached'

        if os.path.isfile(cached_gtf) is False:
            self.annot = gtfparse.read_gtf(self.path)
            self.annot.to_csv(cached_gtf, encoding='utf-8', sep='\t')

        self.annot = pd.read_table(cached_gtf, sep='\t', low_memory=False)

        # 2020-11-06 if Gencode GTF is used, column is transcript_type rather than transcript_biotype.
        # This will make the transcript_biotype column manually

        if not 'transcript_biotype' in self.annot.columns:
            self.annot['transcript_biotype'] = self.annot['transcript_type']

        return True



class ReadGenome(object):

    def __init__(self, f_loc):
        """

        :param f_loc: :param f_loc:   Location of fasta file (.fa or .gz)
        """

        self.f_loc = f_loc
        self.logger = logging.getLogger('jcast.genome')

        self.genome = None
        self._read_fasta()


    def _read_fasta(self):
        """
        Reading genome. Just a wrapper for Biopython SeqIO with gzip if needed.

        :return:
        """

        if self.f_loc.endswith('.gz'):
            with gzip.open(self.f_loc, 'rt') as f:
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                self.logger.info('Read from zipped genome. \n\n')

        else:
            with open(self.f_loc, 'rt') as f:
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                self.logger.info('Read from genome. \n\n')

        return True


class AnnotatedTranscript(object):
    """ Holds the transcripts read from GTF. """

    def __init__(self,
                 transcript_name,
                 protein_id,
                 exons,
                 start_codon,
                 end_codon,
                 starting_translation_phase,
                 ):
        """ init """
        self.transcript_name = transcript_name
        self.protein_id = protein_id
        self.exons = exons
        self.start_codon = start_codon
        self.end_codon = end_codon
        self.starting_translation_phase = starting_translation_phase

    def __repr__(self):
        """ repr """
        return 'Annotated transcript from GTF {0}'.format(self.transcript_name)

    def __str__(self):
        """ str """
        return '{0}'.format(self.transcript_name)

    def __len__(self):
        """ get length of all exons """
        return sum([end-start+1 for (start, end) in self.exons])
