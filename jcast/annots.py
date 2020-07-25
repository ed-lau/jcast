# -*- coding: utf-8 -*-

""" Methods that concern genome annotations. """

import os.path
import logging
import gzip

import gtfparse
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna



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
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta', generic_dna))
                self.logger.info('Read from zipped genome.')

        else:
            with open(self.f_loc, 'rt') as f:
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta', generic_dna))
                self.logger.info('Read from genome.')

        return True


class AnnotatedTranscript(object):

    def __init__(self, name):
        """ init """
        self.transcript_name = name
        self.protein_id = None
        self.exons = {}
        self.start_codon = None
        self.end_codon = None
        self.is_canonical = False

