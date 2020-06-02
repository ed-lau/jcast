# -*- coding: utf-8 -*-

""" Methods to read genome fa file. """

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import logging
import gzip

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
            f = gzip.open(self.f_loc, 'rt')
            self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta', generic_dna))

            self.logger.info('Read from zipped genome')
            f.close()

        else:
            self.genome = SeqIO.to_dict(SeqIO.parse(self.f_loc, 'fasta', generic_dna))

            self.logger.info('Read from genome.')

        return True



