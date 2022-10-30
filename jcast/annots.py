# -*- coding: utf-8 -*-

""" Genome annotations. """
import _io
import argparse
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

    def __init__(self,
                 logger: logging.Logger,
                 path: _io.TextIOWrapper, ):
        """
        :parrm logger: logger object
        :param path: Path to GTF file
        """
        self.logger = logger
        self.path = path

        self.annot = pd.DataFrame()

        self.read_gtf()


    def read_gtf(self) -> None:
        """
        Read gtf file based on the location and name supplied
        :return:
        """

        self.logger.info("Reading GTF file. This could take a minute.")

        # Writes cached gtf if it doesn't exist
        cached_gtf = self.path.name + '.cached'

        if os.path.isfile(cached_gtf) is False:
            self.annot = gtfparse.read_gtf(self.path)
            self.annot.to_csv(cached_gtf, encoding='utf-8', sep='\t')

        self.annot = pd.read_table(cached_gtf, sep='\t', low_memory=False)

        # 2020-11-06 if Gencode GTF is used, column is transcript_type rather than transcript_biotype.
        # This will make the transcript_biotype column manually

        if 'transcript_biotype' not in self.annot.columns:
            self.annot['transcript_biotype'] = self.annot['transcript_type']

        return None


class ReadGenome(object):

    def __init__(self,
                 logger: logging.Logger,
                 path: _io.TextIOWrapper, ):
        """
        :param logger: logger object
        :param path: Path to genome fasta file

        """


        self.path = path
        self.logger = logger

        self.genome = None
        self._read_fasta()


    def _read_fasta(self):
        """
        Reading genome. Just a wrapper for Biopython SeqIO with gzip if needed.

        :return:
        """

        if self.path.name.endswith('.gz'):
            with gzip.open(self.path, 'rt') as f:
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                self.logger.info('Read from zipped genome. \n\n')

        else:
            with open(self.path.name, 'rt') as f:
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
