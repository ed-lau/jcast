# -*- coding: utf-8 -*-

""" Methods that concern genome annotations. """

import os.path
import logging
import gtfparse
import pandas as pd

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
