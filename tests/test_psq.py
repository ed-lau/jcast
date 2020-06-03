# -*- coding: utf-8 -*-

""" Tests """


import unittest
import tempfile
import sys
import os

from io import StringIO

from Bio import SeqIO
from Bio.Alphabet import IUPAC

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util import Retry

from jcast.get_jxn import Junction
from jcast.get_gtf import ReadAnnotations
from jcast.get_seq import Sequence
from jcast.get_rma import RmatsResults
from jcast.read_fa import ReadGenome

class GenomeTest(unittest.TestCase):
    """
    Test cases involving reading genome sequence and references
    """

    def setUp(self):
        """

        :return:
        """

        global gtf, rma, genome, test_data_loc

        test_data_loc = os.path.join('tests', 'data')
        genome_loc = os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
        gtf_loc = os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf')
        rmats_loc = os.path.join(test_data_loc, 'rmats')

        genome = ReadGenome(genome_loc)

        gtf = ReadAnnotations(gtf_loc)
        gtf.read_gtf()

        rmats_results = RmatsResults(rmats_loc)
        rma = rmats_results.__getattribute__('rmats_mxe')

        pass

    def tearDown(self):

        # Remove the cached gtf file
        os.remove(os.path.join(test_data_loc, 'genome', 'Homo_sapiens.GRCh38.89.chromosome.15.gtf.cached'))

        pass

    def test_that_genome_loads(self):

        i = 71 # PKM

        junction = Junction(id=rma.id[i], \
                            gene_id=rma.gene_id[i], \
                            strand=rma.strand[i], \
                            gene_symbol=rma.gene_symbol[i], \
                            chr=rma.chr[i], \
                            anc_es=rma.anc_es[i], \
                            anc_ee=rma.anc_ee[i], \
                            alt1_es=rma.alt1_es[i], \
                            alt1_ee=rma.alt1_ee[i], \
                            alt2_es=rma.alt2_es[i], \
                            alt2_ee=rma.alt2_ee[i], \
                            down_es=rma.down_es[i], \
                            down_ee=rma.down_ee[i], \
                            junction_type=rma.jxn_type[i], \
                            species='human', )

        junction.get_translated_region(gtf)
        junction.get_translated_phase(gtf)
        junction.trim()


        sequence = Sequence(junction, directory_to_write=tempfile.tempdir)
        sequence.make_slice_localgenome(genome.genome)
        sequence.translate(use_phase=True)


        server = 'https://www.ebi.ac.uk'
        ext = '/proteins/api/proteins/Ensembl:' + sequence.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

        # retry 10 times
        retries = Retry(total=10,
                        backoff_factor=0.1,
                        status_forcelist=[500, 502, 503, 504])

        rqs = requests.Session()
        rqs.mount('https://', HTTPAdapter(max_retries=retries))
        ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

        if not ret.ok:
            print("Network still not okay after 10 retries. Quitting.")
            ret.raise_for_status()
            sys.exit()

        fasta_handle = SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein)

        merge_length = 10
        for loop in fasta_handle:
            record = loop[:]  # [:] needed to copy list rather than add new alias

        # Find out where the first (10) amino acids meets the UniProt canonical sequences..
        merge_start1 = record.seq.find(sequence.slice1_aa[:merge_length])
        merge_end1 = record.seq.find(sequence.slice1_aa[-merge_length:])

        merge_start2 = record.seq.find(sequence.slice2_aa[:merge_length])
        merge_end2 = record.seq.find(sequence.slice2_aa[-merge_length:])

        if merge_start1 != -1 and merge_end1 != -1:
            record1 = record[:merge_start1] + sequence.slice1_aa + record[merge_end1 + merge_length:]
        elif merge_start1 != -1 and merge_end1 == -1:
            record1 = record[:merge_start1] + sequence.slice1_aa
        elif merge_start1 == -1 and merge_end1 != -1:
            record1 = record[:merge_start1] + sequence.slice1_aa + record[merge_end1 + merge_length:]
        elif merge_start1 == -1 and merge_end1 == -1:
            record1 = sequence.slice1_aa

        if merge_start2 != -1 and merge_end2 != -1:
            record2 = record[:merge_start2] + sequence.slice2_aa + record[merge_end2 + merge_length:]
        elif merge_start2 != -1 and merge_end2 == -1:
            record2 = record[:merge_start2] + sequence.slice2_aa
        elif merge_start2 == -1 and merge_end2 != -1:
            record2 = record[:merge_start2] + sequence.slice2_aa + record[merge_end2 + merge_length:]
        elif merge_start2 == -1 and merge_end2 == -1:
            record2 = sequence.slice2_aa

        print(record1.seq)
        print(record2.seq)

        return True