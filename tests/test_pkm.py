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

from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence

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

    def test_genome_loads(self):

        # TODO: test negative strand and trimming
        i = 71  # PKM

        junction = Junction(**rma.iloc[i].to_dict())

        junction.trim_cds(gtf)

        sequence = Sequence(junction)
        sequence.make_slice_localgenome(genome.genome)
        sequence.translate(use_phase=True)

        server = 'https://www.ebi.ac.uk'
        ext = '/proteins/api/proteins/Ensembl:' + sequence.j.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

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

        merge_length = 10

        with StringIO(ret.text) as fas:
            parsed = SeqIO.parse(fas, 'fasta', IUPAC.extended_protein)
            for loop in parsed:
                record = loop[:]  # [:] needed to copy list rather than add new alias

        ret.close()

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

        #print(record1.seq)
        #print(record2.seq)

        # Compare with Uniprot retrieved sequences
        retrieved_pkm2 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                          'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                          'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                          'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                          'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                          'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                          'LSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEAS'
                          'FKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                          'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

        retrieved_pkm1 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                          'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                          'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                          'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                          'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                          'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                          'LSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEAS'
                          'YKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                          'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

        self.assertEqual(record1.seq, retrieved_pkm2)
        self.assertEqual(record2.seq, retrieved_pkm1)