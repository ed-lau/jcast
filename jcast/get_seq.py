# -*- coding: utf-8 -*-

""" Methods that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids """


import logging
import os.path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import requests as rq
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from io import StringIO
import sqlite3 as sq

from jcast import helpers as h

# TODO: Sequence should really inherit directly from junction
class Sequence(object):

    def __init__(self, junction, directory_to_write):
        """
        :type junction: object
        :param junction: The splice junction object

        Mostly copying properties of the junction object (already trimmed) to start a new sequence object. This
        sequence object will be used to make the nucleotide slices and translate into protein sequences.

        To do: Change this so that sequence inherits junction class directly (have to think it through)

        """
        self.anc_ee = junction.anc_ee
        self.anc_es = junction.anc_es
        self.alt1_es = junction.alt1_es
        self.alt1_ee = junction.alt1_ee
        self.alt2_es = junction.alt2_es
        self.alt2_ee = junction.alt2_ee
        self.down_es = junction.down_es
        self.down_ee = junction.down_ee
        #self.species = junction.species
        self.gene_id = junction.gene_id
        self.junction_type = junction.junction_type
        self.chr = junction.chr
        self.phase = junction.phase
        self.strand = junction.strand
        self.slice1_nt = ''
        self.slice2_nt = ''
        self.slice1_aa = ''
        self.slice2_aa = ''
        self.frameshift = False
        self.gene_symbol = junction.gene_symbol
        self.name = junction.name
        self.fate = 'Nothing done.'
        self.translated_phase = -1                  # Phase that was actually used for translation.
        self.translated_strand = junction.strand    # Strand that was actually used for translation
        self.min_read_count = junction.min_read_count

        self.logger = logging.getLogger('jcast.seq')


    def __str__(self):
        return 'Sequence object: ' + self.gene_id + ' ' + self.gene_symbol + ' ' + self.name

    # To do: redo this as @ property
    def set_frameshift_to_true(self):
        """
        A setter to mark that there is a frameshift; used to determine whether the slice should be tier 1 or tier 2

        :return:
        """
        self.frameshift = True

    def make_slice_localgenome(self, genome_index):
        """
        This gets the nucleotide sequence from the coordinates, using a local genome
        :return:
        """

        anc_nt = h.get_local_nuc(genome_index, self.chr, self.anc_es, self.anc_ee)
        alt1_nt = h.get_local_nuc(genome_index, self.chr, self.alt1_es, self.alt1_ee)
        alt2_nt = h.get_local_nuc(genome_index, self.chr, self.alt2_es, self.alt2_ee)
        down_nt = h.get_local_nuc(genome_index, self.chr, self.down_es, self.down_ee)

        self.slice1_nt = str((anc_nt + alt1_nt + down_nt).seq)
        self.slice2_nt = str((anc_nt + alt2_nt + down_nt).seq)

        self.logger.info('Retrieved nucleotide for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice1_nt))
        self.logger.info('Retrieved nucleotide for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice2_nt))

    def make_slice(self):
        """
        This gets the nucleotide sequence from the coordinates, using a web API
        :return:
        """

        anc_nt = h.get_nuc(self.species, self.chr, self.anc_es, self.anc_ee)
        alt1_nt = h.get_nuc(self.species, self.chr, self.alt1_es, self.alt1_ee)
        alt2_nt = h.get_nuc(self.species, self.chr, self.alt2_es, self.alt2_ee)
        down_nt = h.get_nuc(self.species, self.chr, self.down_es, self.down_ee)

        self.slice1_nt = anc_nt + alt1_nt + down_nt
        self.slice2_nt = anc_nt + alt2_nt + down_nt

        print(self.slice1_nt)
        print(self.slice2_nt)

    def translate(self, use_phase=True):
        """
        This is the translate function for tier 1/2/3 peptide. Calls make_pep with the terminal option
        which will terminate when it runs into a stop codon. Later on we should make a force_translate that does
        three frame translation and just either return all frames or return the longest.

        :param use_phase: T/F Whether to ue the stored phase or attempt to do three-frame translation
        :return: True
        """

        if self.phase in [0, 1, 2] and use_phase:
            self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, self.phase, terminate=True)
            self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, self.phase, terminate=True)
            self.logger.debug('Used retrieved phase for {0} {1}:'.format(self.name, self.gene_symbol))
            self.translated_phase = self.phase

        else:
            for i in range(3):
                self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, i, terminate=True)
                self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, i, terminate=True)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    self.translated_phase = i
                    break
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice1_aa))
        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice2_aa))

        return True

    def translate_forced(self, slice_to_translate):
        """
        This is the fallback translate function for tier 4/5 peptides. We will start by using the retrieved GTF phase
        and attempt to translate as normal. If only slice_1 returns normal without stop codon but slice_2 runs into
        a stop codon, then I think what we could do is to try all three phases for slice_2 and take the longest peptide.

        This is intended to catch some suspected sequences where there *is* a biological premature stop codon, or that
        the trimming of translation ends did not complete correctly. I am suspecting that this may be more frequent for
        A5SS and RI sequences, which in the first three tiers did not perform as well as MXE and SI sequences.

        :param slice_to_translate:  Int     Should be 1 or 2, depending on which slice we want to force
        :return:
        """

        assert slice_to_translate in [1, 2], 'Forced translation: slice must be either 1 or 2'

        # Determine whether slice 1 or slice 2 nucleotides is being forced-translated
        if slice_to_translate == 1:
            nt_to_translate = self.slice1_nt
        elif slice_to_translate == 2:
            nt_to_translate = self.slice2_nt

        # For each nucleotide sequence, do each of three phases then get the longest translated product
        for i in range(3):
            best_seq = ''
            best_phase = -1
            seq = h.make_pep(nt_to_translate, self.strand, i, terminate=False)

            # Not taking care of equal lengths for now - only taking first phase if there are two with equal length.
            if len(seq) > len(best_seq):
                best_seq = seq
                best_phase = i

        if len(best_seq) > 0:
            self.translated_phase = best_phase
            forced_translated_aa = best_seq

        else:
            forced_translated_aa = ''

        # Determine whether to write the force_translated AA into slice 1 or slice 2
        if slice_to_translate == 1:
            self.slice1_aa = forced_translated_aa
        elif slice_to_translate == 2:
            self.slice2_aa = forced_translated_aa

        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice1_aa))
        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.name, self.gene_symbol, self.slice2_aa))

        return True

    def translate_sixframe(self, given_strand, given_phase):
        """
        This is the most basic translate function which we will use for six-frame translation
        It will just take the phase and strand being sent to it, and it will return the peptide even if it
        runs into a stop codon.

        :param given_strand:        chr strand to use for translation
        :param given_phase:         chr phase to use for translation
        :return:                    T
        """"""
        """

        assert given_strand in ['+', '-'], 'Given strand is not a valid string'
        assert given_phase in [0, 1, 2], 'Given phase is not a valid integer'

        self.slice1_aa = h.make_pep(self.slice1_nt, given_strand, given_phase, terminate=False)
        self.slice2_aa = h.make_pep(self.slice2_nt, given_strand, given_phase, terminate=False)

        self.translated_phase = given_phase
        self.translated_strand = given_strand

        return True

    def write_to_fasta(self, output, suffix):
        """

        WARNING: Deprecated, see the new "extend and write" function below, which will first attempt to harmonize
        the sequence to FASTA before writing the full-length sequence.

        :param output: File name of the .fasta output.
        :return: True

        Create the /out directory if it does not exist, then write the translated splice junctions into the .fasta file

        """

        # Note to self, this condition should probably be moved to main to make this function more reusable.
        # So the idea is to write out different files depending on whether translation was successful

        fa1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                        id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-1-' +
                            self.name + '-' + str(self.phase) + self.strand),
                        name=self.gene_symbol,
                        description='Slice 1')
        fa2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
                        id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-2-' +
                            self.name + '-' + str(self.phase) + self.strand),
                        name=self.gene_symbol,
                        description='Slice 2')

        os.makedirs('../out', exist_ok=True)
        o = os.path.join('../out', output + suffix + '.fasta')

        output_handle = open(o, 'a')
        SeqIO.write([fa1, fa2], output_handle, 'fasta')
        output_handle.close()


        return True

    # To do this needs to be redone
    def extend_and_write(self,
                         # species,
                         output='out',
                         suffix='T0',
                         merge_length=10):
        """
        Given a translated junction sequence, look for the fasta entry that overlaps with it, then return the entry
        and the coordinates. This will be used to extend said junction sequence to encompass  entire protein sequence.

        While this does not help resolve junction sequences or isoforms, I think it is important to match the database
        file as closely as possible to the actual spectra in a mass spectrometry experiment to strengthen the
        hypotheses being tested in the database search.

        Instead of reading from the FASTA file, it should just read fetch the Uniprot directly via API first. However
        if the Uniprot API returns an empty SwissProt result (somehow the Ensembl gene ID linking to a TremBL but not
        Swissprot sequence), then we will fall back to a pre-loaded SwissProt only FASTA file.

        After extension, write the SeqRecord objects created into fasta file...

        :param species:         string  Specices (mouse or human) to determine which fasta to grab
        :param output:          string  Output directory
        :param merge_length:    int     The minimal number of amino acids needed for merging (default 10)
        :param suffix:          string  Additional suffix to add to the end of an output file
        :return:
        """


        assert type(merge_length) is int and merge_length >= 6, 'Merge length must be integer and at least 6'

        #'''
        #Retrieve sequences from Ensembl.



        # server = 'https://www.ebi.ac.uk'
        # ext = '/proteins/api/proteins/Ensembl:' + self.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'
        #
        # self.logger.debug(server + ext)
        #
        #
        # # retry 10 times
        # retries = Retry(total=10,
        #                 backoff_factor=0.1,
        #                 status_forcelist=[500, 502, 503, 504])
        #
        # rqs = rq.Session()
        # rqs.mount('https://', HTTPAdapter(max_retries=retries))
        # ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})
        #
        # if not ret.ok:
        #     self.logger.error("Network still not okay after 10 retries. Quit without writing.")
        #     return True
        #
        # # If ret.text is empty, the fall back is to use a local fasta or exit
        # if len(ret.text) == 0:
        #     return True
        #
        #     # # Load local fasta (for now) based on species
        #     # if species == 'mouse':
        #     #     fasta_handle = SeqIO.parse('data/fasta/20170918_Mm_Sp_16915.fasta', 'fasta',
        #     #                                IUPAC.extended_protein)
        #     # elif species == 'human':
        #     #     fasta_handle = SeqIO.parse('data/fasta/20170918_Hs_Sp_20205.fasta', 'fasta',
        #     #                                IUPAC.extended_protein)
        #
        #
        #
        # # If ret.text is not empty, then get from online record.
        # else:
        #     fasta_handle = SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein)


        # The UniProt API retrieves a retrieval object, with a text field inside ret.text
        # Since Biopython SeqIO only works with file, use io.StringIO to turn the string into a file for parsing.

        #'''

        '''
        # Load local fasta (for now) based on species
        if species == 'mouse':
            fasta_handle = SeqIO.parse('data/fasta/20170918_Mm_Sp_16915.fasta', 'fasta',
                                       IUPAC.extended_protein)
        elif species == 'human':
            fasta_handle = SeqIO.parse('data/fasta/20170918_Hs_Sp_20205.fasta', 'fasta',
                                       IUPAC.extended_protein)
        '''

        cache = self.gene_id

        self.logger.info(cache)

        # Create cache folder if not exists
        if not os.path.exists('cache'):
            os.makedirs('cache')

        con = sq.connect(os.path.join('cache', 'uniprot-cache.db'))
        cur = con.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')

        # try:

        cur.execute('''SELECT id, seq FROM sequences WHERE id=:cache''',
                    {'cache': cache})
        read_fasta = cur.fetchone()

        if read_fasta:
            record = list(SeqIO.parse(StringIO(read_fasta[1]), 'fasta', IUPAC.extended_protein))[0]
            self.logger.info('Locally cached sequence retrieved')
            con.close()
            # return nuc

        else:
            self.logger.info("Sequence not yet cached locally. 0")

            server = 'https://www.ebi.ac.uk'
            ext = '/proteins/api/proteins/Ensembl:' + self.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

            self.logger.info(server + ext)
            retries = Retry(total=15,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            rqs = rq.Session()
            rqs.mount('https://', HTTPAdapter(max_retries=retries))
            ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

            if not ret.ok:
                self.logger.warning("Network still not okay after 10 retries. Skipped protein.")
                return True

            if ret.status_code == 200 and ret.text !='':
                record = list(SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein))[0]

                cur.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''',
                            {'id': cache, 'seq': record.format('fasta')})
                con.commit()
                con.close()
                self.logger.info("Sequence retrieved from Uniprot and written into local cache.")

            elif ret.status_code == 200 and ret.text == '':
                self.logger.info('Retrieved empty fasta from Ensembl. Skipped protein.')
                return True

            elif ret.status_code != 200:
                self.logger.warning('Retrieval of protein sequence failed. Skipped protein.')
                return True

        # Find out where the first (10) amino acids meets the UniProt canonical sequences..
        merge_start1 = record.seq.find(self.slice1_aa[:merge_length])
        merge_end1 = record.seq.find(self.slice1_aa[-merge_length:])

        merge_start2 = record.seq.find(self.slice2_aa[:merge_length])
        merge_end2 = record.seq.find(self.slice2_aa[-merge_length:])

        # Only proceed to write file if we can bridge the first and last 10 amino acids of either
        # Slice 1 or slice 2 to the sequence.

        if (merge_start1 != -1 and merge_end1 != -1) or (merge_start2 != -1 and merge_end2 != -1):

            # Write the UniProt canonical first
            canonical = record[:]  # [:] needed to copy list rather than add new alias

            # Format the name of the slice 1 record
            if merge_start1 != -1 and merge_end1 != -1:
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + merge_length:]
            elif merge_start1 != -1 and merge_end1 == -1:
                record1 = record[:merge_start1] + self.slice1_aa
            elif merge_start1 == -1 and merge_end1 != -1:
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + merge_length:]
            elif merge_start1 == -1 and merge_end1 == -1:
                record1 = record[:0] + self.slice1_aa

            # If the slice is different from the UniProt canonical, then also write it.
            if record.seq.find(self.slice1_aa) == -1:

                record1.id += ('|' + self.gene_id + '|' + self.junction_type + '1|' + self.name + '|'
                               + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                               + '|' + self.translated_strand + str(self.translated_phase) + '|'
                               + 'r' + str(self.min_read_count) + '|' + suffix)


                h.write_seqrecord_to_fasta(record1, output, suffix)

            # If not, then change name of canonical to reflect that it is also slice 1.
            else:
                canonical.id = record1.id
                h.write_seqrecord_to_fasta(canonical, output, 'canonical')

            # Format the name of the slice 2 record
            if merge_start2 != -1 and merge_end2 != -1:
                record2 = record[:merge_start2] + self.slice2_aa + record[merge_end2 + merge_length:]
            elif merge_start2 != -1 and merge_end2 == -1:
                record2 = record[:merge_start2] + self.slice2_aa
            elif merge_start2 == -1 and merge_end2 != -1:
                record2 = record[:merge_start2] + self.slice2_aa + record[merge_end2 + merge_length:]
            elif merge_start2 == -1 and merge_end2 == -1:
                record2 = record[:0] + self.slice2_aa

                # If the slice is not the same as the UniProt canonical, then also write it.
            if record.seq.find(self.slice2_aa) == -1:
                record2.id += ('|' + self.gene_id + '|' + self.junction_type + '2|' + self.name + '|'
                               + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                               + '|' + self.translated_strand + str(self.translated_phase) + '|'
                               + 'r' + str(self.min_read_count) + '|' + suffix)

                h.write_seqrecord_to_fasta(record2, output, suffix)

            # If not, then change name of canonical to reflect that it is also slice 2.
            else:
                canonical.id = record2.id
                h.write_seqrecord_to_fasta(canonical, output, suffix)

            #print(record.id)
            #print(record.seq[:merge_start1] + self.slice1_aa + record.seq[merge_end1 + merge_length:])
            #print(record.seq[:merge_start2] + self.slice2_aa + record.seq[merge_end2 + merge_length:])

            # Once you found a match and wrote the sequence, quit.
            return True

        # If the slice is not matched to any of the FASTA entries, or if the slices are too short,
        # write the slices to an orphan fasta.
        # NOTE: we are separating these out for now because we want to find out why they fall through.

        self.logger.info('Slices are not stitchable to fasta. Writing to orphan file for reference.')
        self.logger.info('Slice 1: {0}'.format(self.slice1_aa))
        self.logger.info('Slice 2: {0}'.format(self.slice2_aa))
        self.logger.info('Canonical: {0}'.format(record.seq))

        # Format the name of the orphan slice 1 record
        orphan_slice1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                                  id=('xx|ORPHN|' + self.gene_symbol + '|'
                                      + self.gene_id + '|' + self.junction_type + '1|' + self.name + '|'
                                      + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                                      + '|' + self.translated_strand + str(self.translated_phase) + '|'
                                      + 'r' + str(self.min_read_count) + '|' + suffix),
                                  name='Protein name here',
                                  description='Description',)

        # Format the name of the orphan slice 2 record
        orphan_slice2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
                                  id=('xx|ORPHN|' + self.gene_symbol + '|'
                                      + self.gene_id + '|' + self.junction_type + '2|' + self.name + '|'
                                      + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                                      + '|' + self.translated_strand + str(self.translated_phase) + '|'
                                      + 'r' + str(self.min_read_count) + '|' + suffix),
                                  name='Protein name here',
                                  description='Description',)

        h.write_seqrecord_to_fasta(orphan_slice1, output, (suffix + '_orphan'))
        h.write_seqrecord_to_fasta(orphan_slice2, output, (suffix + '_orphan'))


        return True

