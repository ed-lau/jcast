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
from jcast.junctions import Junction
from jcast import params


class Sequence(object):

    def __init__(self,
                 junction: Junction,
                 ):
        """
        :type junction: object
        :param junction: the splice junction object

        Mostly copying properties of the junction object (already trimmed) to start a new sequence object. This
        sequence object will be used to make the nucleotide slices and translate into protein sequences.

        """

        self.j = junction
        self.slice1_nt = ''
        self.slice2_nt = ''
        self.slice1_aa = ''
        self.slice2_aa = ''

        self.translated_phase = -1                  # Phase that was actually used for translation.
        self.translated_strand = self.j.strand    # Strand that was actually used for translation

        self.logger = logging.getLogger('jcast.seq')

    def __repr__(self):
        """ repr """
        return 'Sequence object: ' + self.j.gene_id + ' ' + self.j.gene_symbol + ' ' + self.j.name

    def __str__(self):
        """ str """
        return 'Sequence object: ' + self.j.gene_id + ' ' + self.j.gene_symbol + ' ' + self.j.name

    @property
    def frameshift(self) -> bool:
        """
        check for frameshift; used to determine whether the slice should be tier 1 or tier 2
        by checking if slice 1 nucleotides are different in length from slice 2 nucleotide by
        multiples of 3, which probably denotes frame shift (unless there are loose amino acids near the end)?

        :return:
        """

        #
        return ((len(self.slice1_nt) - len(self.slice2_nt)) % 3) != 0


    def make_slice_localgenome(self,
                               genome_index,
                               ):
        """
        This gets the nucleotide sequence from the coordinates, using a local genome

        :param genome_index: read genome file
        :return:
        """

        anc_nt = h.get_local_nuc(genome_index, self.j.chr, self.j.anc_es, self.j.anc_ee)
        alt1_nt = h.get_local_nuc(genome_index, self.j.chr, self.j.alt1_es, self.j.alt1_ee)
        alt2_nt = h.get_local_nuc(genome_index, self.j.chr, self.j.alt2_es, self.j.alt2_ee)
        down_nt = h.get_local_nuc(genome_index, self.j.chr, self.j.down_es, self.j.down_ee)

        if self.j.junction_type in ['MXE', 'SE', 'RI']:
            self.slice1_nt = str((anc_nt + alt1_nt + down_nt).seq)
            self.slice2_nt = str((anc_nt + alt2_nt + down_nt).seq)

        # 2020-07-25 the flanking exons for A5SS and A3SS:
        elif self.j.junction_type == 'A5SS':
            if self.j.strand == '+':
                self.slice1_nt = str((alt1_nt + anc_nt).seq)
                self.slice2_nt = str((alt2_nt + anc_nt).seq)
            elif self.j.strand == '-':
                self.slice1_nt = str((anc_nt + alt1_nt).seq)
                self.slice2_nt = str((anc_nt + alt2_nt).seq)

        elif self.j.junction_type == 'A3SS':
            if self.j.strand == '+':
                self.slice1_nt = str((anc_nt + alt1_nt).seq)
                self.slice2_nt = str((anc_nt + alt2_nt).seq)
            elif self.j.strand == '-':
                self.slice1_nt = str((alt1_nt + anc_nt).seq)
                self.slice2_nt = str((alt2_nt + anc_nt).seq)

        self.logger.info('Retrieved nucleotide for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice1_nt))
        self.logger.info('Retrieved nucleotide for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice2_nt))

    def translate(self,
                  use_phase=True,
                  ):
        """
        This is the translate function for tier 1/2/3 peptide. Calls make_pep with the terminal option
        which will terminate when it runs into a stop codon. Later on we should make a force_translate that does
        three frame translation and just either return all frames or return the longest.

        :param use_phase: T/F Whether to ue the stored phase or attempt to do three-frame translation
        :return: True
        """

        if self.j.phase in [0, 1, 2] and use_phase:
            self.slice1_aa = h.make_pep(self.slice1_nt, self.j.strand, self.j.phase, terminate=True)
            self.slice2_aa = h.make_pep(self.slice2_nt, self.j.strand, self.j.phase, terminate=True)
            self.logger.debug('Used retrieved phase for {0} {1}:'.format(self.j.name, self.j.gene_symbol))
            self.translated_phase = self.j.phase

        else:
            for i in range(3):
                self.slice1_aa = h.make_pep(self.slice1_nt, self.j.strand, i, terminate=True)
                self.slice2_aa = h.make_pep(self.slice2_nt, self.j.strand, i, terminate=True)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    self.translated_phase = i
                    break
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice1_aa))
        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice2_aa))

        return True

    def translate_forced(self,
                         slice_to_translate: int,
                         ):
        """
        This is the fallback translate function for tier 4 peptides. We will start by using the retrieved GTF phase
        and attempt to translate as normal. If only slice_1 returns normal without stop codon but slice_2 runs into
        a stop codon, then force-translate slice_2 and take the peptide if it reaches a certain length.

        This is intended to catch some suspected sequences where there is a biological premature stop codon, or that
        the trimming of translation ends did not complete correctly. This may be more frequent for
        A5SS and RI sequences, which in the first three tiers did not perform as well as MXE and SI sequences.

        :param slice_to_translate:  int     Should be 1 or 2, depending on which slice we want to force
        :return:
        """

        assert slice_to_translate in [1, 2], 'Forced translation: slice must be either 1 or 2'

        # Determine whether slice 1 or slice 2 nucleotides is being forced-translated
        if slice_to_translate == 1:
            nt_to_translate = self.slice1_nt
        elif slice_to_translate == 2:
            nt_to_translate = self.slice2_nt
        else:
            nt_to_translate = None

        # translate without terminating at stop codon
        seq = h.make_pep(nt_to_translate, self.j.strand, self.j.phase, terminate=False)

        if len(seq) > 0:
            self.translated_phase = self.j.phase
            force_translated_aa = seq

        else:
            force_translated_aa = ''

        # Determine whether to write the force_translated AA into slice 1 or slice 2
        if slice_to_translate == 1:
            self.slice1_aa = force_translated_aa
        elif slice_to_translate == 2:
            self.slice2_aa = force_translated_aa

        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice1_aa))
        self.logger.info('Translated AA for {0} {1}: {2}'.format(self.j.name, self.j.gene_symbol, self.slice2_aa))

        return True

    def translate_threeframe(self,
                             given_strand,
                             given_phase):
        """
        This is the most basic translate function which we will use for three-frame translation
        It will just take the phase and strand being sent to it, and it will return the peptide even if it
        runs into a stop codon.

        :param given_strand:        chr strand to use for translation
        :param given_phase:         chr phase to use for translation
        :return:
        """

        assert given_strand in ['+', '-'], 'Given strand is not a valid string'
        assert given_phase in [0, 1, 2], 'Given phase is not a valid integer'

        self.slice1_aa = h.make_pep(self.slice1_nt, given_strand, given_phase, terminate=False)
        self.slice2_aa = h.make_pep(self.slice2_nt, given_strand, given_phase, terminate=False)

        self.translated_phase = given_phase
        self.translated_strand = given_strand

        return True

    # TODO: save all sequences and write once
    def extend_and_write(self,
                         output='out',
                         suffix='T0',
                         canonical_only: bool = False,
                         ):
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

        :param output:          string  output directory
        :param suffix:          string  additional suffix to add to the end of an output file
        :param canonical_only:  bool    writes only the canonical sequence
        :return:
        """

        merge_length = params.stitch_length
        assert type(merge_length) is int and merge_length >= 6, 'Merge length must be integer and at least 6'

        # cache retrieved sequences if available, otherwise retrieve from Ensembl
        cache = self.j.gene_id
        self.logger.info(cache)

        # Create cache folder if not exists
        if not os.path.exists('cache'):
            os.makedirs('cache')

        con = sq.connect(os.path.join('cache', 'uniprot-cache.db'))
        cur = con.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')

        cur.execute('''SELECT id, seq FROM sequences WHERE id=:cache''',
                    {'cache': cache})
        read_fasta = cur.fetchone()

        if read_fasta:
            record = list(SeqIO.parse(StringIO(read_fasta[1]), 'fasta', IUPAC.extended_protein))[0]
            self.logger.info('Locally cached sequence retrieved')
            con.close()

        else:
            self.logger.info("Sequence not yet cached locally. 0")

            server = 'https://www.ebi.ac.uk'
            ext = '/proteins/api/proteins/Ensembl:' + self.j.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

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

        # find out where the first (stitch_length) amino acids meets the UniProt canonical sequences..
        merge_start1 = record.seq.find(self.slice1_aa[:merge_length])
        merge_end1 = record.seq.find(self.slice1_aa[-merge_length:])

        merge_start2 = record.seq.find(self.slice2_aa[:merge_length])
        merge_end2 = record.seq.find(self.slice2_aa[-merge_length:])

        # if meant to salvage canonical only, write the record to the gene_canonical file and exit.
        if canonical_only:
            canonical = record[:]
            h.write_seqrecord_to_fasta(canonical, output, 'gene_canonical')
            return True

        # only proceed to write file if we can bridge the first and last (stitch_length) amino acids of either
        # slice 1 or slice 2 to the sequence.
        if (merge_start1 != -1 and merge_end1 != -1) or (merge_start2 != -1 and merge_end2 != -1):

            # write the UniProt canonical first
            canonical = record[:]  # [:] needed to copy list rather than add new alias

            # format the name of the slice 1 record
            if merge_start1 != -1 and merge_end1 != -1:
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + merge_length:]
            elif merge_start1 != -1 and merge_end1 == -1:
                record1 = record[:merge_start1] + self.slice1_aa
            elif merge_start1 == -1 and merge_end1 != -1:
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + merge_length:]
            elif merge_start1 == -1 and merge_end1 == -1:
                record1 = record[:0] + self.slice1_aa

            # if the slice is different from the UniProt canonical, then also write it.
            # TODO: use format for this, unless redo canonical lookup from gtf
            if record.seq.find(self.slice1_aa) == -1:

                record1.id += ('|' + self.j.gene_id + '|' + self.j.junction_type + '1|' + self.j.name + '|'
                               + str(self.j.chr) + '|' + str(self.j.anc_ee) + '|' + str(self.j.alt1_ee)
                               + '|' + self.translated_strand + str(self.translated_phase) + '|'
                               + 'r' + str(self.j.min_read_count) + '|' + suffix)

                h.write_seqrecord_to_fasta(record1, output, suffix)

            # otherwise change name of canonical to reflect that it is also slice 1.
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

            # if the slice is not the same as the UniProt canonical, then also write it.
            if record.seq.find(self.slice2_aa) == -1:
                record2.id += ('|' + self.j.gene_id + '|' + self.j.junction_type + '2|' + self.j.name + '|'
                               + str(self.j.chr) + '|' + str(self.j.anc_ee) + '|' + str(self.j.alt1_ee)
                               + '|' + self.translated_strand + str(self.translated_phase) + '|'
                               + 'r' + str(self.j.min_read_count) + '|' + suffix)

                h.write_seqrecord_to_fasta(record2, output, suffix)

            # otherwise change name of canonical to reflect that it is also slice 2.
            else:
                canonical.id = record2.id
                h.write_seqrecord_to_fasta(canonical, output, 'canonical')

            #print(record.id)
            #print(record.seq[:merge_start1] + self.slice1_aa + record.seq[merge_end1 + merge_length:])
            #print(record.seq[:merge_start2] + self.slice2_aa + record.seq[merge_end2 + merge_length:])

            # Once you found a match and wrote the sequence, quit.
            return True

        # if the slice is not matched to any of the FASTA entries, or if the slices are too short,
        # write the slices to an orphan fasta.
        self.logger.info('Slices are not stitchable to fasta. Writing to orphan file for reference.')
        self.logger.info('Slice 1: {0}'.format(self.slice1_aa))
        self.logger.info('Slice 2: {0}'.format(self.slice2_aa))
        self.logger.info('Canonical: {0}'.format(record.seq))

        # format the name of the orphan slice 1 record
        orphan_slice1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                                  id=('xx|ORPHN|' + self.j.gene_symbol + '|'
                                      + self.j.gene_id + '|' + self.j.junction_type + '1|' + self.j.name + '|'
                                      + str(self.j.chr) + '|' + str(self.j.anc_ee) + '|' + str(self.j.alt1_ee)
                                      + '|' + self.translated_strand + str(self.translated_phase) + '|'
                                      + 'r' + str(self.j.min_read_count) + '|' + suffix),
                                  name='Protein name here',
                                  description='Description',)

        # format the name of the orphan slice 2 record
        orphan_slice2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
                                  id=('xx|ORPHN|' + self.j.gene_symbol + '|'
                                      + self.j.gene_id + '|' + self.j.junction_type + '2|' + self.j.name + '|'
                                      + str(self.j.chr) + '|' + str(self.j.anc_ee) + '|' + str(self.j.alt1_ee)
                                      + '|' + self.translated_strand + str(self.translated_phase) + '|'
                                      + 'r' + str(self.j.min_read_count) + '|' + suffix),
                                  name='Protein name here',
                                  description='Description',)

        h.write_seqrecord_to_fasta(orphan_slice1, output, (suffix + '_orphan'))
        h.write_seqrecord_to_fasta(orphan_slice2, output, (suffix + '_orphan'))

        return True
