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
        self.logger = logging.getLogger('jcast.seq')

        self.j = junction
        self.slice1_nt = None
        self.slice2_nt = None
        self.slice1_aa = ''
        self.slice2_aa = ''
        self.canonical_aa = self.get_canonical_aa_uniprot()
        self.slice1_stitched = None
        self.slice2_stitched = None

        self.translated_phase = -1                  # Phase that was actually used for translation.
        self.translated_strand = self.j.strand    # Strand that was actually used for translation

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
            self.slice1_nt = anc_nt + alt1_nt + down_nt
            self.slice2_nt = anc_nt + alt2_nt + down_nt

        # 2020-07-25 the flanking exons for A5SS and A3SS:
        elif self.j.junction_type == 'A5SS':
            if self.j.strand == '+':
                self.slice1_nt = alt1_nt + anc_nt
                self.slice2_nt = alt2_nt + anc_nt
            elif self.j.strand == '-':
                self.slice1_nt = anc_nt + alt1_nt
                self.slice2_nt = anc_nt + alt2_nt

        elif self.j.junction_type == 'A3SS':
            if self.j.strand == '+':
                self.slice1_nt = anc_nt + alt1_nt
                self.slice2_nt = anc_nt + alt2_nt
            elif self.j.strand == '-':
                self.slice1_nt = alt1_nt + anc_nt
                self.slice2_nt = alt2_nt + anc_nt

        for i_, s_ in enumerate([self.slice1_nt, self.slice2_nt]):
            self.logger.info('Retrieved nucleotide for {0} {1} {2} slice {3}: {4}'.format(self.j.junction_type,
                                                                                          self.j.name,
                                                                                          self.j.gene_symbol,
                                                                                          i_+1,
                                                                                          s_.seq))

    def translate(self,
                  use_phase: bool = True,
                  log: bool = True,
                  ):
        """
        This is the translate function for tier 1/2/3 peptide. Calls make_pep with the terminal option
        which will terminate when it runs into a stop codon. Later on we should make a force_translate that does
        three frame translation and just either return all frames or return the longest.

        :param use_phase: bool Whether to ue the stored phase or attempt to do three-frame translation
        :param log: bool Whether to write results to log file
        :return: True
        """

        if self.j.phase in [0, 1, 2] and use_phase:
            self.slice1_aa = h.make_pep(self.slice1_nt.seq, self.j.strand, self.j.phase, terminate=True)
            self.slice2_aa = h.make_pep(self.slice2_nt.seq, self.j.strand, self.j.phase, terminate=True)
            self.logger.debug('Used retrieved phase for {0} {1}:'.format(self.j.name, self.j.gene_symbol))
            self.translated_phase = self.j.phase

        else:
            for i in range(3):
                self.slice1_aa = h.make_pep(self.slice1_nt.seq, self.j.strand, i, terminate=True)
                self.slice2_aa = h.make_pep(self.slice2_nt.seq, self.j.strand, i, terminate=True)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    self.translated_phase = i
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        if log:
            for i_, s_ in enumerate([self.slice1_aa, self.slice2_aa]):
                if s_ != '':
                    self.logger.info(
                        'Translated AA for {0} {1} {2} phase {3}{4} slice {5} (use_phase={6}): {7}'.format(
                            self.j.junction_type,
                            self.j.name,
                            self.j.gene_symbol,
                            self.j.strand,
                            self.translated_phase,
                            i_ + 1,
                            use_phase,
                            s_,
                        )
                    )
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

        # Determine whether slice 1 or slice 2 nucleotides is being force-translated
        if slice_to_translate == 1:
            nt_to_translate = self.slice1_nt
        elif slice_to_translate == 2:
            nt_to_translate = self.slice2_nt
        else:
            nt_to_translate = None

        # translate without terminating at stop codon
        seq = h.make_pep(nt_to_translate.seq, self.j.strand, self.j.phase, terminate=False)

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

        self.logger.info('Force-Translated AA for {0} {1} {2} phase {3}{4} slice {5}: {6}'.format(self.j.junction_type,
                                                                                                  self.j.name,
                                                                                                  self.j.gene_symbol,
                                                                                                  self.j.strand,
                                                                                                  self.translated_phase,
                                                                                                  slice_to_translate,
                                                                                                  force_translated_aa,
                                                                                                  )
                         )

        return True

    def get_canonical_aa_uniprot(self,
                                 ) -> SeqRecord:
        """
        get the canonical sequences from Uniprot

        :return: canonical aa seqrecord
        """

        record = SeqRecord(Seq('', IUPAC.extended_protein))

        # cache retrieved sequences if available, otherwise retrieve from Ensembl
        cache = self.j.gene_id

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
            self.logger.info('Locally cached sequence retrieved for {0}.'.format(cache))
            con.close()

        else:

            server = 'https://www.ebi.ac.uk'
            ext = '/proteins/api/proteins/Ensembl:' + self.j.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

            self.logger.info('Sequence not cached locally. Attempting to get from Uniprot: {0}'.format(
                server + ext
            ))

            retries = Retry(total=params.max_retries,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            rqs = rq.Session()
            rqs.mount('https://', HTTPAdapter(max_retries=retries))
            ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

            if not ret.ok:
                self.logger.warning('Failed retrieval for {0} after {1} retries.'.format(self.j.gene_id,
                                                                                         params.max_retries)
                                    )

            if ret.status_code == 200 and ret.text != '':
                record = list(SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein))[0]

                cur.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''',
                            {'id': cache, 'seq': record.format('fasta')})
                con.commit()
                con.close()
                self.logger.info('Sequence retrieved from Uniprot and written into local cache.')

            elif ret.status_code == 200 and ret.text == '':
                self.logger.info('Retrieved empty fasta from Ensembl for {0}'.format(self.j.gene_id))
                # TODO: A known issue where Uniprot sequences do not exist for some Ensembl genes
                # TODO: Rather than fix this we will simply use the GTF in future versions.

            elif ret.status_code != 200:
                self.logger.warning('Retrieval of protein sequence failed.')

        return record

    # TODO: save all sequences and write once
    def stitch_to_canonical(self,
                            slice_to_stitch: int,
                            slice_has_ptc: bool = False,
                            ):
        """
        Given a translated junction sequence, look for the fasta entry that overlaps with it, then return the entry
        and the coordinates. This will be used to extend the junction sequence to encompass entire protein sequence.

        While this does not help resolve junction sequences or isoforms, this helps match the database
        file as closely as possible to the actual spectra in a mass spectrometry experiment to strengthen the
        hypotheses being tested in the database search.

        :param slice_to_stitch: int     1 or 2
        :param slice_has_ptc:   bool    allows stitching to only n-terminus
        :return: True
        """

        canonical = self.canonical_aa[:]
        if len(canonical) == 0:
            self.get_canonical_aa_uniprot()
            canonical = self.canonical_aa[:]

        # stitch length
        stitch_length = params.stitch_length
        assert type(stitch_length) is int and stitch_length >= 6, 'Stitch length must be integer and at least 6'

        # which slice to stitch
        if slice_to_stitch == 1:
            slice_ = self.slice1_aa
        elif slice_to_stitch == 2:
            slice_ = self.slice2_aa
        else:
            raise Exception('slice must be 1 or 2.')

        # if the slice is too short, do not stitch
        if len(slice_) < params.stitch_length:
            return True

        # find out where the first (stitch_length) amino acids meets the UniProt canonical sequences
        merge0, merge1 = canonical.seq.find(slice_[:stitch_length]), canonical.seq.find(slice_[-stitch_length:])

        # slice stitches to canonical at both ends:
        if merge0 != -1 and merge1 != -1:
            stitched = canonical[:merge0] + slice_ + canonical[merge1+stitch_length:]

        # slice has ptc but beginning stitches to canonical
        elif slice_has_ptc and merge0 != -1:
            stitched = canonical[:merge0] + slice_

        else:
            stitched = None

        if stitched:
            if slice_to_stitch == 1:
                self.slice1_stitched = stitched[:]
            elif slice_to_stitch == 2:
                self.slice2_stitched = stitched[:]

        return True

    def write_canonical(self,
                        outdir,
                        ):
        """
        write canonical sequences to fasta file
        :param output: output directory
        """

        canonical = self.canonical_aa[:]
        if len(self.canonical_aa[:]) == 0:
            self.get_canonical_aa_uniprot()
            canonical = self.canonical_aa[:]
        h.write_seqrecord_to_fasta(canonical, outdir, 'canonical')

        return True

    def write_slices(self,
                     outdir,
                     suffix,
                     ):
        """
        write the translated SeqRecord objects created into fasta file

        :param outdir:          string  output directory
        :param suffix:          string  additional suffix to add to the end of an output file

        :return:
        """

        for i, stitched in enumerate([self.slice1_stitched, self.slice2_stitched]):

            # if no stitched aa, write the slice aa to orphan if it exists
            if stitched is None:
                unstitched = [self.slice1_aa, self.slice2_aa][i]
                if len(unstitched) > 0:

                    self.logger.info('Slice {0} does not stitch to canonical sequence. Writing to orphan file.'.format(
                        i+1
                    ))

                    orphan = SeqRecord(Seq(unstitched, IUPAC.extended_protein),
                                       id='xx|ORPHN|{0}|{1}|{2}{3}|{4}|{5}|{6}:{7}|{8}:{9}|{10}{11}|r{12}|{13}'.format(
                                           self.j.gene_symbol,
                                           self.j.gene_id,
                                           self.j.junction_type,
                                           i+1,
                                           self.j.name,
                                           self.j.chr,
                                           self.j.anc_es,
                                           self.j.anc_ee,
                                           self.j.alt1_es,
                                           self.j.alt1_ee,
                                           self.translated_strand,
                                           self.translated_phase,
                                           self.j.min_read_count,
                                           suffix,
                                       ),
                                       name=self.j.gene_symbol,
                                       description='Orphan',
                                       )
                    h.write_seqrecord_to_fasta(orphan, outdir, (suffix + '_orphan'))

            # if the stitched is same as canonical, write canonical
            elif stitched.seq == self.canonical_aa.seq:
                self.write_canonical(outdir=outdir)

            # otherwise write alternative slice to fasta
            elif len(stitched) > 0:
                stitched.id = self.canonical_aa.id + '|{0}|{1}{2}|{3}|{4}|{5}:{6}|{7}:{8}|{9}{10}|r{11}|{12}'.format(
                    self.j.gene_id,
                    self.j.junction_type,
                    i+1,
                    self.j.name,
                    self.j.chr,
                    self.j.anc_es,
                    self.j.anc_ee,
                    self.j.alt1_es,
                    self.j.alt1_ee,
                    self.translated_strand,
                    self.translated_phase,
                    self.j.min_read_count,
                    suffix,
                )
                h.write_seqrecord_to_fasta(stitched, outdir, suffix)

        return True
