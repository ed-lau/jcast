# -*- coding: utf-8 -*-

""" Methods that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids """

import re
import logging
import os.path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC

import requests as rq
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from io import StringIO
import sqlite3 as sq

from jcast import helpers as h
from jcast.junctions import Junction
from jcast import params

from jcast.annots import AnnotatedTranscript


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

        self.canonical_aa = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

        self.slice1_stitched = None
        self.slice2_stitched = None

        self.gtf_canonical_transcript = None
        self.gtf_alternative_transcripts = None

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

    def get_annotated_transcripts_gtf(self, gtf):
        """
        get the annotated transcripts of the junction under question

        :param gtf: genome annotation

        """

        # first, subset the gtf by the gene under question
        gtf_gene = gtf.annot.query('gene_id == @self.j.gene_id')

        # then get a list of all transcripts that are protein_coding, has coding sequences,
        # and has transcript support level equal to or below the threashold

        tsl = params.tsl_threshold
        filter_tx = gtf_gene.query('feature == "CDS" & '
                                            'transcript_biotype == "protein_coding" & '
                                            'transcript_support_level <= @tsl')
        coding_tx = sorted(filter_tx['transcript_name'].drop_duplicates())

        # for each qualifying transcripts, get all the CDS

        annotated_transcripts = []

        for tx in coding_tx:
            tx_tbl = filter_tx.query('transcript_name == @tx').sort_values(by=['exon_number'])
            tx_exons = [ex for ex in zip(tx_tbl['start'], tx_tbl['end'])]
            tx_protein_id = tx_tbl['protein_id'].drop_duplicates()


            try:
                tx_start_codon = gtf_gene.query('feature == "start_codon" & '
                                                'transcript_name == @tx')[['start', 'end']].values.tolist()[0]
                tx_end_codon = gtf_gene.query('feature == "stop_codon" & '
                                                'transcript_name == @tx')[['start', 'end']].values.tolist()[0]
                tx_frame = gtf_gene.query('feature == "start_codon" & '
                                          'transcript_name == @tx')['frame'].drop_duplicates().values[0]

            # some tsl-1 protein_coding CDS (from havana) have no start/stop codons in gtf. skip these for now.
            # note that some have stop codon but not start codon. we should deal with these cases.
            except IndexError:
                # tx_start_codon = [None, None]
                # tx_end_codon = [None, None]
                # tx_frame = None
                continue
            # TODO: these pandas queries are ugly - any tidier way to get the cell values?

            # check if the transcript is already the same as one with exactly the same exon arrangements
            # if not then add the transcript to the pile
            already_exists = False
            if len(annotated_transcripts) > 0:
                for at in annotated_transcripts:
                    if tx_exons == at.exons:
                        already_exists = True
                        at.transcript_name += ', {0}'.format(tx)

            if not already_exists:
                annotated_transcripts += [AnnotatedTranscript(transcript_name=tx,
                                                              protein_id=tx_protein_id,
                                                              exons=tx_exons,
                                                              start_codon=tx_start_codon,
                                                              end_codon=tx_end_codon,
                                                              starting_translation_phase=tx_frame,
                                            )
                ]

        # the longest coding transcript is assumed to be canonical
        if len(annotated_transcripts) > 0:
            longest_tx_length = max([len(tx) for tx in annotated_transcripts])
            canonical_transcript = [tx for tx in annotated_transcripts if len(tx) == longest_tx_length][0]
            alternative_transcripts = [tx for tx in annotated_transcripts if tx != canonical_transcript]

            self.gtf_canonical_transcript = canonical_transcript
            self.gtf_alternative_transcripts = alternative_transcripts

        return True

    def get_canonical_aa(self,
                         gtf,
                         genome_index,
                         ):
        """
        decide whether to get canonical amino acids by GTF or by Uniprot
        :param gtf: gtf annotation
        :return: True
        """

        if params.use_gtf_only:
            self.canonical_aa = self.get_canonical_aa_gtf(gtf,
                                                          genome_index)

            # Uniprot fallback
            # if len(self.canonical_aa) ==0:
            #    self.canonical_aa = self.get_canonical_aa_uniprot()

        else:
            self.canonical_aa = self.get_canonical_aa_uniprot()


        return True


    def translate_annotated_transcriptfrom_gtf(self,
                                               annotated_transcript,
                                               genome_index):
        """
        translate a collection of exons from an annotated transcript
        # TODO: should this be a method instead for AnnotatedTranscripts

        """

        # note that if the strand is '-', the exons should be reversed
        # beware of this when inserting exons into transcripts later ...
        if self.j.strand == '-':
            exons = [ex for ex in reversed(annotated_transcript.exons)]
        else:
            exons = annotated_transcript.exons

        canonical_nt = ''
        for start, end in exons:
            canonical_nt += h.get_local_nuc(genome_index, self.j.chr, start, end)

        canonical_aa = h.make_pep(canonical_nt.seq,
                                  strand=self.j.strand,
                                  phase=self.gtf_canonical_transcript.starting_translation_phase,
                                  terminate=True)



        canonical_aa = SeqRecord(Seq(canonical_aa), annotations={'molecule_type': 'extended_protein'})
        canonical_aa.id = '{0}|{1}'.format(self.j.gene_id,
                                           self.j.gene_symbol)

        return canonical_aa

    def get_canonical_aa_gtf(self,
                             gtf,
                             genome_index,
                             ) -> SeqRecord:
        """
        get the canonical amino acid sequence from gtf

        :return: canonical aa str
        """

        if self.gtf_canonical_transcript is None:
            self.get_annotated_transcripts_gtf(gtf)

        # if nothing from the gtf, return an empty record
        if self.gtf_canonical_transcript is None:
            return SeqRecord(Seq(''), annotations={'molecular_type': 'extended_protein'})

        return self.translate_annotated_transcriptfrom_gtf(annotated_transcript=self.gtf_canonical_transcript,
                                                           genome_index=genome_index)


    def get_canonical_aa_uniprot(self,
                                 ) -> SeqRecord:
        """
        get the canonical sequences from Uniprot

        :return: canonical aa seqrecord
        """

        record = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

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
            record = list(SeqIO.parse(StringIO(read_fasta[1]), 'fasta'))[0]
            self.logger.info('Locally cached sequence retrieved for {0}.'.format(cache))
            con.close()

        else:

            server = 'https://www.ebi.ac.uk'

            # Uniprot does not support transcript version
            ext = '/proteins/api/proteins/Ensembl:' + re.sub('\\..*', '', self.j.gene_id) + \
                  '?offset=0&size=1&reviewed=true&isoform=0'

            self.logger.info('Sequence not cached locally. Attempting to get from Uniprot: {0}'.format(
                server + ext
            ))

            retries = Retry(total=params.uniprot_max_retries,
                            backoff_factor=0.1,
                            status_forcelist=[500, 502, 503, 504])

            rqs = rq.Session()
            rqs.mount('https://', HTTPAdapter(max_retries=retries))
            ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

            if not ret.ok:
                self.logger.warning('Failed retrieval for {0} after {1} retries.'.format(self.j.gene_id,
                                                                                         params.uniprot_max_retries)
                                    )
                con.close() # TODO: change to with statement

            if ret.status_code == 200 and ret.text != '':
                record = list(SeqIO.parse(StringIO(ret.text), 'fasta'))[0]

                # TODO: Catch sqlite3.OperationalError if writing fails, such as in a remote volume.
                cur.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''',
                            {'id': cache, 'seq': record.format('fasta')})
                con.commit()
                con.close()
                self.logger.info('Sequence retrieved from Uniprot and written into local cache.')

            elif ret.status_code == 200 and ret.text == '':
                self.logger.info('Retrieved empty fasta from Ensembl for {0}'.format(self.j.gene_id))
                # TODO: A known issue where Uniprot sequences do not exist for some Ensembl genes
                # TODO: Rather than fix this we will simply use the GTF in future versions.
                con.close()

            elif ret.status_code != 200:
                self.logger.warning('Retrieval of protein sequence failed.')
                con.close()

            else:
                con.close()
                pass

        return record

    def stitch_to_canonical_aa(self,
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

        # don't stitch if there is no canonical
        if len(canonical) == 0:
            return True

        # stitch length
        stitch_length = params.aa_stitch_length
        assert type(stitch_length) is int and stitch_length >= 6, 'Stitch length must be integer and at least 6'

        # which slice to stitch
        if slice_to_stitch == 1:
            slice_ = self.slice1_aa
        elif slice_to_stitch == 2:
            slice_ = self.slice2_aa
        else:
            raise Exception('slice must be 1 or 2.')

        # if the slice is too short, do not stitch
        if len(slice_) < stitch_length:
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

                    orphan = SeqRecord(Seq(unstitched),
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
                                           self.j.sum_sjc,
                                           suffix,
                                       ),
                                       annotations={'molecule_type': 'extended_protein'},
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
                    self.j.sum_sjc,
                    suffix,
                )
                h.write_seqrecord_to_fasta(stitched, outdir, suffix)

        return True
