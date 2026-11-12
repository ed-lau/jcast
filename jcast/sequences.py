# -*- coding: utf-8 -*-

""" Methods that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids """

import re
import logging
import os.path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import requests as rq
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from io import StringIO
import sqlite3 as sq

from jcast import cache
from jcast import helpers as h
from jcast.junctions import Junction
from jcast import params

from jcast.annots import AnnotatedTranscript, gene_rows
from jcast.uniprot import gene_proteins

STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# The UniProt API's answers, by gene id, in the cache directory. 2026-11-12 it was cache/uniprot-cache.db under the
# working directory; see reports/2026-11-12-thresholds-runs-canonical-audit.md §3.5
API_CACHE = 'uniprot-cache.db'
_api_cache_warned = False


def _api_cache_connect():
    """ The API cache, or None (with one warning) if it cannot be opened """
    global _api_cache_warned
    directory = cache.directory()
    if directory is None:
        return None
    try:
        con = sq.connect(os.path.join(directory, API_CACHE))
        con.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')
        return con
    except sq.Error as e:
        if not _api_cache_warned:
            logging.getLogger('jcast').warning(f'Cannot open the UniProt API cache in {directory} ({e}); '
                                               f'running without it.')
            _api_cache_warned = True
        return None


def _api_cache_get(cache_id):
    """ The cached FASTA text for a gene id, or None """
    con = _api_cache_connect()
    if con is None:
        return None
    try:
        row = con.execute('''SELECT id, seq FROM sequences WHERE id=:cache''', {'cache': cache_id}).fetchone()
    except sq.Error:
        row = None
    finally:
        con.close()
    return row[1] if row else None


def _api_cache_put(cache_id, fasta) -> bool:
    """ Cache a gene's FASTA text; False (with one warning) if it cannot be written """
    global _api_cache_warned
    con = _api_cache_connect()
    if con is None:
        return False
    try:
        con.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''', {'id': cache_id, 'seq': fasta})
        con.commit()
        return True
    except sq.Error as e:
        if not _api_cache_warned:
            logging.getLogger('jcast').warning(f'Cannot write to the UniProt API cache ({e}); running without it.')
            _api_cache_warned = True
        return False
    finally:
        con.close()


class Sequence(object):

    def __init__(self,
                 junction: Junction,
                 uniprot=None,
                 ):
        """
        :param junction: the splice junction object
        :param uniprot: UniProtFasta read from --uniprot, or None to get the canonical from the UniProt API

        Mostly copying properties of the junction object (already trimmed) to start a new sequence object. This
        sequence object will be used to make the nucleotide slices and translate into protein sequences.

        """
        self.logger = logging.getLogger('jcast')

        self.j = junction
        self.uniprot = uniprot
        self.slice1_nt = None
        self.slice2_nt = None
        self.slice1_aa = ''
        self.slice2_aa = ''

        self.canonical_aa = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

        self.slice1_stitched = None
        self.slice2_stitched = None

        self.gtf_canonical_transcript = None
        self.gtf_alternative_transcripts = None

        self.translated_phase = -1                  # Phase that was actually used for translation (for logging).
        self.translated_phases = [-1, -1]           # Phase actually used for each slice
        self.alt_suffixes = None                    # Output suffix per slice when routed by alternative termini
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
        by checking whether the two slices' coding lengths (from each slice's own reading frame to its end)
        differ by a non-multiple of 3, i.e. whether the event shifts the frame of the shared downstream sequence.

        2026-11-10 the coding lengths account for each slice's own frame (Junction.slice_phases, PR 4b); comparing
        raw slice lengths misclassified slices that start in different frames. See reports/2026-11-10-alt-termini-audit.md

        :return: True if frameshift, False if not
        """

        offsets = [p if p in [0, 1, 2] else 0 for p in self.j.slice_phases]
        return ((len(self.slice1_nt) - offsets[0]) - (len(self.slice2_nt) - offsets[1])) % 3 != 0

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

        if self.j.junction_type in ['MXE', 'SE']:
            self.slice1_nt = anc_nt + alt1_nt + down_nt
            self.slice2_nt = anc_nt + alt2_nt + down_nt

        # 2026-11-10 the rMATS retained exon (alt1) already spans both flanking exons, so it is slice 1 on its own;
        # previously both flanks were added again. See reports/2026-11-10-ri-slice-audit.md
        elif self.j.junction_type == 'RI':
            self.slice1_nt = alt1_nt
            self.slice2_nt = anc_nt + down_nt

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
                  log: bool = True,
                  ):
        """
        This is the translate function for tier 1/2 peptides. Each slice is translated in its own annotated frame
        (Junction.slice_phases) with make_pep's terminate option, which returns an empty string when it runs into
        a stop codon before the last codon. A slice without an annotated frame is not translated.

        2026-11-10 removed the three-frame fallback (and with it Tier 3): it guessed a frame when the phase was
        missing or did not translate, kept the last frame tried, and its output was not used downstream.
        See reports/2026-11-10-phase-audit.md

        :param log: bool Whether to write results to log file
        :return: True
        """

        phases = self.j.slice_phases
        self.slice1_aa, self.slice2_aa = [h.make_pep(nt.seq, self.j.strand, p, terminate=True)
                                          if p in [0, 1, 2] else ''
                                          for nt, p in zip([self.slice1_nt, self.slice2_nt], phases)]
        self.translated_phases = [p if p in [0, 1, 2] else -1 for p in phases]
        self.translated_phase = self.j.phase_label

        if log:
            for i_, s_ in enumerate([self.slice1_aa, self.slice2_aa]):
                if s_ != '':
                    self.logger.info(
                        'Translated AA for {0} {1} {2} phase {3}{4} slice {5}: {6}'.format(
                            self.j.junction_type,
                            self.j.name,
                            self.j.gene_symbol,
                            self.j.strand,
                            self.translated_phases[i_],
                            i_ + 1,
                            s_,
                        )
                    )
        return True

    def route_alternative_termini(self) -> bool:
        """
        Recognise alternative termini that would otherwise be read as premature stops or as start-less fragments,
        and, if the junction qualifies, translate and join its slices accordingly.

        A slice is C-anchored if its first in-frame stop codon is an annotated protein-coding stop codon of the gene
        (internal or in its last codon): it ends there. A slice is N-anchored if it lacks the gene's start codon
        while the other slice contains it, and it contains another annotated protein-coding start codon: it is
        translated from the first such ATG. A C-anchored slice is joined to the canonical at its N-terminus only,
        an N-anchored one at its C-terminus only, one anchored at both ends not at all.

        The junction qualifies if at least one slice is anchored, both slices reach an acceptable end (no stop
        codon, a stop in the last codon, or an annotated stop) and it is not a frameshift. The frameshift test is
        waived if a slice ends at an internal annotated stop (the slices then end at different stops); otherwise it
        compares the coding lengths from each slice's translation start.

        Each slice's output suffix depends on which terminus differs from the canonical's (first and last
        params.aa_stitch_length residues): T1, T1_alt_cterm, T1_alt_nterm or T1_alt_both; an anchored slice that
        cannot be joined goes to the orphan file of its anchor.

        See reports/2026-11-10-alt-termini-audit.md

        :return: True if the junction was routed (slices, stitched products and alt_suffixes set)
        """

        gene_start = self.j.tx0 if self.j.strand == '+' else self.j.tx1
        positions = self.j.slice_positions()
        nts = [self.slice1_nt, self.slice2_nt]
        seqs = [str(nt.seq).upper() for nt in nts]
        if self.j.strand == '-':
            seqs = [h.get_complementary(s) for s in seqs]
        contains_start = [gene_start is not None and gene_start in set(p) for p in positions]

        plans = []
        for k in range(2):
            s, pos = seqs[k], positions[k]
            phase = self.j.slice_phases[k]
            start, n_anchored = (phase if phase in [0, 1, 2] else None), False

            if not contains_start[k] and contains_start[1 - k]:
                alt = next((i for i, p in enumerate(pos) if p in self.j.start_sites and s[i:i + 3] == 'ATG'), None)
                if alt is not None:
                    start, n_anchored = alt, True

            if start is None:
                return False

            stop = next((i for i in range(start, len(s) - 2, 3) if s[i:i + 3] in STOP_CODONS), None)
            terminal = stop is not None and stop > len(s) - 5
            annotated = stop is not None and pos[stop] in self.j.stop_sites
            if stop is not None and not terminal and not annotated:
                return False

            plans.append({'start': start, 'n': n_anchored, 'c': annotated, 'internal': annotated and not terminal,
                          'length': len(s)})

        if not any(p['n'] or p['c'] for p in plans):
            return False
        if not any(p['internal'] for p in plans) and \
                ((plans[0]['length'] - plans[0]['start']) - (plans[1]['length'] - plans[1]['start'])) % 3 != 0:
            return False

        canonical = str(self.canonical_aa.seq)
        length = params.aa_terminus_length
        suffixes, products, peptides = [], [], []
        for k, p in enumerate(plans):
            aa = h.make_pep(nts[k].seq, self.j.strand, p['start'], terminate=False)
            product = self._join(aa, {(True, True): 'none', (True, False): 'c',
                                      (False, True): 'n', (False, False): 'both'}[(p['n'], p['c'])],
                                 slice_index=k, start=p['start'])
            if product is None:
                suffix = {(True, True): 'T1_alt_both', (True, False): 'T1_alt_nterm',
                          (False, True): 'T1_alt_cterm', (False, False): 'T1'}[(p['n'], p['c'])]
            elif str(product.seq) == canonical:
                suffix = 'T1'
            else:
                n_same = str(product.seq)[:length] == canonical[:length]
                c_same = str(product.seq)[-length:] == canonical[-length:]
                suffix = {(True, True): 'T1', (True, False): 'T1_alt_cterm',
                          (False, True): 'T1_alt_nterm', (False, False): 'T1_alt_both'}[(n_same, c_same)]
            peptides.append(aa)
            products.append(product)
            suffixes.append(suffix)
            self.translated_phases[k] = p['start'] % 3

        self.slice1_aa, self.slice2_aa = peptides
        self.slice1_stitched, self.slice2_stitched = products
        self.alt_suffixes = suffixes
        p1, p2 = self.translated_phases
        self.translated_phase = str(p1) if p1 == p2 else f'{p1}/{p2}'
        return True

    def _join(self, aa: str, mode: str, slice_index: int = None, start: int = None):
        """
        Join a peptide to the canonical: 'both' ends, 'both_or_n' (both ends if possible, else the N-terminus only),
        'n' (N-terminus only), 'c' (C-terminus only) or 'none' (the peptide alone). None if there is no canonical,
        the peptide is shorter than the anchor, or an anchor is not found.

        The peptide is placed by its first and last params.aa_stitch_length residues (the anchors). An anchor is used
        only if it occurs exactly once in the canonical, and a join at both ends needs the start anchor not to be
        after the end anchor (REVIEW.md MED-3). If slice_index and start (the peptide's translation start in that
        slice) are given, an end that is already complete is not joined: a peptide ending at a stop codon has its
        C-terminus, one starting at the gene's start codon its N-terminus.

        See reports/2026-11-11-stitch-anchor-audit.md

        :return: the joined SeqRecord, or None
        """

        canonical = self.canonical_aa
        length = params.aa_stitch_length
        if mode == 'none':
            return canonical[:0] + aa
        if len(canonical) == 0 or len(aa) < length:
            return None

        n_complete = False
        if slice_index is not None and start in range(len(self.slice_nt(slice_index))) and mode != 'c':
            n_complete = self.j.slice_positions()[slice_index][start] == \
                (self.j.tx0 if self.j.strand == '+' else self.j.tx1)
            if self.ends_at_stop(slice_index, start, aa):
                mode = 'n'

        def unique(anchor):
            i = canonical.seq.find(anchor)
            return i if i != -1 and canonical.seq.find(anchor, i + 1) == -1 else -1

        m0, m1 = unique(aa[:length]), unique(aa[-length:])
        before = canonical[:0] if n_complete else (canonical[:m0] if m0 != -1 else None)
        if mode in ('both', 'both_or_n') and before is not None and m1 != -1 and (n_complete or m0 <= m1):
            return before + aa + canonical[m1 + length:]
        if mode in ('n', 'both_or_n') and before is not None:
            return before + aa
        if mode == 'c' and m1 != -1:
            return canonical[:0] + aa + canonical[m1 + length:]
        return None

    def slice_nt(self, slice_index: int) -> str:
        """ a slice's nucleotides in transcript orientation, uppercase """
        nt = str([self.slice1_nt, self.slice2_nt][slice_index].seq).upper()
        return nt if self.j.strand == '+' else h.get_complementary(nt)

    def ends_at_stop(self, slice_index: int, start: int, aa: str) -> bool:
        """ whether the peptide translated from start in the slice is followed by a stop codon """
        end = start + 3 * len(aa)
        return self.slice_nt(slice_index)[end:end + 3] in STOP_CODONS

    def translate_forced(self,
                         slice_to_translate: int,
                         ) -> None:
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
        phase = self.j.slice_phases[slice_to_translate - 1]
        seq = h.make_pep(nt_to_translate.seq, self.j.strand, phase, terminate=False) if phase in [0, 1, 2] else ''

        if len(seq) > 0:
            self.translated_phases[slice_to_translate - 1] = phase
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

        return None

    def get_annotated_transcripts_gtf(self, gtf):
        """
        get the annotated transcripts of the junction under question

        :param gtf: genome annotation

        """

        # first, subset the gtf by the gene under question
        gtf_gene = gene_rows(gtf, self.j.gene_id)

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

        # a supplied UniProt FASTA, falling back to the GTF (reports/2026-11-11-uniprot-fasta-evaluation.md)
        if self.uniprot is not None:
            self.canonical_aa = self.uniprot.canonical(
                self.j.gene_id,
                self.j.gene_symbol,
                proteins=lambda: gene_proteins(gtf, genome_index, self.j.gene_id),
                gtf_canonical=lambda: self.get_canonical_aa_gtf(gtf, genome_index),
            )

        elif params.use_gtf_only:
            self.canonical_aa = self.get_canonical_aa_gtf(gtf,
                                                          genome_index)
            # Uniprot fallback
            # if len(self.canonical_aa) ==0:
            #    self.canonical_aa = self.get_canonical_aa_uniprot()

        # try to get the reviewed canonical entry from uniprot; if no reviewed entry, get un-reviewed entry
        else:
            self.canonical_aa = self.get_canonical_aa_uniprot(reviewed='true')
            if len(self.canonical_aa[:]) == 0:
                self.canonical_aa = self.get_canonical_aa_uniprot(reviewed='false')

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
                                 reviewed: str = 'true',
                                 ) -> SeqRecord:
        """
        get the canonical sequences from Uniprot

        :param reviewed: get only reviewed sequence
        :return: canonical aa seqrecord
        """

        record = SeqRecord(Seq(''), annotations={'molecule_type': 'extended_protein'})

        # cache retrieved sequences if available, otherwise retrieve from Ensembl
        cache_id = self.j.gene_id
        read_fasta = _api_cache_get(cache_id)

        if read_fasta:
            record = list(SeqIO.parse(StringIO(read_fasta), 'fasta'))[0]
            self.logger.info('Locally cached sequence retrieved for {0}.'.format(cache_id))

        else:

            server = 'https://www.ebi.ac.uk'

            # Uniprot does not support transcript version
            ext = '/proteins/api/proteins/Ensembl:' + re.sub('\\..*', '', self.j.gene_id) + \
                  '?offset=0&size=1&reviewed=' + reviewed + '&isoform=0'

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

            if ret.status_code == 200 and ret.text != '':
                record = list(SeqIO.parse(StringIO(ret.text), 'fasta'))[0]
                if _api_cache_put(cache_id, record.format('fasta')):
                    self.logger.info('Sequence retrieved from Uniprot and written into local cache.')

            elif ret.status_code == 200 and ret.text == '':
                self.logger.info('Retrieved empty fasta from Ensembl for {0}'.format(self.j.gene_id))
                # TODO: A known issue where Uniprot sequences do not exist for some Ensembl genes
                # TODO: Rather than fix this we will simply use the GTF in future versions.

            elif ret.status_code != 200:
                self.logger.warning('Retrieval of protein sequence failed.')

        return record

    def stitch_to_canonical_aa(self,
                                slice_to_stitch: int,
                                slice_has_ptc: bool = False,
                                ) -> None:
        """
        Given a translated junction sequence, look for the fasta entry that overlaps with it, then return the entry
        and the coordinates. This will be used to extend the junction sequence to encompass entire protein sequence.

        While this does not help resolve junction sequences or isoforms, this helps match the database
        file as closely as possible to the actual spectra in a mass spectrometry experiment to strengthen the
        hypotheses being tested in the database search.

        The slice is placed by its anchors, and a complete end is not joined (see _join).

        :param slice_to_stitch: int     1 or 2
        :param slice_has_ptc:   bool    allows stitching to only n-terminus
        :return: True
        """


        if slice_to_stitch not in [1, 2]:
            raise Exception('slice must be 1 or 2.')
        assert type(params.aa_stitch_length) is int and params.aa_stitch_length >= 6, \
            'Stitch length must be integer and at least 6'

        k = slice_to_stitch - 1
        stitched = self._join([self.slice1_aa, self.slice2_aa][k], 'both_or_n' if slice_has_ptc else 'both',
                              slice_index=k, start=self.translated_phases[k])

        if stitched:
            if slice_to_stitch == 1:
                self.slice1_stitched = stitched[:]
            elif slice_to_stitch == 2:
                self.slice2_stitched = stitched[:]

        return None

    def write_canonical(self,
                        out_dir,
                        ) -> None:
        """
        write canonical sequences to fasta file
        :param out_dir: output directory
        """

        canonical = self.canonical_aa[:]

        if len(self.canonical_aa[:]) == 0 and self.uniprot is None:
            # try to get the reviewed sequence
            self.get_canonical_aa_uniprot(reviewed='true')
            canonical = self.canonical_aa[:]
            # try to get unreviewed sequence
            if len(self.canonical_aa[:]) == 0:
                self.get_canonical_aa_uniprot(reviewed='false')

        if len(self.canonical_aa[:]) > 0:
            h.write_seqrecord_to_fasta(canonical, out_dir, 'canonical')

        return None

    def write_slices(self,
                     out_dir,
                     suffix,
                     slice_suffixes=None,
                     ) -> None:
        """
        write the translated SeqRecord objects created into fasta file

        :param out_dir:          string  output directory
        :param suffix:          string  additional suffix to add to the end of an output file
        :param slice_suffixes:  list    optional suffix per slice, overriding suffix (alternative termini)

        :return:
        """

        for i, stitched in enumerate([self.slice1_stitched, self.slice2_stitched]):
            slice_suffix = slice_suffixes[i] if slice_suffixes else suffix

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
                                           self.translated_phases[i],
                                           self.j.sum_sjc,
                                           slice_suffix,
                                       ),
                                       annotations={'molecule_type': 'extended_protein'},
                                       name=f'Sequence_{str(self.j.gene_symbol)}',
                                       description='Orphan',
                                       )
                    h.write_seqrecord_to_fasta(orphan, out_dir, (slice_suffix + '_orphan'))

            # if the stitched is same as canonical, write canonical
            elif stitched.seq == self.canonical_aa.seq:
                self.write_canonical(out_dir=out_dir)

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
                    self.translated_phases[i],
                    self.j.sum_sjc,
                    slice_suffix,
                )
                h.write_seqrecord_to_fasta(stitched, out_dir, slice_suffix)

        return None
