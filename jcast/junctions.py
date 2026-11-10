# -*- coding: utf-8 -*-

""" Methods that concern splice junctions - getting their coordinates, transcription starts/ends, phases. """
import collections
import logging
import os.path
import pandas as pd
import numpy as np

from jcast import params


class RmatsResults(object):
    """
    Container to hold the rMATS output folder and create individual objects based on the five splice types.

    Note that each splice event will be defined by the anchor exon (anc), potential alternative exons (alt1, alt2), and
    the downstream exon (down). For some splice type, some of these will not be present, for example, Skipped Exons
    either contain the alt1 exon or no alt1 exon in the two alternative forms, and no alt2 is present.

    """

    def __init__(self,
                 logger,
                 rmats_dir,

                 ):
        """
        :param rmats_dir: directory containing the rMATS output files
        :param logger: logger object
        """

        self.logger = logger

        self.dir = rmats_dir
        self.rmats_mxe = self._read_rmats_mxe()
        self.rmats_se = self._read_rmats_se()
        self.rmats_ri = self._read_rmats_ri()
        self.rmats_a5ss = self._read_rmats_a5ss()
        self.rmats_a3ss = self._read_rmats_a3ss()

        self.rmats_list = [self.rmats_mxe, self.rmats_se, self.rmats_ri, self.rmats_a5ss, self.rmats_a3ss]
        self.sum_sjc_array = None



    def _read_rmats_mxe(self):
        """
        Read input data frame for Mutually Exclusive Exons (MXE)
        The slices should be anc-alt1-down and anc-alt2-down

        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'MXE.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                      'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'anc_es',
                      'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                      'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                      'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'MXE'

        return df

    def _read_rmats_se(self):
        """
        Read input data frame for Skipped Exons (SE)
        The slices should be anc-down and anc-alt1-down

        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'SE.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'SE'
        df['alt2_es'] = -1
        df['alt2_ee'] = -1

        return df

    def _read_rmats_ri(self):
        """
        Read input data frame for Retained Introns (RI)
        The slices should be anc-down and anc-alt1-down

        :return:
        """
        df = pd.read_table(os.path.join(self.dir, 'RI.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                 'alt1_es', 'alt1_ee', 'anc_es',
                                 'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                 'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                 'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'RI'
        df['alt2_es'] = -1
        df['alt2_ee'] = -1

        return df

    def _read_rmats_a5ss(self):
        """
        Read input data frame for Alternative 5' Splice Sites (A5SS)
        Note this splice type is without the 'downstream' exon, but the anchor (flanking) is downstream.
        The slices should be alt1-anc and alt2-anc for
        2020-07-235 A5SS may be treated differently here: I think for (-) strand, the slices are
        anc-alt1 and anc-alt2
        Note if the strand is +, alt1_es and alt2_es should be identical and before anchor in genomic position.
        Note if the strand is -, alt1_ee and alt2_ee should be the same and after anchor in genomic position.
        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'A5SS.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'A5SS'
        df['down_es'] = -1
        df['down_ee'] = -1

        return df

    def _read_rmats_a3ss(self):
        """
        Read input data frame for Alternative 3' Splice Sites (A3SS)
        Note this splice type is without the downstream exon, the anchor is the upstream.
        The slices are anc-alt1 and anc-alt2.
        2020-07-235 A3SS may be treated differently here: I think for (-) strand, the slices are
        alt1-anc and alt2-anc
        Note if the strand is +, alt1_ee and alt2_ee are identical and after anchor in genomic position.
        If the stand is -, alt1_es and alt2_es are the same, and they are both before anchor in genomic position.
        :return:
        """

        df = pd.read_table(os.path.join(self.dir, 'A3SS.MATS.JC.txt'), sep='\t', low_memory=False)
        df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
                                   'alt1_es', 'alt1_ee', 'alt2_es',
                                   'alt2_ee', 'anc_es', 'anc_ee', 'id0', 'ijc_s1', 'sjc_s1',
                                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
                                   'inc_s1', 'inc_s2', 'inc_dif']
        df['jxn_type'] = 'A3SS'
        df['down_es'] = -1
        df['down_ee'] = -1

        return df

    def get_junction_count_array(self):
        """
        Output an array of all the junction SJC sum counts for this set of rMATS results
        :return: True
        """

        tot = pd.concat(self.rmats_list, ignore_index=True)

        junctions = [Junction(**tot.iloc[i].to_dict()) for i in range(len(tot))]

        self.sum_sjc_array = np.array([[j.sum_sjc + 1] for j in junctions])

        return True


class Junction(object):
    """
    Splice junctions (one row in an rMTAS output file), containing gene ID, and exon coordinates. There are five types
    of exon splice junctions, as from the five rMATS results files.

    """

    def __init__(self,
                 **kwargs):
        """
        Initialize a Junction object
        :param logger: logger object
        :param kwargs: keyword arguments, including the following:
        gene_id: gene ID
        gene_symbol: gene symbol
        chr: chromosome
        strand: strand
        alt1_es: alternative 1 exon start
        alt1_ee: alternative 1 exon end
        alt2_es: alternative 2 exon start
        alt2_ee: alternative 2 exon end
        anc_es: anchor exon start
        anc_ee: anchor exon end
        down_es: downstream exon start
        down_ee: downstream exon end
        jxn_type: junction type, one of 'MXE', 'SE', 'RI', 'A5SS', 'A3SS'
        ijc_s1: inclusion junction count sample 1
        sjc_s1: skipped junction count sample 1
        ijc_s2: inclusion junction count sample 2
        sjc_s2: skipped junction count sample 2

        :param logger:
        :param kwargs:
        """
        self.name = str(kwargs['id'])
        self.gene_id = kwargs['gene_id']
        self.strand = kwargs['strand']
        self.chr = kwargs['chr']
        self.anc_es = kwargs['anc_es']+1
        self.anc_ee = kwargs['anc_ee']
        self.alt1_es = kwargs['alt1_es']+1
        self.alt1_ee = kwargs['alt1_ee']
        self.alt2_es = kwargs['alt2_es']+1
        self.alt2_ee = kwargs['alt2_ee']
        self.down_es = kwargs['down_es']+1
        self.down_ee = kwargs['down_ee']
        self.junction_type = kwargs['jxn_type']
        self.gene_symbol = kwargs['gene_symbol']

        # 2022-05-21 If no symbol
        if self.gene_symbol is None:
            self.gene_symbol = f'NoGeneName_{self.gene_id}'

        self.fdr = kwargs['fdr']
        self.sjc_s1 = kwargs['sjc_s1']
        self.sjc_s2 = kwargs['sjc_s2']
        self.ijc_s1 = kwargs['ijc_s1']
        self.ijc_s2 = kwargs['ijc_s2']

        self.tx1 = None
        self.tx0 = None
        self.phase = None
        self.slice_phases = (None, None)
        self.phase_ambiguous = False
        self.num_start_codons = 0

        self.logger = logging.getLogger('jcast')


    def __repr__(self):
        """ repr """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    def __str__(self):
        """ str """
        return 'Splice junction object: ' + self.gene_id + ' ' \
               + self.junction_type + ' ' + self.gene_symbol + ' ' + self.name

    @property
    def sum_sjc(self):
        """
        Returns the sum of all SJCs in the junction. if rMATS was run with one technical replicate,
        the count field is an int, otherwise it is a list. Currently this takes the skipped junction count (SJC)
        as filtering criterion because the majority of translatable events are probably SE (skipped exon).
        Essentially this filters out alternative junctions that are very rarely skipped (high inclusion
        level of the exons) that are not likely to be translatable.

        2022-05-21: This now returns the min of SJC and IJC sums

        :return:
        """
        try:
            sum_sjc_count_sample1 = int(sum([int(x) for x in (str(self.sjc_s1).split(sep=','))]))
            sum_sjc_count_sample2 = int(sum([int(x) for x in (str(self.sjc_s2).split(sep=','))]))

        except ValueError:
            sum_sjc_count_sample1 = 0
            sum_sjc_count_sample2 = 0

        try:
            sum_ijc_count_sample1 = int(sum([int(x) for x in (str(self.ijc_s1).split(sep=','))]))
            sum_ijc_count_sample2 = int(sum([int(x) for x in (str(self.ijc_s2).split(sep=','))]))

        except ValueError:
            sum_ijc_count_sample1 = 0
            sum_ijc_count_sample2 = 0

        return min(sum_sjc_count_sample1 + sum_sjc_count_sample2, sum_ijc_count_sample1 + sum_ijc_count_sample2)

    def _get_translated_region(self,
                               gtf,
                               startsite_index: int = 0,
                               ):
        """
        Read the genome annotation (.gtf) file, find the coding sequences (CDS) that share the gene name and coordinates
        of the anchor exon (anc) supplied, then find out where the annotated translation start and end sites are. If the
        splice junction sequences extend BEYOND the translation starts and ends, trim them to avoid running into stop
        codons.

        :param gtf: genome annotation
        :param startsite_index: the index of which start site to retrieve; default is 0 (most upstream TSS).
        :return: True
        """
        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')

        #
        # 	Get the translation start and end positions
        #
        tsl = params.tsl_threshold

        # 2020-07-25 now getting the start and stop codons of all protein coding transcripts at TSL threshold
        # 2022-05-21 if there is no TSL, do not look up TSL
        # 2026-11-10 the end of translation now comes from stop codons (it was taken from start codons);
        #            see reports/2026-11-03-tx1-translation-bounds-audit.md
        coding = gtf0.query('transcript_biotype == "protein_coding"')
        if 'transcript_support_level' in coding.columns:
            coding = coding.query('transcript_support_level <= @tsl')

        starts = coding.query('feature == "start_codon"')
        stops = coding.query('feature == "stop_codon"')

        # Number of start sites:
        self.num_start_codons = len(starts['start'].drop_duplicates())

        #
        # tx0 and tx1 are the lower and upper genomic bounds of the coding region, as _trim expects.
        # The start site is chosen by startsite_index from the most upstream (longest transcript);
        # the end is the most downstream stop codon. GTF codon coordinates are inclusive, so no offset is needed.
        # Each bound stays None (no trimming on that side) if it cannot be determined.
        #
        if self.strand == '+':
            five_prime = sorted(starts['start'].drop_duplicates())
            self.tx0 = five_prime[startsite_index] if len(five_prime) > startsite_index else None
            self.tx1 = int(stops['end'].max()) if len(stops) > 0 else None

        else:
            five_prime = sorted(starts['end'].drop_duplicates(), reverse=True)
            self.tx1 = five_prime[startsite_index] if len(five_prime) > startsite_index else None
            self.tx0 = int(stops['start'].min()) if len(stops) > 0 else None

        self.logger.debug('Chosen coding region is {0} to {1}; tsl is {2}.'.format(self.tx0,
                                                                                     self.tx1,
                                                                                     tsl,
                                                                                     )
                          )
        return True

    def get_translated_phase(self, gtf):
        """
        Get the annotated translation phase of each slice from the GTF file. Call after trim_cds.

        The phase is read at the 5' end of each slice (lowest coordinate on +, highest on -), i.e. at the start
        of its first retained exon, from the protein-coding CDS rows covering that position: first those at
        the TSL threshold used for the coding bounds, then any TSL. Where transcripts disagree on the frame,
        the most frequent is used and phase_ambiguous is set. A slice whose 5' end no CDS covers gets None.

        2026-11-10 replaces a lookup that required one fixed exon (e.g. the anchor) to match a CDS row exactly
        at both ends, which failed whenever that exon had been trimmed or removed, and which assumed both slices
        start in the same exon. See reports/2026-11-10-phase-audit.md

        :param gtf:  genome annotation
        :return: None
        """

        # Subset the gtf file
        gtf0 = gtf.annot.query('gene_id == @self.gene_id')
        coding = gtf0.query('feature == "CDS" & transcript_biotype == "protein_coding"')

        tsl = params.tsl_threshold
        supported = coding
        if 'transcript_support_level' in coding.columns:
            supported = coding.query('transcript_support_level <= @tsl')

        phases, self.phase_ambiguous = [], False
        for exons in self.slice_exons():
            if not exons:
                phases.append(None)
                continue

            five_prime = min(es for es, ee in exons) if self.strand == '+' else max(ee for es, ee in exons)
            frames = self._frames_at(supported, five_prime) or self._frames_at(coding, five_prime)
            if not frames:
                phases.append(None)
                continue

            counts = collections.Counter(frames)
            self.phase_ambiguous = self.phase_ambiguous or len(counts) > 1
            phases.append(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0])

        self.slice_phases = tuple(phases)

        # a single phase only when both slices share it
        self.phase = phases[0] if phases[0] == phases[1] else None

        self.logger.info('Transcription start: {0} Transcript end: {1}'.format(self.tx0,
                                                                               self.tx1))
        self.logger.info('Retrieved phases: {0}'.format(self.slice_phases))

    def _frames_at(self, cds, position: int) -> list:
        """
        The reading frame at a genomic position, from every CDS row covering it

        GTF frame is the number of bases to skip from the 5' end of the CDS row to the first complete codon,
        so the frame at a position inside the row is offset by its distance from that end.

        :param cds: CDS rows of the gene
        :param position: genomic position, the 5' end of a slice
        :return: list of frames (0, 1 or 2), one per covering CDS row
        """

        rows = cds[(cds['start'] <= position) & (cds['end'] >= position)]
        frames = rows['frame'].astype(int)

        if self.strand == '+':
            return [int(x) for x in (rows['start'] + frames - position) % 3]
        return [int(x) for x in (frames - (rows['end'] - position)) % 3]

    @property
    def phase_label(self) -> str:
        """ the retrieved phase for logging: one value, or slice1/slice2 when the slices differ """
        if self.slice_phases[0] == self.slice_phases[1]:
            return str(self.slice_phases[0])
        return '{0}/{1}'.format(*self.slice_phases)

    def trim_cds(self, gtf):
        """
        Wrapper to trimming the exons by CDS; first retrieve the coordinates of the translation starts and ends,
        then
        """

        if self.tx0 is None:
            self._get_translated_region(gtf=gtf,
                                        startsite_index=0)

        self._trim()

        # 2026-11-10 removed the retry with the second start site when the trimmed anchor was not a CDS exon.
        # It re-trimmed to a narrower, gene-wide bound, and once exons outside the bounds are removed this
        # deleted coding exons of the longest isoform; see reports/2026-11-10-utr-exon-trimming-audit.md


    def _trim(self):
        """
        Core trim function that shaves off the UTRs
        :return: True

        Trims every exon of the junction to the coding region [tx0, tx1]: an exon containing a bound is cut at
        that bound, and an exon lying entirely outside the coding region is removed (coordinates set to -1).
        The bounds are genomic, so the same rule applies to every junction type on either strand.
        A bound of None means it could not be determined, and that side is not trimmed.

        2026-11-10 replaces the per-type rules, which never removed whole UTR exons and on the - strand trimmed
        A5SS/A3SS exons at the wrong bound; see reports/2026-11-10-utr-exon-trimming-audit.md
        """

        for exon in ['anc', 'alt1', 'alt2', 'down']:
            es, ee = getattr(self, exon + '_es'), getattr(self, exon + '_ee')

            # skip exons that are absent for this junction type or already removed
            if es <= 0 or ee <= 0:
                continue

            if (self.tx0 is not None and ee < self.tx0) or (self.tx1 is not None and es > self.tx1):
                es, ee = -1, -1

            else:
                if self.tx0 is not None and es < self.tx0 <= ee:
                    es = self.tx0
                if self.tx1 is not None and es <= self.tx1 < ee:
                    ee = self.tx1

            setattr(self, exon + '_es', es)
            setattr(self, exon + '_ee', ee)

        return True

    @property
    def event_outside_cds(self) -> bool:
        """
        Whether trimming has left both slices with the same exons, i.e. the splice event lies entirely
        outside the coding region and does not change the protein. Only meaningful after trim_cds.
        The slices are composed as in Sequence.make_slice_localgenome.

        :return: True if the two slices are identical
        """

        slice1, slice2 = self.slice_exons()
        return slice1 == slice2

    def slice_exons(self):
        """
        The retained exons of each slice as sorted (start, end) tuples, composed as in
        Sequence.make_slice_localgenome.

        :return: (slice 1 exons, slice 2 exons)
        """

        if self.junction_type in ['MXE', 'SE']:
            slice1, slice2 = ['anc', 'alt1', 'down'], ['anc', 'alt2', 'down']
        elif self.junction_type == 'RI':
            slice1, slice2 = ['alt1'], ['anc', 'down']
        else:
            slice1, slice2 = ['anc', 'alt1'], ['anc', 'alt2']

        def exons(names):
            return sorted((getattr(self, n + '_es'), getattr(self, n + '_ee')) for n in names
                          if getattr(self, n + '_es') > 0 and getattr(self, n + '_ee') > 0)

        return exons(slice1), exons(slice2)
