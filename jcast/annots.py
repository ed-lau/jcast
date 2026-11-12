# -*- coding: utf-8 -*-

""" Genome annotations. """
import os
import glob
import hashlib
import logging
import gzip
import tempfile

import gtfparse
import pandas as pd
from Bio import SeqIO
#from Bio.Alphabet import generic_dna

from jcast import cache


class ReadAnnotations(object):
    """
    Class holds the name and location of the GTF file, plus a pandas dataframe
    """

    def __init__(self,
                 logger: logging.Logger,
                 path: str, ):
        """
        :parrm logger: logger object
        :param path: Path to GTF file
        """
        self.logger = logger
        self.path = os.fspath(path)

        self.annot = pd.DataFrame()
        self.by_gene = {}
        self.gene_name_counts = {}

        self.read_gtf()


    def read_gtf(self) -> None:
        """
        Read gtf file based on the location and name supplied, from the cache if it was read before

        The table is the GTF as gtfparse reads it, written as a tab-separated file and read back. That round trip
        sets the column types the rest of JCAST relies on: exon_number is numeric (exons are sorted by it),
        transcript_support_level 'NA' is missing, frame is an integer. The cache keeps that table as a pickle in
        the cache directory, keyed by the GTF's path, size and modification time and the pandas version, and is
        renamed into place once written, so a run never reads a partial cache.

        2026-11-12 the cache was <gtf>.cached next to the GTF, keyed by its name only and written in place: stale
        after the GTF changed, partial while another run wrote it, impossible in a read-only directory. See
        reports/2026-11-12-thresholds-runs-canonical-audit.md §3.5

        :return:
        """

        self.logger.info("Reading GTF file. This could take a minute.")

        cache_dir = cache.directory('gtf')
        cached_gtf = _cache_path(self.path, cache_dir) if cache_dir else None

        self.annot = None
        if cached_gtf and os.path.isfile(cached_gtf):
            try:
                self.annot = pd.read_pickle(cached_gtf)
                self.logger.info(f'Read the GTF from the cache {cached_gtf}.')
            except Exception as e:
                self.logger.warning(f'Cannot read the GTF cache {cached_gtf} ({e}); reading the GTF again.')

        if self.annot is None:
            self.annot = _read_gtf_table(self.path, cache_dir)
            if cached_gtf:
                _write_cache(self.annot, cached_gtf, self.logger)

        # 2020-11-06 if Gencode GTF is used, column is transcript_type rather than transcript_biotype.
        # This will make the transcript_biotype column manually

        if 'transcript_biotype' not in self.annot.columns:
            self.annot['transcript_biotype'] = self.annot['transcript_type']

        # every lookup filters by gene first: index the rows by gene once instead of scanning the whole GTF each
        # time (2026-11-12; reports/2026-11-12-parallel-audit.md)
        self.by_gene = {g: rows for g, rows in self.annot.groupby('gene_id', sort=False)}

        # gene name -> number of genes with it: names are not unique (e.g. two GENCODE v50 genes are PINX1), so a name
        # identifies a gene's UniProt entry only if no other gene has it (reports/2026-11-12-uniprot-mapping-audit.md)
        if 'gene_name' in self.annot.columns:
            genes = self.annot[['gene_id', 'gene_name']].drop_duplicates('gene_id')
            self.gene_name_counts = genes['gene_name'].value_counts().to_dict()

        return None


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def _cache_path(path: str, cache_dir: str) -> str:
    """ <cache dir>/<GTF name>.<digest of its path>.<digest of its path, size, mtime and the pandas version>.pkl """
    path = os.path.abspath(path)
    stat = os.stat(path)
    key = f'{path}\t{stat.st_size}\t{stat.st_mtime_ns}\tpandas {pd.__version__}'
    return os.path.join(cache_dir, f'{os.path.basename(path)}.{_digest(path)}.{_digest(key)}.pkl')


def _read_gtf_table(path: str, tmp_dir: str = None) -> pd.DataFrame:
    """ The GTF through gtfparse, then written as a tab-separated file and read back, which sets the column types """
    annot = gtfparse.read_gtf(path)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as d:
        tsv = os.path.join(d, 'gtf.tsv')
        annot.to_csv(tsv, encoding='utf-8', sep='\t', index=False)
        return pd.read_table(tsv, sep='\t', low_memory=False)


def _write_cache(annot: pd.DataFrame, cached_gtf: str, logger: logging.Logger) -> None:
    """ Write the cache under a temporary name and rename it into place; then remove this GTF's older caches """
    tmp = f'{cached_gtf}.{os.getpid()}.tmp'
    try:
        annot.to_pickle(tmp)
        os.replace(tmp, cached_gtf)
    except OSError as e:
        logger.warning(f'Cannot write the GTF cache {cached_gtf} ({e}); the next run reads the GTF again.')
        if os.path.exists(tmp):
            os.remove(tmp)
        return
    logger.info(f'Wrote the GTF cache {cached_gtf}.')

    same_gtf = glob.escape(cached_gtf.rsplit('.', 2)[0]) + '.*.pkl'
    for old in glob.glob(same_gtf):
        if old != cached_gtf:
            try:
                os.remove(old)
            except OSError:
                pass


def gene_name_is_unique(gtf, gene_name) -> bool:
    """ False if more than one gene in the GTF has this name; True for a name the GTF does not have """
    return getattr(gtf, 'gene_name_counts', {}).get(gene_name, 0) <= 1


def gene_rows(gtf, gene_id) -> pd.DataFrame:
    """
    The GTF rows of one gene, in file order: from the per-gene index if the annotation has one, else by filtering.

    :param gtf: ReadAnnotations, or any object with an annot data frame
    :param gene_id: gene id
    :return: data frame
    """
    by_gene = getattr(gtf, 'by_gene', None)
    if by_gene:
        rows = by_gene.get(gene_id)
        return gtf.annot.iloc[:0] if rows is None else rows
    return gtf.annot.query('gene_id == @gene_id')


class ReadGenome(object):

    def __init__(self,
                 logger: logging.Logger,
                 path: str, ):
        """
        :param logger: logger object
        :param path: Path to genome fasta file, optionally gzipped

        """


        self.path = os.fspath(path)
        self.logger = logger

        self.genome = None
        self._read_fasta()


    def _read_fasta(self):
        """
        Reading genome. Just a wrapper for Biopython SeqIO with gzip if needed.

        :return:
        """

        if self.path.endswith('.gz'):
            with gzip.open(self.path, 'rt') as f:
                self.genome = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                self.logger.info('Read from zipped genome. \n\n')

        else:
            with open(self.path, 'rt') as f:
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
