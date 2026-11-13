# -*- coding: utf-8 -*-

""" Helper functions that translate nucleotides to amino acids, etc. """

import os.path

from Bio import SeqIO
from Bio.Seq import reverse_complement

from jcast import constants


def get_local_nuc(genome_index,
                  chromosome,
                  es: int,
                  ee: int,
                  ):
    """

    :param genome_index: genome index
    :param chromosome: chromosome
    :param es: exon start
    :param ee: exon end
    :return: sequence

    """

    # try Gencode chromosome nomenclature scheme ('chr1', etc.)
    try:
        genome_index[chromosome][0:1]
    except KeyError:
        # get numeric portion of chromosome if using Ensembl (removing chr)
        chromosome = chromosome[3:]  # skip 'chr'

    # Skip empty exons
    if es <= 0 or ee <= 0:
        return genome_index[chromosome][0:0]  # 2020-07-21 return SeqRecord rather than Seq

    # soft-masked (lowercase) bases are ordinary sequence (reports/2026-11-11-masked-bases-audit.md)
    return genome_index[chromosome][es-1:ee].upper()


def get_complementary(nt: str,
                      ) -> str:
    """
    :param nt: Nucleotide string
    :return: Complementary sequence

    Function to get the complementary coding sequence on the reverse (-) DNA
    strand of the coordinates being give. Ambiguous (IUPAC) bases are complemented
    and case is kept.

    >>> get_complementary('ATGCAA')
    'TTGCAT'

    >>> get_complementary('acgR')
    'Ycgt'

    """

    return reverse_complement(str(nt))


def make_pep(nt: str,
             strand: str,
             phase: int,
             terminate: bool = True,
             ) -> str:
    """

    :param nt: Nucleotide sequence
    :param strand: Strand (+ or -)
    :param phase: Translation frame (0, 1, of 2)
    :param terminate: T/F End if running into premature termination codons

    :return: Amino acid sequence

    Function to translate a nucleotide sequence into peptide,
    taking into account the nucleotide sequence, strand, and phase.

    If terminate flag is True (default), return empty string is there is a stop codon.
    Otherwise return peptides up till codon.

        >>> rna_sequence = 'CTTAAATCCAGCCACTGCTCCAGACTGCAAATCGGATTCAATATTTCCTGGAT\\
        CCAGGTAGGCAATGCTCATAAGAAAACCTGGTCCGGTGAAAGCCCAGAGTTTACGAAAGCTAAAACAAGAGT\\
        ACTCCTCCTCAGGAATGGAGATCTTCTCATTAAAGTAAGTGGCGAAGTACTCCTCTGAGTCCCCAGGGGACATCT\\
        GACATCTTCTGTTCAGGACCCAGCACCATGGTGGATACCTGAGTGGCTGAGTTCTTAGAATATGATT'
        >>> make_pep(rna_sequence, '+', 0, False)
        'LKSSHCSRLQIGFNISWIQVGNAHKKTWSGESPEFTKAKTRVLLLRNGDLLIKVSGEVLL'

    Phase shifts the starting position according to Ensembl convention.

        >>> make_pep(rna_sequence, '+', 1, False)
        'LNPATAPDCKSDSIFPGSR'

        >>> make_pep(rna_sequence, '+', 2, True)
        ''

    If the phase of a sequence is on the negative strand, get
    the complementary sequence then translate

        >>> make_pep(rna_sequence, '-', 1, False)
        'IIF'

    Lowercase (soft-masked) bases are translated as uppercase. A codon with any other base (N or another
    IUPAC code) cannot be translated, so the whole sequence is untranslatable

        >>> make_pep('atgaaatag', '+', 0)
        'MK'

        >>> make_pep('ATGRAATAG', '+', 0, False)
        ''

    """

    #
    # dictionary for genetic code
    #
    code = constants.genetic_code

    #
    # if rMATS says the strand is negative, get the complementary sequence, which is the Watson-Crick pairing of the
    # original sequence in anti-parallel direction (5'-to-3' becomes 3'-to-5')
    #
    nt = str(nt).upper()
    if strand == '-':
        nt = get_complementary(nt)

    #
    # get the starting position to translate, based on the phase.
    #
    pos0 = phase  # (3 - phase) % 3
    pep = ''

    #
    # after finding the phase, loop through every 3 nucleotides from the start
    # then use the dictionary to find the amino acid sequence.
    #
    for i in range(pos0, len(nt) - 2, 3):

        aa = code.get(nt[i:i + 3])

        # a codon with an ambiguous base: untranslatable
        if aa is None:
            return ''

        if aa == 'X':

            # if the terminate flag is set to true and this is not the last codon, return an empty sequence
            if terminate:
                if (len(nt)-2) - i < 3:
                    return pep
                else:
                    return ''

            # If terminate is set to false, return present sequence
            else:
                return pep

        pep += aa

    return pep


# output file -> (its size after our last write, the sequences in it), so a write need not re-read the file
_written = {}

# In a worker process (-p > 1), records are collected here instead of written; the main process writes them in
# junction order (main._translate_in_worker)
_collected = None


def write_seqrecord_to_fasta(seqrecord,
                             output,
                             suffix,
                             ):
    """
    Write Biopython SeqRecord to the fasta file unless a record with the same sequence is already in it.

    The sequences in each file are kept in memory. The file is read again only if its size differs from that after
    our last write (first use, a file from an earlier run, or one changed by something else). 2026-11-12: it used to
    be re-read before every write; see reports/2026-11-12-parallel-audit.md

    :param seqrecord:   Biopython SeqRecord object
    :param output:
    :param suffix:
    :return: True if the record was written (or collected, see _collected), False if its sequence was already in the
        file
    """

    if _collected is not None:
        _collected.append((output, suffix, seqrecord))
        return True

    outfile = os.path.join(output, 'psq_' + suffix + '.fasta')

    size = os.path.getsize(outfile) if os.path.exists(outfile) else None
    cached = _written.get(outfile)
    if cached is not None and cached[0] == size:
        sequences = cached[1]
    else:
        sequences = {str(r.seq) for r in SeqIO.parse(outfile, 'fasta')} if size is not None else set()
        _written[outfile] = (size, sequences)

    # Test if the slice is already in the fasta, then do not write the new sequence into the fasta file.
    if str(seqrecord.seq) in sequences:
        return False

    with open(outfile, 'a') as output_handle:
        SeqIO.write(seqrecord, output_handle, 'fasta')
    sequences.add(str(seqrecord.seq))
    _written[outfile] = (os.path.getsize(outfile), sequences)

    return True
