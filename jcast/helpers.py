# -*- coding: utf-8 -*-

""" Helper functions that translate nucleotides to amino acids, etc. """

import os.path
import sqlite3 as sq

import requests as rq

from Bio import SeqIO
# from Bio.Alphabet import IUPAC

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

    return genome_index[chromosome][es-1:ee]


def get_rest_nuc(species,
                 chr,
                 es,
                 ee,
                 ):
    """
    :param species: Species, e.g., 'mouse
    :param chr: Chromosome, e.g., 1
    :param es: Exon start, e.g., 10000000
    :param ee: Exon end, e.g., 10000100
    :return:

    Function to get nucleotide sequences by REST API
    (anchor, alternative-1, alternative-2, downstream)

    If the coordinates are minus 1 such as in the case when a certain segment is not present in a particular splice
    type, this should return an empty string

    To reduce loading time, it will try to cache the sequence coordinates in a local SQLite DB.

    To-do: note that if the network retrieval is affected by some redirecting so we are retrieving a HTML 200 result
    that does not contain sequence, this could corrupt the database and require manual removal of the cached sequence.

    >>> get_nuc('mouse', 1, 10000000, 10000050)
    'GTTTTCAATGCAGGAAATGCAATTGTTCTGTAGGTACAAGTGGGTCAGATT'

    >>> get_nuc('mouse', 1, -1, -1)
    ''

    """

    if es <= 0 or ee <= 0:
        return ''

    cache = species + '-' + str(chr) + '-' + str(es) + '-' + str(ee)

    con = sq.connect('seq-cache.db')
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')

    try:
        cur.execute('''SELECT id, seq FROM sequences WHERE id=:cache''',
                    {'cache': cache})
        sequence = cur.fetchone()

        if len(sequence) > 0:
            nuc = sequence[1]
            print('Locally cached sequence retrieved')
            con.close()
            return nuc

        else:
            print("Sequence not yet cached locally. 0")

    except:
        print("Sequence not yet cached locally. 1")

    server = "http://rest.ensembl.org"
    ext = "/sequence/region/" + species + "/" + str(chr) + ":" + str(es) + ".." + str(ee) + ":1?"

    print(server+ext)

    try:
        ret = rq.get(server + ext, headers={"Content-Type": "text/plain"})

        if ret.status_code == 200:
            nuc = ret.text
            cur.execute('''INSERT INTO sequences(id, seq) VALUES(:id, :seq)''', {'id': cache, 'seq': nuc})
            print("Sequence retrieved from Uniprot and written into local cache.")
            con.commit()
            con.close()

        else:
            nuc = ''

    except:
        nuc = ''

    return nuc


def get_complementary(nt: str,
                      ) -> str:
    """
    :param nt: Nucleotide string
    :return: Complementary sequence

    Function to get the complementary coding sequence on the reverse (-) DNA
    strand of the coordinates being give. First reverses the input nucleotide
    sequence, then get base-pairing nucleotide

    >>> get_complementary('ATGCAA')
    'TTGCAT'

    """

    nt = nt[::-1]
    comp = ''

    bp = {'A': 'T',
          'T': 'A',
          'G': 'C',
          'C': 'G'}

    for n in nt:
        try:
            comp += bp[n]
        except KeyError:
            comp += 'X'

    return comp


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

    """

    #
    # dictionary for genetic code
    #
    code = constants.genetic_code

    #
    # if rMATS says the strand is negative, get the complementary sequence, which is the Watson-Crick pairing of the
    # original sequence in anti-parallel direction (5'-to-3' becomes 3'-to-5')
    #
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

        aa = code[nt[i:i + 3]]

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


def write_seqrecord_to_fasta(seqrecord,
                             output,
                             suffix,
                             ):
    """
    Write Biopython SeqRecord to the fasta file after checking whether the SeqRecord is already inside the file.

    :param seqrecord:   Biopython SeqRecord object
    :param output:
    :param suffix:
    :return:
    """

    outfile = os.path.join(output, 'psq_' + suffix + '.fasta')

    # if the file already exists, open it and amend that record.
    existing_records = []
    if os.path.exists(outfile):
        for existing_record in SeqIO.parse(outfile, 'fasta'):
            existing_records.append(existing_record)

    # Test if the slice is already in the fasta, then do not write the new sequence into the fasta file.
    for read_record in existing_records:
        if read_record.seq == seqrecord.seq:
            return True

    with open(outfile, 'a') as output_handle:
        SeqIO.write(seqrecord, output_handle, 'fasta')

    return True
