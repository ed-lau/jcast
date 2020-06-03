# -*- coding: utf-8 -*-

""" Helper functions that translate nucleotides to amino acids, etc. """

import requests as rq
import sqlite3 as sq
import os.path

from Bio import Seq
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, DNAAlphabet


def get_local_nuc(genome_index, chromosome, es, ee):
    """

    :param genome_index: genome index
    :param chr: chromosome
    :param es: exon start
    :param ee: exon end
    :return: sequence

    >>> idx = SeqIO.index("/Volumes/Gingko_Data/Data/Exosome_smRNAseq/ReferenceGenomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa", 'fasta')
    >>> len(get_local_nuc(idx, "chr17", 44091700, 44091831))
    132
    """

    # Skip empty exons
    if es <= 0 or ee <= 0:
        #print("Skipping empty exon.")
        return Seq('', DNAAlphabet())

    # Get numeric portion of chromosome
    chromosome = chromosome[3:]

    return genome_index[chromosome][es-1:ee]


def get_nuc(species, chr, es, ee):
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
        #print("Skipping empty exon.")
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


def get_complementary(nt):
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
        except:
            comp += 'X'

    return comp


def make_pep(nt, strand, phase, terminate=True):
    """

    :param nt: Nucleotide sequence
    :param strand: Strand (+ or -)
    :param phase: Translation frame (0, 1, of 2)
    :param terminate: T/F End if running into premature termination codons

    :return: Amino acid sequence

    Function to translate a nucleotide sequence into peptide,
    taking into account the nucleotide sequence, strand, and phase.

    At the moment this will stop when you run into a stop codon (X)!
    To-do: create a variant where the sequence till before the stop codon is returned.

        >>> rna_sequence = 'CTTAAATCCAGCCACTGCTCCAGACTGCAAATCGGATTCAATATTTCCTGGAT\\
        CCAGGTAGGCAATGCTCATAAGAAAACCTGGTCCGGTGAAAGCCCAGAGTTTACGAAAGCTAAAACAAGAGT\\
        ACTCCTCCTCAGGAATGGAGATCTTCTCATTAAAGTAAGTGGCGAAGTACTCCTCTGAGTCCCCAGGGGACATCT\\
        GACATCTTCTGTTCAGGACCCAGCACCATGGTGGATACCTGAGTGGCTGAGTTCTTAGAATATGATT'
        >>> make_pep(rna_sequence, '+', 0, False)
        Translating on the + strand with phase 0
        Stop codon encountered
        'LKSSHCSRLQIGFNISWIQVGNAHKKTWSGESPEFTKAKTRVLLLRNGDLLIKVSGEVLL'

    Phase shifts the starting position according to Ensembl convention.

        >>> make_pep(rna_sequence, '+', 1, False)
        Translating on the + strand with phase 1
        Stop codon encountered
        'LNPATAPDCKSDSIFPGSR'


        >>> make_pep(rna_sequence, '+', 2, True)
        Translating on the + strand with phase 2
        Stop codon encountered
        ''

    If the phase of a sequence is on the negative strand, get
    the complementary sequence then translate

        >>> make_pep(rna_sequence, '-', 1, False)
        Translating on the - strand with phase 1
        Stop codon encountered
        'IIF'

    """

    #
    # Dictionary for genetic code
    #
    code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    #
    # If rMATS says the strand is negative, get the complementary sequence, which is the Watson-Crick pairing of the
    # original sequence in anti-parallel direction (5'-to-3' becomes 3'-to-5')
    #
    if strand == '-':
        nt = get_complementary(nt)

    #
    # Get the starting position to translate, based on the phase. Will have to check this carefully.
    #
    pos0 = phase # (3 - phase) % 3

    #print('Translating on the ' + strand + ' strand with phase ' + str(phase))

    pep = ''

    has_ptc = False
    #
    # After finding the phase, loop through every 3 nucleotides from the start
    # then use the dictionary to find the amino acid sequence.
    #
    for i in range(pos0, len(nt) - 2, 3):

        aa = code[nt[i:i + 3]]

        if aa == 'X':
            #print('Stop codon encountered')
            has_ptc = True

            # If the terminate flag is set to true, return an empty sequence
            if terminate:
                pep = ''
                

            # If terminate is set to false, return present sequence
            break

        pep += aa

    return pep #[pep, has_ptc] 2018-09-08 not sure why this returns a list. Forgot what I wanted to do



def write_seqrecord_to_fasta(seqrecord, output, suffix):
    """
    Write Biopython SeqRecord to the fasta file after checking whether the SeqRecord is already inside the file.

    :param seqrecord:   Biopython SeqRecord object
    :param sequence:    Splice Sequence object
    :param output:
    :param suffix:
    :return:
    """


    o = os.path.join(output, 'psq_' + suffix + '.fasta')
    #print(o)

    # If the file already exists, open it and amend that record.
    existing_records = []
    if os.path.exists(o):
        for existing_record in SeqIO.parse(o, 'fasta', IUPAC.extended_protein):
            existing_records.append(existing_record)

    # Test if the slice is already in the fasta, then do not write the new sequence into the fasta file.
    for read_record in existing_records:
        if read_record.seq == seqrecord.seq:
            #print(seqrecord.seq)
            #print("Already in fasta file - we will consider modifying the fasta entry name later on to reflect this")
            return True

    output_handle = open(o, 'a')
    SeqIO.write(seqrecord, output_handle, 'fasta')
    output_handle.close()

    return True
