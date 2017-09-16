#
#   Classes that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids
#


def get_nuc(species, chr, es, ee):
    """
    :param species: Species, e.g., 'mouse
    :param chr: Chromosome, e.g., 1
    :param es: Exon start, e.g., 10000000
    :param ee: Exon end, e.g., 10000100
    :return:

    Function to get nucleotide sequences by REST API
    (anchor, alternative-1, alternative-2, downstream)

    To cut down time, it will try to cache the sequence coordinates.

    >>> get_nuc('mouse', 1, 10000000, 10000050)
    'GTTTTCAATGCAGGAAATGCAATTGTTCTGTAGGTACAAGTGGGTCAGATT'

    """

    import requests as rq
    import sqlite3 as sq

    if es <= 0 or ee <= 0:
        print("Skipping empty exon.")
        return ''

    cache = species + '-' + str(chr) + '-' + str(es) + '-' + str(ee)


    con = sq.connect('seq-cache.db')
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS sequences(pk INTEGER PRIMARY KEY, id TEXT, seq TEXT)''')

    try:
        cur.execute('''SELECT id, seq FROM sequences WHERE id=:cache''',
                    {'cache': cache})
        sequence = cur.fetchone()

        if len(sequence) > 0 :
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

    bp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    for n in nt:
        try:
            comp += bp[n]
        except:
            comp += 'X'

    return comp


def make_pep(nt, strand, phase):
    """
    :param nt: Nucleotide sequence
    :param strand: Strand (+ or -)
    :param phase: Translation frame (0, 1, of 2)
    :return: Amino acid sequence

    Function to translate a nucleotide sequence into peptide,
    taking into account the nucleotide sequence, strand, and phase.

        >>> make_pep('ATTTTGCTT', '+', 0)
        Translating on the + strand with phase 0
        'ILL'

    Phase shifts the starting position according to Ensembl convention.

        >>> make_pep('ATTTTGCTT', '+', 1)
        Translating on the + strand with phase 1
        'FA'

        >>> make_pep('ATTTTGCTT', '+', 2)
        Translating on the + strand with phase 2
        'FC'

    If the phase of a sequence is on the negative strand, get
    the complementary sequence then translate

        >>> make_pep('ATTTTGCTT', '-', 0)
        Translating on the - strand with phase 0
        'KQN'
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

    if strand == '-':
        nt = get_complementary(nt)

    #
    # Get the starting position to translate, based on the phase
    #
    pos0 = (3 - phase) % 3

    print('Translating on the ' + strand + ' strand with phase ' + str(phase))

    pep = ''

    for i in range(pos0, len(nt) - 2, 3):

        aa = code[nt[i:i + 3]]

        if aa == 'X':
            print('Stop codon encountered')
            pep = ''
            break

        pep += aa

    return pep



class Sequence(object):

    def __init__(self, junction):
        """
        :type junction: object
        :param junction: The splice junction object

        Mostly copying properties of the junction object (already trimmed)

        """
        self.anc_ee = junction.anc_ee
        self.anc_es = junction.anc_es
        self.alt1_es = junction.alt1_es
        self.alt1_ee = junction.alt1_ee
        self.alt2_es = junction.alt2_es
        self.alt2_ee = junction.alt2_ee
        self.down_es = junction.down_es
        self.down_ee = junction.down_ee
        self.species = junction.species
        self.gene_id = junction.gene_id
        self.junction_type = junction.junction_type
        self.chr = junction.chr
        self.phase = junction.phase
        self.strand = junction.strand
        self.slice1_nt = ''
        self.slice2_nt = ''
        self.slice1_aa = ''
        self.slice2_aa = ''
        self.gene_symbol = junction.gene_symbol
        self.name = junction.name

    def make_slice(self):

        anc_nt = get_nuc(self.species, self.chr, self.anc_es, self.anc_ee)
        alt1_nt = get_nuc(self.species, self.chr, self.alt1_es, self.alt1_ee)
        alt2_nt = get_nuc(self.species, self.chr, self.alt2_es, self.alt2_ee)
        down_nt = get_nuc(self.species, self.chr, self.down_es, self.down_ee)

        self.slice1_nt = anc_nt + alt1_nt + down_nt
        self.slice2_nt = anc_nt + alt2_nt + down_nt

        print(self.slice1_nt)
        print(self.slice2_nt)

    def translate(self):
        """

        :return: True
        """


        if self.phase in [0, 1, 2]:
            self.slice1_aa = make_pep(self.slice1_nt, self.strand, self.phase)
            self.slice2_aa = make_pep(self.slice2_nt, self.strand, self.phase)
            print("Used Retrieved Phase")

        else:
            for i in range(3):
                self.slice1_aa = make_pep(self.slice1_nt, self.strand, i)
                self.slice2_aa = make_pep(self.slice2_nt, self.strand, i)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    break
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        print(self.slice1_aa)
        print(self.slice2_aa)

        return True

    def write_to_fasta(self, output):
        """
        :param output: File name of the .fasta output.
        :return: True

        Create the /out directory if it does not exist, then write the translated splice junctions into the .fasta file
        if both slices are translated.

        """


        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import os
        import os.path

        if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:

            fa1 = SeqRecord(Seq(self.slice1_aa, IUPAC.protein),
                            id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-1-' +
                                self.name + '-' + str(self.phase) + self.strand),
                            name=self.gene_symbol,
                            description='Slice 1')
            fa2 = SeqRecord(Seq(self.slice2_aa, IUPAC.protein),
                            id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-2-' +
                                self.name + '-' + str(self.phase) + self.strand),
                            name=self.gene_symbol,
                            description='Slice 2')

            from Bio import SeqIO

            os.makedirs('out', exist_ok=True)
            o = os.path.join('out', output + '.fasta')

            output_handle = open(o, 'a')
            SeqIO.write([fa1, fa2], output_handle, 'fasta')
            output_handle.close()


        return True



#
# For doctest
#
if __name__ == '__main__':
    import doctest
    doctest.testmod()