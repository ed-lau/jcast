#
#   Classes that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids
#



def find_in_fasta(slice, fasta):
    """
    Given a translated junction sequence, look for the fasta entry that overlaps with it, then return the entry
    and the coordinates. This will be used to extend said junction sequence to encompass the entire protein sequence.

    Instead of reading from the fasta file, it should just read fetch the Uniprot directly via API.


    :param seq:
    :param fasta:
    :return:
    """
    from Bio import SeqIO
    import request as rq
    from io import StringIO
    import sys

    server = 'https://www.ebi.ac.uk'
    ext = '/proteins/api/proteins/Ensembl:' + sequence.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

    print(server + ext)

    ret = rq.get(server + ext, headers={"Accept": "text/x-fasta"})

    if not ret.ok:
        ret.raise_for_status()
        sys.exit()

    # The UniProt API retrieves a retrieval object, with a text field inside ret.text
    # Since Biopython SeqIO only works with file, use io.StringIO to turn the string into a file for parsing.
    for record in SeqIO.parse(StringIO(ret.text), 'fasta'):

        # Find out where the first 10 amino acids meets the UniProt canonical sequences..
        merge_start1 = record.seq.find(sequence.slice1_aa[:10])
        merge_end1 = record.seq.find(sequence.slice1_aa[-10:])

        merge_start2 = record.seq.find(sequence.slice2_aa[:10])
        merge_end2 = record.seq.find(sequence.slice2_aa[-10:])

        if (merge_start1 != -1 and merge_end1 != -1) or (merge_start2 != -1 and merge_start2 != -1):

            # Write the UniProt canonical first
            canonical = record

            # If the slice is not the same as the UniProt canonical, then also write it.
            if not record.seq.find(sequence.slice1_aa) == -1:
                record1 = record[:merge_start1] + sequence.slice1_aa + record[merge_end1+10:]

            else:
                # Change name of canonical to reflect that it is also slice 1.

            # If the slice is not the same as the UniProt canonical, then also write it.
                if not record.seq.find(sequence.slice2_aa) == -1:
                    record1 = record[:merge_start2] + sequence.slice2_aa + record[merge_end2 + 10:]





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

            print(record.id)
            print(record.seq[:merge_start1] + sequence.slice1_aa + record.seq[merge_end1+10:])
            print(record.seq[:merge_start2] + sequence.slice2_aa + record.seq[merge_end2+10:])

    print(next(fasta).seq)




    # for edvelopment
    # slice = sequence.slice1_aa
    fasta_dct = {}
    #for record in SeqIO.parse('data/fasta/20170918_Mm_Sp_16915.fasta', 'fasta'):
    for record in SeqIO.parse('data/fasta/20170918_Hs_Sp_20205.fasta', 'fasta'):
        #fasta_dct[record.id] = record.seq # id, name, seq, description
        merge_start1 = record.seq.find(sequence.slice1_aa[:10])
        merge_end1 = record.seq.find(sequence.slice1_aa[-10:])

        merge_start2 = record.seq.find(sequence.slice2_aa[:10])
        merge_end2 = record.seq.find(sequence.slice2_aa[-10:])

        # see if you can find a FASTA entry that shares the first 10 amino acids of the slice, and the last 10
        if (merge_start1 != -1 and merge_end1 != -1) or (merge_start2 != -1 and merge_start2 != -1):
            print(record.id)
            print(record.seq[:merge_start1] + sequence.slice1_aa + record.seq[merge_end1+10:])
            print(record.seq[:merge_start2] + sequence.slice2_aa + record.seq[merge_end2+10:])

        # We have to make it so that if the canonical is neither slice 1 or slice 2 (e.g., because it might contain
        # both of the mutually exclusive exons, then the fasta should contain all three.

        # Also we still need to combine all identical entries in the fasta file.

    #slice1 = 'MLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEAAAVGAVEASFKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDAVLNAWAEDVDLRVNLAMDV'
    # fasta_dct['sp|P52480|KPYM_MOUSE'].find(slice)
    return True


class Sequence(object):

    def __init__(self, junction):
        """
        :type junction: object
        :param junction: The splice junction object

        Mostly copying properties of the junction object (already trimmed) to start a new sequence object. This
        sequence object will be used to make the nucleotide slices and translate into protein sequences.

        Change this so that sequence inherits junction class directly (have to think it through)
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
        self.frameshift = False
        self.gene_symbol = junction.gene_symbol
        self.name = junction.name

    def __str__(self):
        return "sequence object" + self.name

    def set_frameshift_to_true(self):
        """
        Mark that there is a frameshift

        :return:
        """
        self.frameshift = True

    def make_slice(self):

        import helpers as h

        anc_nt = h.get_nuc(self.species, self.chr, self.anc_es, self.anc_ee)
        alt1_nt = h.get_nuc(self.species, self.chr, self.alt1_es, self.alt1_ee)
        alt2_nt = h.get_nuc(self.species, self.chr, self.alt2_es, self.alt2_ee)
        down_nt = h.get_nuc(self.species, self.chr, self.down_es, self.down_ee)

        self.slice1_nt = anc_nt + alt1_nt + down_nt
        self.slice2_nt = anc_nt + alt2_nt + down_nt

        print(self.slice1_nt)
        print(self.slice2_nt)

    def translate(self):
        """

        :return: True
        """

        import helpers as h


        if self.phase in [0, 1, 2]:
            self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, self.phase)
            self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, self.phase)
            print("Used Retrieved Phase")

        else:
            for i in range(3):
                self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, i)
                self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, i)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    break
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        print(self.slice1_aa)
        print(self.slice2_aa)

        return True


    def write_to_fasta(self, output, suffix):
        """


        :param output: File name of the .fasta output.
        :return: True

        Create the /out directory if it does not exist, then write the translated splice junctions into the .fasta file

        """

        from Bio.Seq import Seq
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import os.path

        # Note to self, this condition should probably be moved to main to make this function more reusable.
        # So the idea is to write out different files depending on whether translation was successful


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

        os.makedirs('out', exist_ok=True)
        o = os.path.join('out', output + suffix + '.fasta')

        output_handle = open(o, 'a')
        SeqIO.write([fa1, fa2], output_handle, 'fasta')
        output_handle.close()


        return True

    def write_to_fasta_old(self, output, suffix):
        """

        DEPRECATED - see write_to_fasta above

        :param output: File name of the .fasta output.
        :return: True

        Create the /out directory if it does not exist, then write the translated splice junctions into the .fasta file
        if both slices are translated.

        """

        from Bio.Seq import Seq
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import os.path

        # Note to self, this condition should probably be moved to main to make this function more reusable.
        # So the idea is to write out different files depending on whether translation was successful
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

            os.makedirs('out', exist_ok=True)
            o = os.path.join('out', output + suffix + '.fasta')

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