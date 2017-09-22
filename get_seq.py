#
#   Classes that concern sequences - retrieving and cacheing nucleotide sequences, translating into amino acids
#


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

    # def __str__(self):
    #     return "sequence object" + self.name

    def set_frameshift_to_true(self):
        """
        A setter to mark that there is a frameshift; used to determine whether the slice should be tier 1 or tier 2

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
        This is the translate function for tier 1 peptide. Calls make_pep with the to-be-implemented "forced" flag
        which will terminate when it runs into a stop codon. Later on we should make a force_translate that does
        three frame translation and just either return all frames or return the longest.

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

        Deprecated, see the new monstrous "extend and write" function

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


        fa1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                        id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-1-' +
                            self.name + '-' + str(self.phase) + self.strand),
                        name=self.gene_symbol,
                        description='Slice 1')
        fa2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
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

    def extend_and_write(self, species, output, suffix):
        """
        Given a translated junction sequence, look for the fasta entry that overlaps with it, then return the entry
        and the coordinates. This will be used to extend said junction sequence to encompass the entire protein sequence.

        Instead of reading from the fasta file, it should just read fetch the Uniprot directly via API first. However
        if the Uniprot API returns an empty SwissProt result (somehow the Ensembl gene ID linking to a trembl rather than
        a Swissprot sequence), then we will fall back to a pre-loaded SwissProt only FASTA file.

        After extension, write the SeqRecord objects created into fasta file...

        :param seq:
        :param fasta:
        :return:
        """
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import requests as rq
        from io import StringIO
        import sys

        server = 'https://www.ebi.ac.uk'
        ext = '/proteins/api/proteins/Ensembl:' + self.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

        print(server + ext)

        ret = rq.get(server + ext, headers={"Accept": "text/x-fasta"})

        if not ret.ok:
            ret.raise_for_status()
            sys.exit()

        # If ret.text is empty, the fall back is to use a local fasta
        if len(ret.text) == 0:
            # Load local fasta (for now) based on species
            if species == 'mouse':
                fasta_handle = SeqIO.parse('data/fasta/20170918_Mm_Sp_16915.fasta', 'fasta', IUPAC.extended_protein)
            elif species == 'human':
                fasta_handle = SeqIO.parse('data/fasta/20170918_Hs_Sp_20205.fasta', 'fasta', IUPAC.extended_protein)
        else:
            fasta_handle = SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein)


        # The UniProt API retrieves a retrieval object, with a text field inside ret.text
        # Since Biopython SeqIO only works with file, use io.StringIO to turn the string into a file for parsing.

        for record in fasta_handle:

            # Find out where the first 10 amino acids meets the UniProt canonical sequences..
            merge_start1 = record.seq.find(self.slice1_aa[:10])
            merge_end1 = record.seq.find(self.slice1_aa[-10:])

            merge_start2 = record.seq.find(self.slice2_aa[:10])
            merge_end2 = record.seq.find(self.slice2_aa[-10:])

            # Only proceed to write file if we can bridge the first and last 10 amino acids of either
            # Slice 1 or slice 2 to the sequence.

            # Later on we should catch whether the first 10 aa is matched to multiple entries if using
            # the fallback protein fasta.
            if (merge_start1 != -1 and merge_end1 != -1) or (merge_start2 != -1 and merge_start2 != -1):

                # Write the UniProt canonical first
                canonical = record

                # Format the name of the slice 1 record
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + 10:]
                record1.id += ('|' + self.gene_id + '|' + self.junction_type + '1|'
                              + str(self.chr) + '_' + str(self.anc_ee) + '_' + str(self.alt1_ee)
                              + '_' + str(self.alt2_ee) + '|' + self.strand + str(self.phase))

                # If the slice is not the same as the UniProt canonical, then also write it.
                if record.seq.find(self.slice1_aa) == -1:
                    write_seqrecord_to_fasta(record1, output, suffix)

                # If not, then change name of canonical to reflect that it is also slice 1.
                else:
                    canonical.id = record1.id
                    write_seqrecord_to_fasta(canonical, output, suffix)

                # Format the name of the slice 2 record
                record2 = record[:merge_start2] + self.slice2_aa + record[merge_end2 + 10:]
                record2.id += ('|' + self.gene_id + '|' + self.junction_type + '2|'
                              + str(self.chr) + '_' + str(self.anc_ee) + '_' + str(self.alt2_ee) + '|'
                              + '_' + str(self.alt2_ee) + self.strand + str(self.phase))

                # If the slice is not the same as the UniProt canonical, then also write it.
                if record.seq.find(self.slice2_aa) == -1:
                    write_seqrecord_to_fasta(record2, output, suffix)

                # If not, then change name of canonical to reflect that it is also slice 2.
                else:
                    canonical.id = record2.id
                    write_seqrecord_to_fasta(canonical, output, suffix)


                print(record.id)
                print(record.seq[:merge_start1] + self.slice1_aa + record.seq[merge_end1 + 10:])
                print(record.seq[:merge_start2] + self.slice2_aa + record.seq[merge_end2 + 10:])

            # To do : if the slice is not matched to any of the FASTA entries,
            # We want to find out why. Also, write the slices to an orphan fasta
            else:
                print("==== SLICE IS NOT FOUND IN THE FASTA ==== ")
                print(self.slice1_aa)
                print(self.slice2_aa)
                print(record.seq)

                orphan_slice1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                                id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-1-' +
                                    self.name + '-' + str(self.phase) + self.strand),
                                name=self.gene_symbol,
                                description='Orphan Slice 1')
                orphan_slice2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
                                id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-2-' +
                                    self.name + '-' + str(self.phase) + self.strand),
                                name=self.gene_symbol,
                                description='Orphan Slice 2')

                write_seqrecord_to_fasta(orphan_slice1, output, suffix + '_orphan')
                write_seqrecord_to_fasta(orphan_slice2, output, suffix + '_orphan')


        return True

def write_seqrecord_to_fasta(seqrecord, output, suffix):
    """
    Write Biopython SeqRecord to the fasta file after checking whether the SeqRecord is already inside the file.

    :param seqrecord:   Biopython SeqRecord object
    :param sequence:    Splice Sequence object
    :param output:
    :param suffix:
    :return:
    """
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    import os.path


    os.makedirs('out', exist_ok=True)
    o = os.path.join('out', output + '_' + suffix + '.fasta')


    # If the file already exists, open it and amend that record.
    existing_records = []
    if os.path.exists(o):
        for existing_record in SeqIO.parse(o, 'fasta', IUPAC.extended_protein):
            existing_records.append(existing_record)

        # Test if the slice is already in the fasta, then do not write the new sequence into the fasta file.
        for existing_record in existing_records:
            if existing_record.seq == seqrecord.seq:
                print(seqrecord.seq)
                print("Already in fasta file - we will consider modifying the fasta entry name later on to reflect this")
                #return True

    output_handle = open(o, 'a')
    SeqIO.write(seqrecord, output_handle, 'fasta')
    output_handle.close()

    return True

#
# For doctest
#
if __name__ == '__main__':
    import doctest
    doctest.testmod()