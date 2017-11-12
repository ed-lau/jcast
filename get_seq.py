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
        self.fate = 'Nothing done.'
        self.translated_phase = -1      # Phase that was actually used for translation.

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

    def translate(self, use_phase=True):
        """
        This is the translate function for tier 1/2 peptide. Calls make_pep with the terminal option
        which will terminate when it runs into a stop codon. Later on we should make a force_translate that does
        three frame translation and just either return all frames or return the longest.

        :param use_phase: T/F Whether to ue the stored phase or attempt to do three-frame translation
        :return: True
        """

        import helpers as h

        if self.phase in [0, 1, 2] and use_phase:
            self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, self.phase, terminate=True)
            self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, self.phase, terminate=True)
            print("Used Retrieved Phase")
            self.translated_phase = self.phase

        else:
            for i in range(3):
                self.slice1_aa = h.make_pep(self.slice1_nt, self.strand, i, terminate=True)
                self.slice2_aa = h.make_pep(self.slice2_nt, self.strand, i, terminate=True)
                if len(self.slice1_aa) > 0 and len(self.slice2_aa) > 0:
                    self.translated_phase = i
                    break
                else:
                    self.slice1_aa = ''
                    self.slice2_aa = ''

        print(self.slice1_aa)
        print(self.slice2_aa)

        return True


    def write_to_fasta(self, output, suffix):
        """

        WARNING: Deprecated, see the new "extend and write" function below, which will first attempt to harmonize
        the sequence to FASTA before writing the full-length sequence.

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
        and the coordinates. This will be used to extend said junction sequence to encompass  entire protein sequence.

        While this does not help resolve junction sequences or isoforms, I think it is important to match the database
        file as closely as possible to the actual spectra in a mass spectrometry experiment to strengthen the
        hypotheses being tested in the database search.

        Instead of reading from the FASTA file, it should just read fetch the Uniprot directly via API first. However
        if the Uniprot API returns an empty SwissProt result (somehow the Ensembl gene ID linking to a TremBL but not
        Swissprot sequence), then we will fall back to a pre-loaded SwissProt only FASTA file.

        After extension, write the SeqRecord objects created into fasta file...

        :param species: string  Specices (mouse or human) to determine which fasta to grab
        :param output:  string  Output directory
        :param: suffix: string  Additional suffix to add to the end of an output file
        :return:
        """
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        import requests as rq
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry
        from io import StringIO
        import sys
        import time
        import helpers as h


        server = 'https://www.ebi.ac.uk'
        ext = '/proteins/api/proteins/Ensembl:' + self.gene_id + '?offset=0&size=1&reviewed=true&isoform=0'

        print(server + ext)



        # retry 10 times

        retries = Retry(total=10,
                        backoff_factor=0.1,
                        status_forcelist=[500, 502, 503, 504])

        rqs = rq.Session()
        rqs.mount('https://', HTTPAdapter(max_retries=retries))
        ret = rqs.get(server + ext, headers={"Accept": "text/x-fasta"})

        if not ret.ok:
            print("Network still not okay after 10 retries. Quitting.")
            ret.raise_for_status()
            sys.exit()

        # If ret.text is empty, the fall back is to use a local fasta
        if len(ret.text) == 0:
            # Load local fasta (for now) based on species
            if species == 'mouse':
                fasta_handle = SeqIO.parse('data/fasta/20170918_Mm_Sp_16915.fasta', 'fasta',
                                           IUPAC.extended_protein)
            elif species == 'human':
                fasta_handle = SeqIO.parse('data/fasta/20170918_Hs_Sp_20205.fasta', 'fasta',
                                           IUPAC.extended_protein)
                # Don't save sequence if it is from the fall-back fasta file.


        # If ret.text is not empty, then get from online record.
        else:
            fasta_handle = SeqIO.parse(StringIO(ret.text), 'fasta', IUPAC.extended_protein)


        # The UniProt API retrieves a retrieval object, with a text field inside ret.text
        # Since Biopython SeqIO only works with file, use io.StringIO to turn the string into a file for parsing.

        for loop in fasta_handle:

            record = loop[:]

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
                canonical = record[:]

                # Format the name of the slice 1 record
                record1 = record[:merge_start1] + self.slice1_aa + record[merge_end1 + 10:]
                record1.id += ('|' + self.gene_id + '|' + self.junction_type + '1|' + self.name + '|'
                               + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                               + '|' + self.strand + str(self.phase)) + '|' + suffix

                # If the slice is different from the UniProt canonical, then also write it.
                if record.seq.find(self.slice1_aa) == -1:
                    h.write_seqrecord_to_fasta(record1, output, suffix)

                # If not, then change name of canonical to reflect that it is also slice 1.
                else:
                    canonical.id = record1.id
                    h.write_seqrecord_to_fasta(canonical, output, suffix)

                # Format the name of the slice 2 record
                record2 = record[:merge_start2] + self.slice2_aa + record[merge_end2 + 10:]
                record2.id += ('|' + self.gene_id + '|' + self.junction_type + '2|' + self.name + '|'
                               + str(self.chr) + '|' + str(self.anc_ee) + '|' + str(self.alt1_ee)
                               + '|' + self.strand + str(self.phase)) + '|' + suffix

                # If the slice is not the same as the UniProt canonical, then also write it.
                if record.seq.find(self.slice2_aa) == -1:
                    h.write_seqrecord_to_fasta(record2, output, suffix)

                # If not, then change name of canonical to reflect that it is also slice 2.
                else:
                    canonical.id = record2.id
                    h.write_seqrecord_to_fasta(canonical, output, suffix)


                print(record.id)
                print(record.seq[:merge_start1] + self.slice1_aa + record.seq[merge_end1 + 10:])
                print(record.seq[:merge_start2] + self.slice2_aa + record.seq[merge_end2 + 10:])



                # Once you found a match and wrote the sequence, quit.
                return True

        # If the slice is not matched to any of the FASTA entries,
        # write the slices to an orphan fasta.
        # NOTE: we are separating these out for now because we want to find out why they fall through.

        print("==== SLICE IS NOT FOUND IN THE FASTA ==== ")
        print(self.slice1_aa)
        print(self.slice2_aa)
        print(record.seq)

        orphan_slice1 = SeqRecord(Seq(self.slice1_aa, IUPAC.extended_protein),
                        id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-1-' +
                            self.name + '-' + str(self.phase) + '|' + self.strand),
                        name=self.gene_symbol,
                        description='Orphan Slice 1')
        orphan_slice2 = SeqRecord(Seq(self.slice2_aa, IUPAC.extended_protein),
                        id=(self.gene_symbol + '-' + self.gene_id + '-' + self.junction_type + '-2-' +
                            self.name + '-' + str(self.phase) + '|' + self.strand),
                        name=self.gene_symbol,
                        description='Orphan Slice 2')

        h.write_seqrecord_to_fasta(orphan_slice1, output, (suffix + '_orphan'))
        h.write_seqrecord_to_fasta(orphan_slice2, output, (suffix + '_orphan'))


        return True

    def write_fate(self, fate, output):
        """
        Write out the outcome of the attempt to translate each junction into a reort file

        :param fate:    int         Code for message to be writtebn
        :param output:  string      Output directory
        :return:
        """

        import os.path

        os.makedirs('out', exist_ok=True)
        o = os.path.join('out', output + '_' + 'fate' + '.txt')
        print(o)

        # Set the stored junction fate as the message
        self.fate = fate

        assert self.fate in [0, 1, 2, 3, 4], 'Junction fate code error.'

        if self.fate == 0:
            msg = ''

        elif self.fate == 1:
            msg = "SUCCESS 1. Retrieved phase: " + str(
                        self.phase) + " Used phase: " + str(self.translated_phase) + ". No Frameshift."

        elif self.fate == 2:
            msg = "SUCCESS 2. Retrieved phase: " + str(
                        self.phase) + " Used phase: " + str(self.translated_phase) + ". Frameshift."

        elif self.fate == 3:
            msg = "SUCCESS 3. The GTF frame appears to be wrong. Retrieved phase: " + str(
            self.phase) + " Used phase: " + str(self.translated_phase)

        elif self.fate == 4:
            msg = 'FAILURE. No translation was done. At least one PTC at each frame.'

        else:
            raise AssertionError

        f = open(o, 'a')
        f.write(self.junction_type + ' ' + self.name + ' ' + self.gene_symbol + ' ' + msg + '\n')
        f.close()

        return True



#
# For doctest
#
if __name__ == '__main__':
    import doctest
    doctest.testmod()