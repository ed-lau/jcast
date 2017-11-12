#
# Splice Junction Mapper
# Edward Lau
#


from get_jxn import Junction, Annotation
from get_seq import Sequence
from get_rma import RmatsResults


def psqM(args):
    """
    Main loop for ProteoSeqM that controls logic flow.
    python main.py human data/encode_human_pancreas/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_pancreas_extended_retry


    :param args:
    :return:
    """

    import os

    #
    # Read arguments from ARg Parser
    #
    rmats_folder = args.rmats_folder
    gtf_loc = args.gtf_file
    out_file = args.out
    species = args.species

    #
    # Open the rMATS output file (MXE) here, rename the columns
    #

    assert os.path.exists(os.path.join(rmats_folder, 'MXE.MATS.JC.txt')), 'rMATS files not found, check directory.'

    rmats_results = RmatsResults(rmats_folder)

    #
    # Read the gtf file using the gtfpase package.
    # Then write as a pandas data frame.
    #

    gtf = Annotation(gtf_loc)
    gtf.read_gtf()


    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #

    for rmats_result in rmats_results.__dict__:
        # rmats_mxe, rmats_se, rmats_ri, rmats_a5ss, rmats_a3ss

        rma = rmats_results.__getattribute__(rmats_result)
        for i in range(len(rma)):

            # Set i to 1993 for rmats_mxe for PKM; PKM wasn't translated because the slice should have phase 0.
            # So phase detection actually failed.
            # Set i to 53 or 428 for rmats_mxe for an orphan slice

            junction = Junction(id=rma.id[i],
                                gene_id=rma.gene_id[i],
                                strand=rma.strand[i],
                                gene_symbol=rma.gene_symbol[i],
                                chr=rma.chr[i],
                                anc_es=rma.anc_es[i],
                                anc_ee=rma.anc_ee[i],
                                alt1_es=rma.alt1_es[i],
                                alt1_ee=rma.alt1_ee[i],
                                alt2_es=rma.alt2_es[i],
                                alt2_ee=rma.alt2_ee[i],
                                down_es=rma.down_es[i],
                                down_ee=rma.down_ee[i],
                                junction_type=rma.jxn_type[i],
                                species=species,)

            print('Analyzing' + rmats_result + str(i) + ' of ' + str(len(rma)))


            #
            # Subset the gtf file by the current gene_id
            #

            junction.get_translated_region(gtf)

            #
            # Get translation phase from GTF file.
            # If there is no phase found in the GTF, use phase -1 for now.
            # To do: look more closely into GTF file, or try translating from all frames
            #
            junction.get_translated_phase(gtf)

            #
            # Trim slice coordinates by translation starts and ends
            #
            junction.trim()

            #
            # Initiate a sequence object that copies most of the junction information
            #
            sequence = Sequence(junction)

            #
            # Get nucleotide sequences of all slices by REST API
            # (anchor, alternative-1, alternative-2, downstream)
            # Conjoin alternative exons to make slice 1 and 2,
            #
            sequence.make_slice()

            ## I think there should be a check here to see if there is a frameshift.
            ## See if slice 1 nucleotides are different in length from slice 2 nucleotide by
            ## multiples of 3
            if (len(sequence.slice1_nt) - len(sequence.slice2_nt)) % 3 != 0:
                print("Suspecting frameshifts between the two slices.")
                sequence.set_frameshift_to_true()

                # Note it looks like some frameshift skipped exon peptides could nevertheless come back in frame

            # We should only consider those without frameshift as tier 1.

            #
            # Translate into peptides
            #
            sequence.translate(use_phase=True)

            #
            # Write the tier 1 results into fasta
            #
            if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:

                # Tier 1: both translated without stop codon, no frameshift
                if not sequence.frameshift:

                    # Do a function like this to extend with fasta, and then write if necessary.
                    sequence.extend_and_write(species=species,
                                            output=out_file,
                                            suffix='T1')

                    fate_code = 1

                #
                # Tier 2: both translated without stop codon, but with frameshift
                #
                elif sequence.frameshift:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T2')

                    fate_code = 2

            else:
                #
                # Tier 3 - retrieved phase is wrong.
                #
                sequence.translate(use_phase=False)
                if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T3')

                    fate_code = 3


            # Tier 4: One hits stop codon (select one that is longest, use semi-supervised learning later).
                else:
                    fate_code = 4

                    pass

            sequence.write_fate(fate=fate_code, output=out_file)

    return True


#
# Code for running main with parsed arguments from command line
#

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='ProteoSeqM retrieves splice junction information'
                                                 'and translates into amino acid')




    # Create a "translate_rmats" subparser and house its specific arguments
    parser.add_argument('species', help='species (mouse or human)',
                        choices=['mouse', 'human'],
                        default='human')
    parser.add_argument('rmats_folder', help='path to folder storing rMATS output')
    parser.add_argument('gtf_file', help='path to ENSEMBL GTF file')

    parser.add_argument('-o', '--out', help='name of the output files',
                              default='psq_rmats')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose error messages.')

    parser.set_defaults(func=psqM)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)


