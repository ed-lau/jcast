"""
Splice Junction Mapper v.1.0.0
Molecular Proteomics Laboratory
http://maggielab.org

Edward Lau 2017-2018
lau1@stanford.edu

"""

# Import classes
from get_jxn import Junction, Annotation
from get_seq import Sequence
from get_rma import RmatsResults
import numpy as np


def psqM(args):
    """
    Main loop for ProteoSeqM that controls logic flow.

    Usage:
    python main.py human data/encode_human_pancreas/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_pancreas_extended_retry
    python main.py human data/encode_human_liver/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_liver_extended_v3 -r
    python main.py human data/encode_human_heart/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_heart_v5 -r -f
    python main.py human data/encode_human_adrenalgland/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_adrenalgland_extended_v4 -r -f
    python main.py human data/encode_human_transversecolon/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_transversecolon_extended_v4 -r -f
    python main.py human data/encode_human_testis/ data/gtf/Homo_sapiens.GRCh38.89.gtf -o psqnew_encode_human_testis_extended_v4 -r -f

    >>> gtf = Annotation('data/gtf/Homo_sapiens.GRCh38.89.gtf')
    >>> gtf.read_gtf()
    True
    >>> rmats_results = RmatsResults('data/encode_human_heart/')
    >>> rma = rmats_results.__getattribute__('rmats_mxe')
    >>> i = 1993
    >>> junction = Junction(id=rma.id[i],\
                                gene_id=rma.gene_id[i],\
                                strand=rma.strand[i],\
                                gene_symbol=rma.gene_symbol[i],\
                                chr=rma.chr[i],\
                                anc_es=rma.anc_es[i],\
                                anc_ee=rma.anc_ee[i],\
                                alt1_es=rma.alt1_es[i],\
                                alt1_ee=rma.alt1_ee[i],\
                                alt2_es=rma.alt2_es[i],\
                                alt2_ee=rma.alt2_ee[i],\
                                down_es=rma.down_es[i],\
                                down_ee=rma.down_ee[i],\
                                junction_type=rma.jxn_type[i],\
                                species='human',)
    >>> junction.get_translated_region(gtf)
    Anchor exon start: 72200474 Anchor exon end: 72200655
    Transcription start: 72199652 Transcription end:72219097
    True
    >>> junction.trim()
    True


    :param args:
    :return:
    """

    import os

    #
    # Read arguments from ArgParser
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
    # Attempt to read the fate file and then, if resume flag is on, skip the rmats_result that are already done
    #
    os.makedirs('out', exist_ok=True)
    list_of_done_junctions = []

    try:
        rf = os.path.join('out', out_file + '_' + 'fate' + '.txt')
        import csv
        existing_fate_text = csv.reader(open(rf, 'r'), delimiter='\t')

        for row in existing_fate_text:
            list_of_done_junctions.append(row[0] + '_' + row[1])

    except FileNotFoundError:
        print('No existing partial fate file found.')

    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #

    for rmats_result in rmats_results.__dict__:
        # rmats_mxe, rmats_se, rmats_ri, rmats_a5ss, rmats_a3ss

        rma = rmats_results.__getattribute__(rmats_result)
        for i in range(len(rma)):


            # To access with pandas, rma.ix[:,'sjc_s1'], etc., rma.ix[i]

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

            #
            # Code for checking and executing whether to resume run
            #
            if (str(rma.jxn_type[i]) + '_' + str(rma.id[i]) in list_of_done_junctions) and args.resume:
                print('resume 1: skipping existing junction' + rmats_result + str(i) + '.')
                continue
            else:
                print('Analyzing ' + rmats_result + str(i) + ' of ' + str(len(rma)))

            #
            # Code for filtering by rMATS results
            #
            if args.filter:

                # Discard this junction if the read count is below 10 in both splice junctions.
                # Maybe should change this to discard everything with read counts below 10 on EITHER junction
                # This is intended to remove junctions that are very low in abundance.

                # If rMATS was run with one technical replicate, the count field is an int, otherwise it is a list
                # The following should take care of both single integer and list of integers.
                try:
                    mean_count_sample1 = int(np.mean([int(x) for x in (str(rma.sjc_s1[i]).split(sep=','))]))
                    mean_count_sample2 = int(np.mean([int(x) for x in (str(rma.sjc_s2[i]).split(sep=','))]))

                except ValueError:
                    mean_count_sample1 = 4
                    mean_count_sample2 = 4
                    if args.verbose:
                        print('verbose 1: discarding sequence since unable to find read counts.')

                if mean_count_sample1 < 10 or mean_count_sample2 < 10:
                    fate_code = -2
                    sequence = Sequence(junction)
                    # Initiate a dummy sequence just to write to the fate file
                    # We probably want a better solution for this (decorator?)

                    sequence.write_fate(fate=fate_code, output=out_file)
                    if args.verbose:
                        print('verbose 1: sequence discarded due to low coverage.')
                    continue

                # Discard this junction if the corrected P value of this read count is < 0.01
                # This is intended to remove junctions that aren't found on both replicates.
                # This might not be a good idea, however.
                if rma.fdr[i] < 0.01:
                    fate_code = -1

                    # Initiate a dummy sequence just to write to the fate file
                    # We probably want a better solution for this (decorator?)
                    sequence = Sequence(junction)
                    sequence.write_fate(fate=fate_code, output=out_file)
                    if args.verbose:
                        print('verbose 1: sequence discarded due to inconsistency across replicates.')
                    continue

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

            #
            # The next section is the six-frame translational by-pass. If the --sixframe flag is on,
            # Then do six-frame with all the qualifying junctions instead
            #

            if args.sixframe:
                if args.verbose:
                    print('verbose 1: using six-frame translation instead of reading frames from gtf.')

                for sixframe_strand in ['+', '-']:
                    for sixframe_phase in [0, 1, 2]:

                        # Do six frame translation to get peptide
                        sequence.translate_sixframe(sixframe_strand, sixframe_phase)

                        # Check that the amino acid slices are at least as long as 10 amino acids)
                        # Otherwise there is no point doing the merging
                        if len(sequence.slice1_aa) >= 10 and len(sequence.slice2_aa) >= 10:

                            # Extend with fasta, and then write if necessary.
                            sequence.extend_and_write(species=species,
                                                      output=out_file,
                                                      suffix='6F',
                                                      merge_length=10)

                fate_code = 7
                sequence.write_fate(fate=fate_code, output=out_file)
                continue

            # I think there should be a check here to see if there is a frameshift.
            # See if slice 1 nucleotides are different in length from slice 2 nucleotide by
            # multiples of 3, which probably denotes frame shift (unless there are loose amino acids near the end)?
            if (len(sequence.slice1_nt) - len(sequence.slice2_nt)) % 3 != 0:
                sequence.set_frameshift_to_true()
                if args.verbose:
                    print("verbose 1: suspecting frameshifts between the two slices.")

                # Note it looks like some frameshift skipped exon peptides could nevertheless come back in frame
                # We should only consider those without frameshift as tier 1.

            #
            # Translate into peptides
            #
            sequence.translate(use_phase=True)

            #
            # Write the Tier 1 and Tier 2 results into fasta file
            #
            if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:

                # Tier 1: both translated without stop codon, no frameshift
                if not sequence.frameshift:

                    # Do a function like this to extend with fasta, and then write if necessary.
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T1',
                                              merge_length=10)

                    fate_code = 1
                    sequence.write_fate(fate=fate_code, output=out_file)
                    continue

                #
                # Tier 2: both translated without stop codon, but with frameshift
                #
                elif sequence.frameshift:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T2',
                                              merge_length=10)

                    fate_code = 2
                    sequence.write_fate(fate=fate_code, output=out_file)
                    continue

            #
            # Tier 3 - retrieved phase is wrong.
            #

            else:
                sequence.translate(use_phase=False)

            #
            # After Tier 3 translation, check if both slices are good
            #
                if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T3',
                                              merge_length=10)

                    fate_code = 3
                    sequence.write_fate(fate=fate_code, output=out_file)
                    continue

            #
            # If sequence is still not good, do Tier 4: One of the two slices hits stop codon.
            # (select one that is longest, use semi-supervised learning later).
            #

            # Translate again after Tier 3 to reset to Tier 1/2 translation state
            sequence.translate(use_phase=True)

            # Do this if slice 2 hits PTC:
            if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) == 0:
                sequence.translate_forced(slice_to_translate=2)

                if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T4',
                                              merge_length=10)

                    fate_code = 4
                    sequence.write_fate(fate=fate_code, output=out_file)
                    continue

            # Do this if slice 1 hits PTC:
            elif len(sequence.slice2_aa) > 0 and len(sequence.slice1_aa) == 0:

                sequence.translate_forced(slice_to_translate=1)

                if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
                    sequence.extend_and_write(species=species,
                                              output=out_file,
                                              suffix='T4',
                                              merge_length=10)

                    fate_code = 5
                    sequence.write_fate(fate=fate_code, output=out_file)
                    continue

            #
            # If nothing works, write FAILURE fate
            #
            if len(sequence.slice1_aa) == 0 and len(sequence.slice2_aa) == 0:
                fate_code = 6
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
    parser.add_argument('-r', '--resume', action='store_true', help='attempt to resume interrupted runs')
    parser.add_argument('-f', '--filter', action='store_true', help='filter junctions based on read counts')
    parser.add_argument('-s', '--sixframe', action='store_true',
                        help='do six-frame translation instead with the junctions')

    parser.add_argument('-v', '--verbose', action='store_true', help='verbose error messages')


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
