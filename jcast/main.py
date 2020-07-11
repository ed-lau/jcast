# -*- coding: utf-8 -*-


"""jcast.main: Main function."""


import os

import datetime
import logging

import tqdm
import numpy as np

from jcast.get_jxn import Junction
from jcast.get_gtf import ReadAnnotations
from jcast.get_seq import Sequence
from jcast.get_rma import RmatsResults
from jcast.read_fa import ReadGenome

from jcast import __version__



def psqM(args):
    """
    Main loop for JCast that controls logic flow.

    Usage:
    python -m jcast /path/to/rmats_folder/ path/to/gtf path/to/genome.fa  -o path/to/out -r 5


    >>> gtf = ReadAnnotations('../data/gtf/Homo_sapiens.GRCh38.89.gtf')
    >>> gtf.read_gtf()
    True
    >>> rmats_results = RmatsResults(rmats_dir='../data/encode_human_heart/')
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

    # Get timestamp for out files
    now = datetime.datetime.now()

    directory_to_write = os.path.join(args.out, 'jcast_' + now.strftime('%Y%m%d%H%M%S'))
    os.makedirs(directory_to_write, exist_ok=True)

    # Main logger setup
    main_log = logging.getLogger('jcast')
    main_log.propagate = False
    main_log.setLevel(logging.INFO)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(directory_to_write, 'jcast_main.log'))
    fh.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)

    # add the handlers to the logger
    main_log.addHandler(fh)

    #
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    ch.setFormatter(formatter)
    main_log.addHandler(ch)

    main_log.info(args)
    main_log.info(__version__)

    #
    # Read arguments from ArgParser
    #
    # species = args.species

    #
    # Open the rMATS output file (MXE) here, rename the columns
    #

    assert os.path.exists(os.path.join(args.rmats_folder, 'MXE.MATS.JC.txt')), 'rMATS files not found, check directory.'
    rmats_results = RmatsResults(rmats_dir=args.rmats_folder)

    #
    # Read the gtf file using the gtfpase package.
    # Then write as a pandas data frame.
    #
    gtf = ReadAnnotations(args.gtf_file)
    gtf.read_gtf()

    #
    # Read genome file into memory
    #
    genome = ReadGenome(args.genome)

    # For the decommissioned resume hack
    #
    # list_of_done_junctions = []
    # try:
    #     rf = os.path.join('out', out_file + '_' + 'fate' + '.txt')
    #     import csv
    #     existing_fate_text = csv.reader(open(rf, 'r'), delimiter='\t')
    #
    #     for row in existing_fate_text:
    #         list_of_done_junctions.append(row[0] + '_' + row[1])
    #
    # except FileNotFoundError:
    #    print('No existing partial fate file found.')

    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #

    for rma in [rmats_results.rmats_mxe,
                rmats_results.rmats_se,
                rmats_results.rmats_ri,
                rmats_results.rmats_a5ss,
                rmats_results.rmats_a3ss,
                ]:

        for i in tqdm.tqdm(range(len(rma)),
                           desc='Processing {0} Junctions'.format(rma.jxn_type[0])):

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
                                #species=species,

                                )

            #
            # Code for checking and executing whether to resume run
            #
            # if (str(rma.jxn_type[i]) + '_' + str(rma.id[i]) in list_of_done_junctions) and args.resume:
            #     print('resume 1: skipping existing junction' + rmats_result + str(i) + '.')
            #     continue
            # else:
            #     print('Analyzing ' + rmats_result + str(i) + ' of ' + str(len(rma)))

            main_log.info('>>>>>> Now doing junction {0} for gene {1}'.format(junction.name,
                                                                            junction.gene_symbol))

            #
            # Code for filtering by rMATS results
            #

            # Discard this junction if the read count is below threshold in both splice junctions.
            # Maybe should change this to discard everything with read counts below threshold on EITHER junction
            # This is intended to remove junctions that are very low in abundance.

            # If rMATS was run with one technical replicate, the count field is an int, otherwise it is a list
            # The following should take care of both single integer and list of integers.

            # First we mandate that there is a read spanning the junction (taking the rMTAS JC rather than JCEC
            # files.

            # We are taking the skipped junction count (SJC) as filtering criterion for now
            # because the majority of translatable events are probably SE (skipped exon)
            # Essentially this filters out alternative junctions that are very rarely skipped
            # (high inclusion level of the exons) that are not likely to be translatable.
            try:
                mean_count_sample1 = int(np.mean([int(x) for x in (str(rma.sjc_s1[i]).split(sep=','))]))
                mean_count_sample2 = int(np.mean([int(x) for x in (str(rma.sjc_s2[i]).split(sep=','))]))

            except ValueError:
                mean_count_sample1 = 0
                mean_count_sample2 = 0

                main_log.info('SKIPPED. Discarding sequence since unable to find read counts. \n\n')

            junction.set_min_read_count(min([mean_count_sample1, mean_count_sample2]))

            # Filter by minimal read counts
            if mean_count_sample1 < args.read or mean_count_sample2 < args.read:
                main_log.info('SKIPPED. Sequence discarded due to low coverage. \n\n')
                continue

            # Discard this junction if the corrected P value of this read count is < 0.01
            # This is intended to remove junctions that aren't found on both replicates.
            # This might not be a good idea, however.
            if rma.fdr[i] < args.pvalue:
                # Initiate a dummy sequence just to write to the fate file
                # We probably want a better solution for this (decorator?)
                main_log.info('SKIPPED. Sequence discarded due to difference across replicates. \n\n')

            #
            # Subset the gtf file by the current gene_id
            #
            junction.get_translated_region(gtf)
            main_log.info('Anchor exon start: ' + str(junction.anc_es) + ' Anchor exon end: ' + str(junction.anc_ee))

            #
            # Get translation phase from GTF file.
            # If there is no phase found in the GTF, use phase -1 for now.
            # To do: look more closely into GTF file, or try translating from all frames
            #
            junction.get_translated_phase(gtf)
            main_log.info('Transcription start: ' + str(junction.tx0) + ' Transcription end:' + str(junction.tx1))
            main_log.info('Retrieved phase: ' + str(junction.phase))
            #
            # Trim slice coordinates by translation starts and ends
            #
            junction.trim()

            #
            # Initiate a sequence object that copies most of the junction information
            #
            sequence = Sequence(junction,
                                directory_to_write=directory_to_write,
                                )

            #
            # Get nucleotide sequences of all slices using genome in memory
            # (anchor, alternative-1, alternative-2, downstream)
            # Conjoin alternative exons to make slice 1 and 2,
            #
            sequence.make_slice_localgenome(genome.genome)

            #
            # The next section is the six-frame translational by-pass. If the --sixframe flag is on,
            # Then do six-frame with all the qualifying junctions instead
            #
            # if args.sixframe:
            #     if args.verbose:
            #         print('verbose 1: using six-frame translation instead of reading frames from gtf.')
            #
            #     for sixframe_strand in ['+', '-']:
            #         for sixframe_phase in [0, 1, 2]:
            #
            #             # Do six frame translation to get peptide
            #             sequence.translate_sixframe(sixframe_strand, sixframe_phase)
            #
            #             # Check that the amino acid slices are at least as long as 10 amino acids)
            #             # Otherwise there is no point doing the merging
            #             if len(sequence.slice1_aa) >= 10 and len(sequence.slice2_aa) >= 10:
            #                 # Extend with fasta, and then write if necessary.
            #                 sequence.extend_and_write(species=species,
            #                                           output=out_file,
            #                                           suffix='6F',
            #                                           merge_length=10)
            #
            #     fate_code = 7
            #     sequence.write_fate(fate=fate_code, output=out_file)
            #     continue

            # I think there should be a check here to see if there is a frameshift.
            # See if slice 1 nucleotides are different in length from slice 2 nucleotide by
            # multiples of 3, which probably denotes frame shift (unless there are loose amino acids near the end)?
            if (len(sequence.slice1_nt) - len(sequence.slice2_nt)) % 3 != 0:
                sequence.set_frameshift_to_true()
                main_log.info("Frame-shift between the two slices (length difference not multiples of 3).")

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
                    sequence.extend_and_write(#species=species,
                                              output=directory_to_write,
                                              suffix='T1',
                                              merge_length=10)

                    main_log.info("SUCCESS 1. Retrieved phase: {0} \n\n"
                                  "Used phase: {1}. No frameshift.".format(sequence.phase,
                                                                           sequence.translated_phase))

                    continue

                #
                # Tier 2: both translated without stop codon, but with frameshift
                #
                elif sequence.frameshift:
                    sequence.extend_and_write(#species=species,
                                              output=directory_to_write,
                                              suffix='T2',
                                              merge_length=10)

                    main_log.info("SUCCESS 2. Retrieved phase: {0} \n\n"
                                  "Used phase: {1}. Frameshift.".format(sequence.phase,
                                                                        sequence.translated_phase))
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
                    sequence.extend_and_write(#species=species,
                                              output=directory_to_write,
                                              suffix='T3',
                                              merge_length=10)

                    main_log.info("SUCCESS 3. GTF phase mismatch. Retrieved phase: {0} \n\n"
                                  "Used phase: {1}".format(sequence.phase,
                                                           sequence.translated_phase))
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
                    sequence.extend_and_write(#species=species,
                                              output=directory_to_write,
                                              suffix='T4',
                                              merge_length=10)

                    main_log.info('PARTIAL 4.  Slice 2 hit a stop codon. Used longest phase.\n\n')
                    continue

            # Do this if slice 1 hits PTC:
            elif len(sequence.slice2_aa) > 0 and len(sequence.slice1_aa) == 0:


                sequence.translate_forced(slice_to_translate=1)

                if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
                    sequence.extend_and_write(#species=species,
                                              output=directory_to_write,
                                              suffix='T4',
                                              merge_length=10)

                    main_log.info('PARTIAL 5.  Slice 1 hit a stop codon. Used longest phase.\n\n')
                    continue

            #
            # If nothing works, write FAILURE fate
            #

            elif len(sequence.slice1_aa) == 0 and len(sequence.slice2_aa) == 0:
                main_log.info('FAILURE 6.  No translation was done. At least one PTC at each frame.\n\n')

    return True


#
# Code for running main with parsed arguments from command line
#

def main():

    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Jcast retrieves splice junction information'
                                                 'and translates into amino acid')

    parser.add_argument('rmats_folder', help='path to folder storing rMATS output')
    parser.add_argument('gtf_file', help='path to ENSEMBL GTF file')
    parser.add_argument('genome', help='path to Genome file')

    parser.add_argument('-o', '--out', help='name of the output files [default: psq_out]',
                        default='out')

    parser.add_argument('-r', '--read'
                        , help='minimum read counts to consider [default: 1]',
                        default=1,
                        type=int)

    parser.add_argument('-p', '--pvalue'
                        , help='discard junctions with rMATS pvalue below this threshold [default: 0.01]',
                        default=0.01,
                        type=float)


    # parser.add_argument('-s', '--sixframe', action='store_true',
    #                     help='do six-frame translation instead with the junctions')

    parser.set_defaults(func=psqM)

    # Print help message if no arguments are given

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
