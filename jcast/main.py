# -*- coding: utf-8 -*-


"""jcast.main: Main function."""


import os
import datetime
import logging
import statistics
from functools import partial
import multiprocessing
import concurrent.futures

import tqdm

from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence

from jcast import __version__



def jcast(args):
    """
    main look for jcast flow.

    :param args: parsed arguments
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

    # number of threads for concurrency
    # cap it at cpu count -1 for now
    threads = min(args.num_threads, multiprocessing.cpu_count()-1)

    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #
    for rma in [rmats_results.rmats_mxe,
                rmats_results.rmats_se,
                rmats_results.rmats_ri,
                rmats_results.rmats_a5ss,
                rmats_results.rmats_a3ss,
                ]:

        junctions = [Junction(**rma.iloc[i].to_dict()) for i in range(len(rma))]

        translate_one_partial = partial(_translate_one,
                                        gtf=gtf,
                                        genome=genome,
                                        args=args,
                                        directory_to_write=directory_to_write,
                                        )

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            for i, f in enumerate(tqdm.tqdm(pool.map(
                    translate_one_partial,
                    junctions,
            ),
                    total=len(junctions),
                    desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
            )):
                main_log.info('>>>>>> Now doing junction {0} for gene {1}'.format(junctions[i].name,
                                                                                  junctions[i].gene_symbol))
                main_log.info(f)

    return True


def _translate_one(junction,
                   gtf,
                   genome,
                   args,
                   directory_to_write,
                    ):
    """ get coordinate and translate one junction; arguments are passed through partial from main"""
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


    #
    # Filter by minimal read counts
    #
    if junction.min_read_count < args.read:
        callback_ = ('SKIPPED. Sequence discarded due to low coverage. \n\n')
        return callback_

    # Discard this junction if the corrected P value of this read count is < 0.01
    # This is intended to remove junctions that aren't found on both replicates.
    # This might not be a good idea, however.
    if junction.fdr < args.pvalue:
        callback_ = 'SKIPPED. Sequence discarded due to difference across replicates. \n\n'
        return callback_

    #
    # Trim slice coordinates by translation starts and ends
    #
    junction.trim_cds(gtf)

    #
    # Get translated phase from GTF. Note this should be done after trimming to get the
    # right frame in case the exon in question is trimmed by the coding start
    #
    junction.get_translated_phase(gtf)

    #
    # Initiate a sequence object that copies most of the junction information
    #
    sequence = Sequence(junction)

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

    if args.sixframe:
        for strand, phase in [(strand, phase) for strand in ['+', '-'] for phase in [0, 1, 2]]:

            # Do six frame translation to get peptide
            sequence.translate_sixframe(strand, phase)

            # Check that the amino acid slices are at least as long as 10 amino acids)
            # Otherwise there is no point doing the merging
            if len(sequence.slice1_aa) >= 10 and len(sequence.slice2_aa) >= 10:
                # Extend with fasta, and then write if necessary.
                sequence.extend_and_write(output=directory_to_write,
                                          suffix='sixframe',
                                          merge_length=10)


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
            # TODO: instead of using SwissProt we should get the canonical exons from the GTF directly
            sequence.extend_and_write(
                output=directory_to_write,
                suffix='T1',
                merge_length=10,
            )

            callback_ = ('SUCCESS 1. Retrieved phase: {0} \n\n'
                         'Used phase: {1}. No frameshift.'.format(sequence.j.phase,
                                                                  sequence.translated_phase))

            return callback_

        #
        # Tier 2: both translated without stop codon, but with one frameshift
        #
        elif sequence.frameshift:
            sequence.extend_and_write(  # species=species,
                output=directory_to_write,
                suffix='T2',
                merge_length=10)

            callback_ = ('SUCCESS 2. Retrieved phase: {0} \n\n'
                          'Used phase: {1}. Frameshift.'.format(sequence.j.phase,
                                                                sequence.translated_phase))
            return callback_

    #
    # Tier 3 - retrieved phase is different from PTC-free frame.
    #
    else:
        sequence.translate(use_phase=False)

        #
        # After Tier 3 translation, check if both slices are good
        #

        if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
            sequence.extend_and_write(output=directory_to_write,
                                      suffix='T3',
                                      merge_length=10)

            callback_ = ('SUCCESS 3. GTF phase mismatch. Retrieved phase: {0} \n\n'
                          'Used phase: {1}'.format(sequence.j.phase,
                                                   sequence.translated_phase))
            return callback_

    #
    # If sequence is still not good, do Tier 4: One of the two slices hits stop codon.
    # (select one that is longest and translate if prior to the stop codon the short
    # slice is at least half as long as the long slice)
    # TODO: Determine likely translated frame, or keep frame based on PTC location from protein end
    #

    # Translate again after Tier 3 to reset to Tier 1/2 translation state
    sequence.translate(use_phase=True)

    # Force-translate slice 2 if slice 2 hits PTC:
    if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) == 0:

        sequence.translate_forced(slice_to_translate=2)

        if len(sequence.slice2_aa) / len(sequence.slice1_aa) >= 0.5:
            sequence.extend_and_write(output=directory_to_write,
                                      suffix='T4',
                                      merge_length=10)

            callback_ = 'PARTIAL 4.  Slice 2 hit a stop codon. Used longest phase.\n\n'
            return callback_

    # Force-translate slice 1 if slice 1 hits PTC:
    elif len(sequence.slice2_aa) > 0 and len(sequence.slice1_aa) == 0:

        sequence.translate_forced(slice_to_translate=1)

        if len(sequence.slice1_aa) / len(sequence.slice2_aa) >= 0.5:
            sequence.extend_and_write(output=directory_to_write,
                                      suffix='T4',
                                      merge_length=10)

            callback_ = 'PARTIAL 5.  Slice 1 hit a stop codon. Used longest phase.\n\n'
            return callback_

    #
    # If nothing works, write FAILURE fate
    #
    else:
        #
        # Salvage the canonical sequence in the long slice if it matches Sp exactly.
        # Note that this means if we identify a gene in RNA-seq, we will append the canonical
        # Sp to the gene_canonical output even if none of the transcript slices are stitchable
        # back to the canonical protein. This is to avoid not having any protein level representation
        # of a gene potentially in the proteome.
        #
        if args.canonical:
            sequence.extend_and_write(output=directory_to_write,
                                      merge_length=10,
                                      canonical_only=True)
        callback_ = 'FAILURE 6.  No alternative translation. \n\n'
        return callback_

    return True


def main():
    """ running main with parsed arguments from command line """

    import argparse
    import sys

    parser = argparse.ArgumentParser(description='jcast retrieves transcript splice junctions'
                                                 'and translates them into amino acid sequences')

    parser.add_argument('rmats_folder', help='path to folder storing rMATS output')
    parser.add_argument('gtf_file', help='path to Ensembl gtf file')
    parser.add_argument('genome', help='path to genome file')

    parser.add_argument('-n', '--num_threads', help='number of threads for concurrency [default: 6]',
                        default=6,
                        type=int)

    parser.add_argument('-o', '--out', help='name of the output files [default: psq_out]',
                        default='out')

    parser.add_argument('-r', '--read'
                        , help='minimum read counts to consider [default: 1]',
                        default=1,
                        type=int)

    parser.add_argument('-c', '--canonical', help='write out canonical protein sequence even if transcript'
                                   'slices are untranslatable [default: True]',
                        default=True,
                        type=bool)

    parser.add_argument('-p', '--pvalue'
                        , help='discard junctions with rMATS pvalue below this threshold [default: 0.01]',
                        default=0.01,
                        type=float)


    parser.add_argument('-s', '--sixframe', action='store_true',
                        help='also do six-frame translation instead with the junctions [default: False]',
                        )

    parser.set_defaults(func=jcast)

    # Print help message if no arguments are given

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
