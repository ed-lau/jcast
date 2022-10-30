# -*- coding: utf-8 -*-

"""jcast.main: Main function."""

import os
import datetime
from functools import partial
import argparse
import tqdm

from jcast import params, fates, model
from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence
from jcast import __version__
from jcast.logger import get_logger


def run_jcast(args) -> None:
    """
    main look for jcast flow.

    :param args: parsed arguments
    :return: None
    """


    # ---- Main logger setup ----
    logger = get_logger('jcast', args.out)
    logger.info(args)
    logger.info(__version__)

    # ---- Read in rMATS files ----
    rmats_results = RmatsResults(logger=logger,
                                 rmats_dir=args.rmats_folder,
                                 )
    # TODO: only read in the splice types that are needed

    # ---- Read the gtf file using  gtfpase then write as a pandas data frame. ----
    gtf = ReadAnnotations(logger=logger,
                          path=args.gtf_file,
                          )

    # ---- Read genome file into memory ----
    genome = ReadGenome(logger=logger,
                        path=args.genome,
                        )

    #
    # Model read count cutoff.
    # TODO: move this to a separate class
    #
    if args.model:

        logger.info('The -m flag is set. The modeled read count will override -r --read values.')

        # Make a numpy array of all junction SJC sum counts
        rmats_results.get_junction_count_array()

        ln_sjc, best_mix_model, min_count = model.general_mixture_model(sum_sjc_array=rmats_results.sum_sjc_array)
        model.plot_general_mixture_model(
            ln_sjc,
            best_mix_model,
            min_count,
            write_dir=args.out,
            filename='model',
            )
        
        # Gaussian mixture model implementation
        """
        pt, gmm, min_count = model.gaussian_mixture(sum_sjc_array=rmats_results.sum_sjc_array)

        # Plot out the model
        model.plot_model(sum_sjc_array=rmats_results.sum_sjc_array,
                         pt=pt,
                         gmm=gmm,
                         min_count=min_count,
                         write_dir=args.out,
                         filename='model',
                         )
        """
    # If the m flag is not set, use the r argument value as min count
    else:
        min_count = args.read
    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #
    for splice_type, rma in zip(["MXE", "SE", "RI", "A5SS", "A3SS",],
                                [rmats_results.rmats_mxe, rmats_results.rmats_se, rmats_results.rmats_ri,
                                 rmats_results.rmats_a5ss, rmats_results.rmats_a3ss, ]):

        if splice_type not in args.splice_type:
            logger.info(f'Skipping {splice_type} because it was not specified in the -s argument.')
            continue
        # TODO: allow skipping of certain splice types

        junctions = [Junction(**rma.iloc[i].to_dict()) for i in range(len(rma))]

        translate_one_partial = partial(_translate_one,
                                        gtf=gtf,
                                        genome=genome,
                                        args=args,
                                        pred_bound=min_count,
                                        )

        #
        # Concurrent futures
        #
        # import concurrent.futures
        # with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()-1) as pool:
        #     for i, f in enumerate(tqdm.tqdm(pool.map(
        #             translate_one_partial,
        #             junctions,
        #     ),
        #             total=len(junctions),
        #             desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
        #     )):
        #         logger.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(junctions[i].junction_type,
        #                                                                               junctions[i].name,
        #                                                                               junctions[i].gene_symbol,
        #                                                                               junctions[i].gene_id,
        #                                                                           ))
        #         logger.info(f)

        #
        # Single threaded for-loop
        #
        for jx in tqdm.tqdm(junctions,
                            total=len(junctions),
                            desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
                            ):

            logger.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(jx.junction_type,
                                                                                  jx.name,
                                                                                  jx.gene_symbol,
                                                                                  jx.gene_id,
                                                                                  ))
            logger.info(translate_one_partial(jx))

    return True


def _translate_one(junction: Junction,
                   gtf: ReadAnnotations,
                   genome: ReadGenome,
                   args: argparse.Namespace,
                   pred_bound: int,
                   ):
    """
    Get coordinates and translate one junction object.
    :param junction: Junction object
    :param gtf: ReadAnnotations object
    :param genome: ReadGenome object
    :param args: parsed arguments
    :param pred_bound: read count cutoff
    :return:
    """

    #
    # trim slice coordinates by translation starts and ends
    #
    junction.trim_cds(gtf)

    #
    # get translated phase from GTF. Note this should be done after trimming to get the
    # right frame in case the exon in question is trimmed by the coding start
    #
    junction.get_translated_phase(gtf)

    #
    # initiate a sequence object that copies most of the junction information
    #
    sequence = Sequence(junction=junction)

    #
    # get nucleotide sequences of all slices using genome in memory
    # (anchor, alternative-1, alternative-2, downstream)
    # conjoin alternative exons to make slice 1 and 2,
    #
    sequence.make_slice_localgenome(genome.genome)

    # 2022-10-28: temporary fix to avoid hard masking
    if args.mask:
        for m in args.mask:
            if m in sequence.slice1_nt or m in sequence.slice2_nt:

                # Write canonical anyhow if the canonical flag is set.
                if args.canonical:
                    sequence.write_canonical(outdir=args.out)

                return fates.SKIPPED_MASKED

    #
    # translate to peptides
    #
    sequence.get_canonical_aa(gtf=gtf, genome_index=genome.genome)
    sequence.translate(use_phase=True)

    #
    # filter by junction read counts - discard junction if the min read count is below threshold
    #

    # If the -r argument is set directly and the -m flag is not, use the -r integer for count filtering
    # If the -m flag is set, use the modeled count for filtering
    if (not args.model and junction.sum_sjc <= args.read) or (args.model and junction.sum_sjc <= pred_bound):
        #
        # If the canonical flag is set, append the canonical
        # Sp to the gene_canonical output even if none of the transcript slices are stitchable
        # back to the canonical protein. This avoids not having any protein level representation
        # of a gene potentially in the proteome.
        #
        if args.canonical:
            sequence.write_canonical(outdir=args.out)

        return fates.SKIPPED_LOW_COUNT

    #
    # discard junction if the corrected P value of this read count is < threshold
    # this removes junctions that are inconsistently found on both replicates.
    #
    q_lo, q_hi = args.qvalue
    if not q_lo <= junction.fdr <= q_hi:

        # Write canonical anyhow if the canonical flag is set.
        if args.canonical:
            sequence.write_canonical(outdir=args.out)

        return fates.SKIPPED_P_VALUE

    #
    # write the Tier 1 and Tier 2 results into fasta file
    #
    if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:

        # Tier 1: both translated without stop codon, no frameshift
        if not sequence.frameshift:

            # Do a function like this to extend with fasta, and then write if necessary.
            # TODO: instead of using Uniprot we should get the canonical exons from the GTF directly

            for slice_ in [1, 2]:
                sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                slice_has_ptc=False)

            sequence.write_slices(
                outdir=args.out,
                suffix='T1',
            )

            return fates.TIER1.format(sequence.j.phase,
                                      sequence.translated_phase,
                                      )

        #
        # Tier 2: both translated without stop codon, but with one frameshift
        #
        elif sequence.frameshift:

            for slice_ in [1, 2]:
                sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                slice_has_ptc=False)

                # 2020-07-30 if slice runs into a frame shift,
                # allows the opportunity to stitch N-terminus only
                if [sequence.slice1_stitched, sequence.slice2_stitched][slice_-1] is None:
                    sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                    slice_has_ptc=True)

            sequence.write_slices(
                outdir=args.out,
                suffix='T2',
            )

            return fates.TIER2.format(sequence.j.phase,
                                      sequence.translated_phase,
                                      )

    #
    # Tier 3 - retrieved phase is different from PTC-free frame.
    #
    else:
        sequence.translate(use_phase=False)
        # after tier 3 translation, check if both slices are good
        if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
            for slice_ in [1, 2]:
                sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                slice_has_ptc=False)

            sequence.write_slices(outdir=args.out,
                                  suffix='T3',
                                  )

            return fates.TIER3.format(sequence.j.phase,
                                      sequence.translated_phase,
                                      )

    #
    # Tier 4: if sequence is still not good, do Tier 4: One of the two slices hits stop codon.
    # write out the slice if it is at least a certain proportion (params.ptc_threshold) as long as the long slice.
    #

    # translate again after tier 3 to reset to tier 1/2 translation state (using retrieved phase)
    sequence.translate(use_phase=True,
                       log=False,
                       )
    # TODO: we should avoid translating twice.

    # force-translate through slice 2 if slice 2 hits PTC:
    if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) == 0:
        forced_slice = 2
        sequence.stitch_to_canonical_aa(slice_to_stitch=1)
        sequence.translate_forced(slice_to_translate=forced_slice)

        if len(sequence.slice2_aa) / len(sequence.slice1_aa) >= params.ptc_threshold:
            sequence.stitch_to_canonical_aa(slice_to_stitch=2,
                                            slice_has_ptc=True)

        sequence.write_slices(outdir=args.out,
                              suffix='T4',
                              )

        return fates.TIER4.format(forced_slice)

    # force-translate through slice 1 if slice 1 hits PTC:
    elif len(sequence.slice2_aa) > 0 and len(sequence.slice1_aa) == 0:
        forced_slice = 1
        sequence.stitch_to_canonical_aa(slice_to_stitch=2)
        sequence.translate_forced(slice_to_translate=1)

        if len(sequence.slice1_aa) / len(sequence.slice2_aa) >= params.ptc_threshold:
            sequence.stitch_to_canonical_aa(slice_to_stitch=1,
                                            slice_has_ptc=True)

        sequence.write_slices(outdir=args.out,
                              suffix='T4',
                              )

        return fates.TIER4.format(forced_slice)

    #
    # if nothing works, write FAIL fate
    #
    else:
        #
        # salvage the canonical sequence in the long slice if it matches Sp exactly.
        # note that this means if we identify a gene in RNA-seq, we will append the canonical
        # Sp to the gene_canonical output even if none of the transcript slices are stitchable
        # back to the canonical protein. This is to avoid not having any protein level representation
        # of a gene potentially in the proteome.
        #
        if args.canonical:
            sequence.write_canonical(outdir=args.out)

        return fates.FAIL


# Check the q value ranges are between 0 and 1 and the second value is larger than the first.
class CheckQValueRange(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not 0 <= values[0] <= 1:
            parser.error("qvalue lower bound must be between 0 and 1")
        if not 0 <= values[1] <= 1:
            parser.error("qvalue upper bound must be between 0 and 1")
        if values[0] >= values[1]:
            parser.error("qvalue lower bound must be smaller than upper bound")
        setattr(namespace, self.dest, values)


def main():
    """ running main with parsed arguments from command line """


    import sys

    parser = argparse.ArgumentParser(description='jcast retrieves transcript splice junctions'
                                                 'and translates them into amino acid sequences')

    parser.add_argument('rmats_folder',
                        help='path to folder storing rMATS output',
                        )
    parser.add_argument('gtf_file',
                        help='path to Ensembl gtf file',
                        type=argparse.FileType('r'),
                        )
    parser.add_argument('genome',
                        help='path to genome file',
                        type=argparse.FileType('r'),
                        )

    # parser.add_argument('-n', '--num_threads', help='number of threads for concurrency [default: 6]',
    #                     default=6,
    #                     type=int)

    parser.add_argument('-o', '--out',
                        help='name of the output files [default: psq_out]',
                        default='out')

    parser.add_argument('-r', '--read',
                        help='the lowest skipped junction read count for a junction to be translated [default: 1]',
                        default=1,
                        type=int,
                        )

    parser.add_argument('-m', '--model',
                        help='models junction read count cutoff using a Gaussian mixture model [default: False]',
                        action='store_true',
                        default=False,
                        #type=bool,
                        )

    parser.add_argument('-c', '--canonical', help='write out canonical protein sequence even if transcript'
                                   'slices are untranslatable [default: False]',
                        default=False,
                        action='store_true',
                        # type=bool,
                        )

    parser.add_argument('--mask',
                        help='masked nucleotides; slices will not be translated if '
                             'containing one of these letters [default: NnXx]',
                        default='NnXx',
                        )

    parser.add_argument('-s', '--splice_type',
                       help='splice type to be translated [default: MXE SE RI A5SS A3SS]',
                       default=set(['MXE', 'SE', 'RI', 'A5SS', 'A3SS']),
                       choices=['MXE', 'SE', 'RI', 'A5SS', 'A3SS'],
                       nargs='+',
                       )

    parser.add_argument('-q', '--qvalue',
                        help='take junctions with rMATS fdr within this threshold [default: 0 1]',
                        metavar=('q_lo', 'q_hi'),
                        nargs=2,
                        default=[0, 1],
                        type=float,
                        action=CheckQValueRange,
                        )

    # DEVELOPMENT ARGUMENT, Should be deleted when a default distribution - Gamma or LogNorm - is determined.
    parser.add_argument("--g_or_ln",
                        help="Switch on distribution to use for low end of histogram, 0 for Gamma, anything else for LogNorm",
                        default=0,
                        type=int)

    parser.set_defaults(func=run_jcast)

    # print help message if no arguments are given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # parse all the arguments
    args = parser.parse_args()

    # Check if the rMATS directory contains the expected files
    for splice_type in args.splice_type:
        assert os.path.exists(os.path.join(args.rmats_folder, f'{splice_type}.MATS.JC.txt')), \
            'rMATS files not found, check directory.'

    # run the function in the argument
    args.func(args)
