# -*- coding: utf-8 -*-

"""jcast.main: Main function."""

import os
import datetime
import logging
from functools import partial

import tqdm
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import PowerTransformer
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

from jcast import params, fates
from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence
from jcast import __version__


def runjcast(args):
    """
    main look for jcast flow.

    :param args: parsed arguments
    :return:
    """

    # Get timestamp for out files
    now = datetime.datetime.now()

    write_dir = os.path.join(args.out, 'jcast_' + now.strftime('%Y%m%d%H%M%S'))
    os.makedirs(write_dir, exist_ok=True)

    # Main logger setup
    main_log = logging.getLogger('jcast')
    main_log.propagate = False
    main_log.setLevel(logging.INFO)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(write_dir, 'jcast_main.log'))
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

    # Model read count cutoff

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

    #
    # Model read count cutoff.
    # TODO: move this to a separate class
    #
    if args.model:

        main_log.info('The -m flag is set. The modeled read count will override -r --read values.')

        tot = rmats_results.rmats_mxe.append(rmats_results.rmats_se).append(rmats_results.rmats_ri).append(
            rmats_results.rmats_a5ss).append(rmats_results.rmats_a3ss).copy()
        junctions = [Junction(**tot.iloc[i].to_dict()) for i in range(len(tot))]

        # Array of junction sum read counts
        matrix_X = np.array([[j.sum_sjc + 1] for j in junctions])
        pt = PowerTransformer(method='box-cox')
        pt.fit(matrix_X)
        matrix_X_trans = pt.transform(matrix_X)

        gmm = GaussianMixture(n_components=2,
                              covariance_type='diag',
                              random_state=1,
                              ).fit(matrix_X_trans)

        # Sort the classification components so the high-read distribution appears second
        order = gmm.means_.argsort(axis=0)[:, 0]
        gmm.means_ = gmm.means_[order]
        gmm.covariances_ = gmm.covariances_[order]
        gmm.weights_ = gmm.weights_[order]
        gmm.precisions_ = gmm.precisions_[order]
        gmm.precisions_cholesky_ = gmm.precisions_cholesky_[order]

        weights = gmm.weights_
        means = gmm.means_
        covars = gmm.covariances_

        # Get the decision boundary (minimum count required to be predicted as second predicted class)
        # TODO: catch errors where model may fail and min_count remains 1
        min_count = 1
        # Loop only from 1 count to the mean of the second distribution
        for i in range(1, round(pt.inverse_transform([gmm.means_[1]])[0][0])):
            dist = gmm.predict(pt.transform(np.array([[i]])))
            if dist == [1]:
                min_count = i
                break

        # Plot out the model figure and decision boundary
        fig, ax = plt.subplots(constrained_layout=True)
        ax.hist(matrix_X_trans, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
        plt.style.use('seaborn-white')
        ax.set_xlabel('Box-Cox transformed total skipped junction counts')
        ax.set_ylabel('Frequency')
        ax.set_title('Skipped junction count distributions')
        tcks, values = plt.xticks()
        new_tcks = np.int_(np.sort(pt.inverse_transform(np.array([[t] for t in tcks])).ravel()))

        # Second axis for untransformed reads
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(tcks[1:-1])
        ax2.set_xlabel('Untransformed total splice junction read counts')
        ax2.set_xticklabels(new_tcks[1:-1])
        ax2.axvline(linewidth=4, ls='--', color='blue', x=pt.transform(np.array([[min_count]])))

        f_axis = matrix_X_trans.copy().ravel()
        f_axis.sort()
        ax.plot(f_axis, weights[0] * stats.norm.pdf(f_axis, means[0], np.sqrt(covars[0])).ravel(), c='red')
        ax.plot(f_axis, weights[1] * stats.norm.pdf(f_axis, means[1], np.sqrt(covars[1])).ravel(), c='blue')

        plt.savefig(fname=os.path.join(write_dir, 'model.png'))

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
                                        write_dir=write_dir,
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
        #         main_log.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(junctions[i].junction_type,
        #                                                                               junctions[i].name,
        #                                                                               junctions[i].gene_symbol,
        #                                                                               junctions[i].gene_id,
        #                                                                           ))
        #         main_log.info(f)

        #
        # Single threaded for-loop
        #
        for jx in tqdm.tqdm(junctions,
                            total=len(junctions),
                            desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
                            ):

            main_log.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(jx.junction_type,
                                                                                  jx.name,
                                                                                  jx.gene_symbol,
                                                                                  jx.gene_id,
                                                                                  ))
            main_log.info(translate_one_partial(jx))

    return True


def _translate_one(junction,
                   gtf,
                   genome,
                   args,
                   write_dir,
                   pred_bound,
                   ):
    """ get coordinate and translate one junction; arguments are passed through partial from main"""

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
    sequence = Sequence(junction)

    #
    # get nucleotide sequences of all slices using genome in memory
    # (anchor, alternative-1, alternative-2, downstream)
    # conjoin alternative exons to make slice 1 and 2,
    #
    sequence.make_slice_localgenome(genome.genome)

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
    if (not args.model and junction.sum_sjc < args.read) or (args.model and junction.sum_sjc < pred_bound):
        #
        # If the canonical flag is set, append the canonical
        # Sp to the gene_canonical output even if none of the transcript slices are stitchable
        # back to the canonical protein. This avoids not having any protein level representation
        # of a gene potentially in the proteome.
        #
        if args.canonical:
            sequence.write_canonical(outdir=write_dir)

        return fates.skipped_low

    #
    # discard junction if the corrected P value of this read count is < threshold
    # this removes junctions that are inconsistently found on both replicates.
    #
    q_lo, q_hi = args.qvalue
    if not q_lo < junction.fdr < q_hi:

        # Write canonical anyhow if the canonical flag is set.
        if args.canonical:
            sequence.write_canonical(outdir=write_dir)

        return fates.skipped_low

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
                outdir=write_dir,
                suffix='T1',
            )

            return fates.tier1.format(sequence.j.phase,
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
                outdir=write_dir,
                suffix='T2',
            )

            return fates.tier2.format(sequence.j.phase,
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

            sequence.write_slices(outdir=write_dir,
                                  suffix='T3',
                                  )

            return fates.tier3.format(sequence.j.phase,
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

        sequence.write_slices(outdir=write_dir,
                              suffix='T4',
                              )

        return fates.tier4.format(forced_slice)

    # force-translate through slice 1 if slice 1 hits PTC:
    elif len(sequence.slice2_aa) > 0 and len(sequence.slice1_aa) == 0:
        forced_slice = 1
        sequence.stitch_to_canonical_aa(slice_to_stitch=2)
        sequence.translate_forced(slice_to_translate=1)

        if len(sequence.slice1_aa) / len(sequence.slice2_aa) >= params.ptc_threshold:
            sequence.stitch_to_canonical_aa(slice_to_stitch=1,
                                            slice_has_ptc=True)

        sequence.write_slices(outdir=write_dir,
                              suffix='T4',
                              )

        return fates.tier4.format(forced_slice)

    #
    # if nothing works, write FAILURE fate
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
            sequence.write_canonical(outdir=write_dir)

        return fates.fail


def main():
    """ running main with parsed arguments from command line """

    import argparse
    import sys

    parser = argparse.ArgumentParser(description='jcast retrieves transcript splice junctions'
                                                 'and translates them into amino acid sequences')

    parser.add_argument('rmats_folder',
                        help='path to folder storing rMATS output',
                        )
    parser.add_argument('gtf_file',
                        help='path to Ensembl gtf file',
                        )
    parser.add_argument('genome',
                        help='path to genome file',
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

    parser.add_argument('-q', '--qvalue',
                        help='take junctions with rMATS fdr within this threshold [default: 0 1]',
                        metavar=('q_lo', 'q_hi'),
                        nargs=2,
                        default=[0, 1],
                        type=float)

    parser.set_defaults(func=runjcast)

    # print help message if no arguments are given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # parse all the arguments
    args = parser.parse_args()

    # run the function in the argument
    args.func(args)
