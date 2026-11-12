# -*- coding: utf-8 -*-

"""jcast.main: Main function."""

import os
import collections
import glob
import argparse
import logging
import multiprocessing
import traceback
import tqdm

from jcast import params, fates, model, cache
from jcast import helpers as h
from jcast.junctions import Junction, RmatsResults
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence
from jcast.uniprot import UniProtFasta, MERGED_FILE, ISOFORMS_FILE, CANONICALS_FILE
from jcast import __version__
from jcast.logger import get_logger, close_logger

# The files a run writes into its output directory. A directory holding any of them is used again only with
# --overwrite, which removes them first (reports/2026-11-12-thresholds-runs-canonical-audit.md §3.2)
OUTPUT_FILES = ['logfile.log', 'model.pdf', MERGED_FILE, ISOFORMS_FILE, CANONICALS_FILE]
OUTPUT_PATTERN = 'psq_*.fasta'


def existing_outputs(out_dir) -> list:
    """ JCAST output files already in a directory """
    if not os.path.isdir(out_dir):
        return []
    found = [f for f in OUTPUT_FILES if os.path.isfile(os.path.join(out_dir, f))]
    found += [os.path.basename(f) for f in glob.glob(os.path.join(glob.escape(out_dir), OUTPUT_PATTERN))]
    return sorted(found)


def run_jcast(args) -> None:
    """
    Main function for jcast

    The output directory is created if needed. If it holds output of an earlier run, the run is refused unless
    args.overwrite is set, which removes that output first: appending to it mixed the two runs (REVIEW.md MED-12).

    :param args: parsed arguments
    :return: None
    """

    os.makedirs(args.out, exist_ok=True)
    found = existing_outputs(args.out)
    if found and not getattr(args, 'overwrite', False):
        raise FileExistsError(f'{args.out} already holds JCAST output ({", ".join(found)}); '
                              f'use --overwrite to replace it, or another output directory.')
    for f in found:
        os.remove(os.path.join(args.out, f))

    cache.set_directory(getattr(args, 'cache_dir', None))

    # ---- Main logger setup ----
    logger = get_logger('jcast', args.out)
    try:
        if found:
            logger.info(f'--overwrite: removed the output of an earlier run: {", ".join(found)}.')
        return _run(args, logger)
    finally:
        close_logger(logger)


def _run(args, logger) -> bool:
    """ The run, logging to logger """
    logger.info(args)
    logger.info(__version__)

    if getattr(args, 'mask', None) is not None:
        logger.warning('--mask is deprecated and ignored: slices containing any base other than A, C, G or T are '
                       'always skipped. It will be removed in the next release.')

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

    # ---- Read the UniProt FASTA if supplied; otherwise canonicals come from the UniProt API ----
    uniprot = UniProtFasta(logger=logger, path=args.uniprot) if getattr(args, 'uniprot', None) else None

    #
    # Model read count cutoff.
    # TODO: move this to a separate class
    #
    if args.model:

        logger.info('The -m flag is set. The modeled read count will override -r --read values.')

        # Make a numpy array of all junction SJC sum counts
        rmats_results.get_junction_count_array()

        ln_sjc, mix_params, min_count = model.general_mixture_model(min_jc_array=rmats_results.min_jc_array)

        if min_count is None:
            logger.warning(f'Mixture model fit was degenerate; falling back to -r {args.read}.')
            min_count = args.read
        else:
            logger.info('Mixture model parameters (gamma shape, gamma rate, normal mean, normal sd, weights): '
                        f'{[round(float(p), 3) for p in mix_params]}; modeled read count cutoff {min_count}.')
            model.plot_general_mixture_model(
                ln_sjc,
                mix_params,
                min_count,
                write_dir=args.out,
                filename='model',
                )

        # Gaussian mixture model implementation
        """
        pt, gmm, min_count = model.gaussian_mixture(min_jc_array=rmats_results.min_jc_array)

        # Plot out the model
        model.plot_model(min_jc_array=rmats_results.min_jc_array,
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

    logger.info(f'Translating junctions with a read count of at least {min_count} (-{"m" if args.model else "r"}): '
                f'the smaller of the summed inclusion and skipping junction counts. Others are skipped.')

    # ---- Worker processes (-p). Started before any progress bar, so no thread exists when they are forked ----
    pool = None
    processes = _number_of_processes(args, uniprot, logger)
    if processes > 1:
        _WORKER.update(gtf=gtf, genome=genome, args=args, pred_bound=min_count, uniprot=uniprot, reported=set())
        pool = multiprocessing.get_context('fork').Pool(processes, initializer=_init_worker)
        logger.info(f'Translating junctions in {processes} worker processes.')

    try:
        failed = _translate_all(args, logger, rmats_results, gtf, genome, min_count, uniprot, pool)
    except BaseException:
        if pool is not None:
            pool.terminate()
        raise
    if pool is not None:
        pool.close()
        pool.join()

    if uniprot is not None:
        uniprot.write_outputs(out_dir=args.out)

    if failed:
        logger.warning(f'{sum(failed.values())} junction(s) failed with an exception and wrote nothing (fate FAILURE '
                       f'(exception)). Their tracebacks are in {os.path.join(args.out, "logfile.log")}.')
        if failed['UniProtAPIError']:
            logger.warning(f'{failed["UniProtAPIError"]} of them failed on a UniProt API request, after retries. Rerun, '
                           f'or use --uniprot with a UniProt FASTA, which needs no network.')

    return True


def _translate_all(args, logger, rmats_results, gtf, genome, min_count, uniprot, pool) -> int:
    """
    Translate every junction of the selected splice types, in file order, in this process or in the pool

    :return: the junctions that failed with an exception, counted by the exception's type
    """
    failed = collections.Counter()

    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #
    for splice_type, rma in zip(["MXE", "SE", "RI", "A5SS", "A3SS"], rmats_results.rmats_list):

        # ---- Skip if the splice type is not in the list of splice types to run ----
        if splice_type not in args.splice_type:
            logger.info(f'Skipping {splice_type} because it was not specified in the -s argument.')
            continue

        # ---- Get a list of junctions from each splice type ----
        junctions = [Junction(**rma.iloc[i].to_dict()) for i in range(len(rma))]

        # ---- Worker processes translate; this process logs and writes each junction's result in file order, so
        # the output and the log are those of a single process (reports/2026-11-12-parallel-audit.md) ----
        if pool is not None:
            results = pool.imap(_translate_in_worker, junctions, chunksize=JUNCTIONS_PER_TASK)
            for jx, result in tqdm.tqdm(zip(junctions, results),
                                        total=len(junctions),
                                        desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
                                        ):
                logger.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(jx.junction_type,
                                                                                      jx.name,
                                                                                      jx.gene_symbol,
                                                                                      jx.gene_id,
                                                                                      ))
                fate = _write_worker_result(result, logger, uniprot)
                logger.info(fate)
                _count_failure(jx, fate, failed, logger)
            continue

        # ---- Single threaded for loop ----
        for jx in tqdm.tqdm(junctions,
                            total=len(junctions),
                            desc='Processing {0} Junctions'.format(rma.jxn_type[0]),
                            ):

            logger.info('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(jx.junction_type,
                                                                                  jx.name,
                                                                                  jx.gene_symbol,
                                                                                  jx.gene_id,
                                                                                  ))
            fate, records = _translate_safely(jx, gtf, genome, args, pred_bound=min_count, uniprot=uniprot)
            for output, suffix, seqrecord in records:
                h.write_seqrecord_to_fasta(seqrecord, output, suffix)
            logger.info(fate)
            _count_failure(jx, fate, failed, logger)

    return failed


def _count_failure(junction, fate, failed, logger) -> None:
    """ Count a junction that failed with an exception, by its type. The first is also shown on the console """
    if not fate.startswith(fates.EXCEPTION.split('.')[0]):
        return
    if not failed:
        logger.warning('Junction {0} {1} of gene {2} {3} failed: {4} The run continues; the traceback of each failed '
                       'junction is in the log.'.format(junction.junction_type, junction.name, junction.gene_symbol,
                                                        junction.gene_id, fate.split('. ', 1)[1].strip()))
    failed[fate.split('. ', 1)[1].split(':', 1)[0]] += 1


def _translate_safely(junction, gtf, genome, args, pred_bound, uniprot):
    """
    Translate one junction, keeping its records to be written once it completes, in this process or a worker.

    An exception fails this junction only: its traceback is logged (to the log file), its fate records it, and none
    of its records are written. 2026-11-12 one exception used to abort the run; see
    reports/2026-11-12-thresholds-runs-canonical-audit.md §3.3

    :return: the fate, and the records to write as (output directory, file suffix, SeqRecord)
    """
    h._collected = []
    try:
        fate = _translate_one(junction, gtf, genome, args, pred_bound=pred_bound, uniprot=uniprot)
        records = h._collected
    except Exception as e:
        logging.getLogger('jcast').error('Junction {0} {1} of gene {2} {3} failed:\n{4}'.format(
            junction.junction_type, junction.name, junction.gene_symbol, junction.gene_id,
            traceback.format_exc().rstrip()), extra={'jcast_console': False})
        fate, records = fates.EXCEPTION.format(type(e).__name__, ' '.join(str(e).split())), []
    finally:
        h._collected = None
    return fate, records


# Junctions sent to a worker at a time
JUNCTIONS_PER_TASK = 20

# tqdm's monitor thread would still be running when a later run in the same process forks its workers
tqdm.tqdm.monitor_interval = 0

# Shared with the worker processes, which are forked after it is filled
_WORKER = {}


def _number_of_processes(args, uniprot, logger) -> int:
    """
    The number of processes to translate in: -p, but one without --uniprot (the UniProt API calls and the local
    cache cannot be shared between processes) or where processes cannot be forked (Windows).
    """
    processes = getattr(args, 'processes', 1) or 1
    if processes > 1 and uniprot is None:
        logger.warning('-p/--processes > 1 needs --uniprot: canonicals from the UniProt API are fetched and cached '
                       'by one process only. Running in 1 process.')
        return 1
    if processes > 1 and 'fork' not in multiprocessing.get_all_start_methods():
        logger.warning('-p/--processes > 1 is not supported on this platform (no fork). Running in 1 process.')
        return 1
    return processes


class _CaptureLog(logging.Handler):
    """ In a worker: keep the junction's log records, to be logged by the main process """

    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        record.msg, record.args, record.exc_info, record.exc_text = record.getMessage(), None, None, None
        self.records.append(record)


def _init_worker() -> None:
    """ In each new worker: capture the jcast logger instead of writing to the log file and the console """
    logger = logging.getLogger('jcast')
    _WORKER['capture'] = _CaptureLog()
    logger.handlers = [_WORKER['capture']]


def _translate_in_worker(junction: Junction):
    """
    Translate one junction in a worker process.

    :return: the fate; the records to write, as (output directory, file suffix, SeqRecord); the log records; and
        the canonical this worker chose for the junction's gene, if not yet reported
    """
    _WORKER['capture'].records = []
    fate, records = _translate_safely(junction, _WORKER['gtf'], _WORKER['genome'], _WORKER['args'],
                                      pred_bound=_WORKER['pred_bound'], uniprot=_WORKER['uniprot'])

    chosen = {}
    uniprot = _WORKER['uniprot']
    if junction.gene_id in uniprot.chosen and junction.gene_id not in _WORKER['reported']:
        chosen[junction.gene_id] = uniprot.chosen[junction.gene_id]
        _WORKER['reported'].add(junction.gene_id)

    return fate, records, _WORKER['capture'].records, chosen


def _write_worker_result(result, logger, uniprot) -> str:
    """
    In the main process: log and write one junction's result from a worker, as a single process would have.

    :return: the fate
    """
    fate, records, log_records, chosen = result

    # a gene's canonical is logged when first chosen; another worker may have chosen it again
    for record in log_records:
        if getattr(record, 'jcast_canonical_for', None) not in uniprot.chosen:
            logger.handle(record)
    for gene_id, choice in chosen.items():
        uniprot.chosen.setdefault(gene_id, choice)

    for output, suffix, seqrecord in records:
        h.write_seqrecord_to_fasta(seqrecord, output, suffix)

    return fate


def _translate_one(junction: Junction,
                   gtf: ReadAnnotations,
                   genome: ReadGenome,
                   args: argparse.Namespace,
                   pred_bound: int,
                   uniprot: UniProtFasta = None,
                   ):
    """
    Get coordinates and translate one junction object.

    :param junction: Junction object
    :param gtf: ReadAnnotations object
    :param genome: ReadGenome object
    :param args: parsed arguments
    :param pred_bound: read count cutoff: a junction with fewer reads (Junction.min_jc) is skipped
    :param uniprot: UniProtFasta from --uniprot, or None to use the UniProt API
    :return: the junction's fate
    """

    sequence = Sequence(junction=junction, uniprot=uniprot)
    fate = _translate(junction, sequence, gtf, genome, args, pred_bound)

    # ---- -c: the gene's canonical is written whatever the junction's fate. This avoids not having any protein level
    # representation of a gene potentially in the proteome. 2026-11-12 the T1, T2 and T4 paths wrote it only when a
    # slice joined to it exactly, and the masked path before it was looked up (REVIEW.md MED-14;
    # reports/2026-11-12-thresholds-runs-canonical-audit.md §4.1) ----
    if args.canonical:
        if not sequence.canonical_looked_up:
            sequence.get_canonical_aa(gtf=gtf, genome_index=genome.genome)
        sequence.write_canonical(out_dir=args.out)

    return fate


def _translate(junction: Junction,
               sequence: Sequence,
               gtf: ReadAnnotations,
               genome: ReadGenome,
               args: argparse.Namespace,
               pred_bound: int,
               ) -> str:
    """
    The steps of _translate_one up to the junction's fate: filter, trim, translate, tier, and write the slices.

    :return: the junction's fate
    """

    # ---- Filters first: the read count, then the FDR. A skipped junction is not translated. 2026-11-12 the filters ran after slicing, translation and the canonical lookup, and
    # most of a run went on junctions it then skipped (reports/2026-11-12-log-audit.md §3). A junction is kept if its
    # read count is at least the cutoff; it used to need one read more
    # (reports/2026-11-12-thresholds-runs-canonical-audit.md §2) ----
    q_lo, q_hi = args.qvalue
    if junction.min_jc < pred_bound:
        skipped = fates.SKIPPED_LOW_COUNT
    elif not q_lo <= junction.fdr <= q_hi:
        skipped = fates.SKIPPED_P_VALUE
    else:
        skipped = None

    if skipped:
        return skipped

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
    # get nucleotide sequences of all slices using genome in memory
    # (anchor, alternative-1, alternative-2, downstream)
    # conjoin alternative exons to make slice 1 and 2,
    #
    sequence.make_slice_localgenome(genome.genome)

    # skip slices with any base outside ACGT: N (hard-masked or a gap) or another IUPAC code. The genome sequence is
    # uppercased on retrieval, so soft-masked bases are kept (reports/2026-11-11-masked-bases-audit.md)
    if not set(str(sequence.slice1_nt.seq) + str(sequence.slice2_nt.seq)) <= set('ACGT'):
        return fates.SKIPPED_MASKED

    # ---- Translate to peptides ----
    sequence.get_canonical_aa(gtf=gtf, genome_index=genome.genome)
    sequence.translate()

    # Discard junction if trimming to the coding region left both slices identical: the splice event
    # lies in the UTR and does not change the protein.
    if junction.event_outside_cds:
        return fates.SKIPPED_OUTSIDE_CDS

    # ---- Alternative termini: a slice that ends at an annotated stop codon, or that lacks the gene's start codon
    # but contains another annotated one, is complete at that end rather than a fragment.
    # 2026-11-10; see reports/2026-11-10-alt-termini-audit.md ----
    if sequence.route_alternative_termini():
        sequence.write_slices(out_dir=args.out,
                              suffix='T1',
                              slice_suffixes=sequence.alt_suffixes,
                              )
        alt_files = [s for s in sequence.alt_suffixes if s.startswith('T1_alt')]
        if alt_files:
            return fates.TIER1_ALT_TERMINUS.format(sequence.j.phase_label,
                                                   sequence.translated_phase,
                                                   ', '.join(alt_files),
                                                   )
        return fates.TIER1.format(sequence.j.phase_label,
                                  sequence.translated_phase,
                                  )

    # ---- Write the Tier 1 and Tier 2 results into fasta file ----
    if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) > 0:
        # Tier 1: both translated without stop codon, no frameshift
        if not sequence.frameshift:
            # Do a function like this to extend with fasta, and then write if necessary.
            # TODO: instead of using Uniprot, have the option to get the canonical exons from the GTF directly
            for slice_ in range(1, 3):
                sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                slice_has_ptc=False)
            sequence.write_slices(out_dir=args.out,
                                  suffix='T1',
                                  )
            return fates.TIER1.format(sequence.j.phase_label,
                                      sequence.translated_phase,
                                      )

        #
        # Tier 2: both translated without stop codon, but with one frameshift
        #
        elif sequence.frameshift:

            for slice_ in range(1, 3):
                sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                slice_has_ptc=False)

                # 2020-07-30 if slice runs into a frame shift,
                # allows the opportunity to stitch N-terminus only
                if [sequence.slice1_stitched, sequence.slice2_stitched][slice_-1] is None:
                    sequence.stitch_to_canonical_aa(slice_to_stitch=slice_,
                                                    slice_has_ptc=True)

            sequence.write_slices(
                out_dir=args.out,
                suffix='T2',
            )

            return fates.TIER2.format(sequence.j.phase_label,
                                      sequence.translated_phase,
                                      )

    #
    # 2026-11-10 Tier 3 (re-translating in a guessed frame when the annotated frame failed) was removed;
    # see reports/2026-11-10-phase-audit.md
    #
    # Tier 4: if sequence is still not good, do Tier 4: One of the two slices hits stop codon.
    # write out the slice if it is at least a certain proportion (params.ptc_threshold) as long as the long slice.
    #

    # force-translate through slice 2 if slice 2 hits PTC:
    if len(sequence.slice1_aa) > 0 and len(sequence.slice2_aa) == 0:
        forced_slice = 2
        sequence.stitch_to_canonical_aa(slice_to_stitch=1)
        sequence.translate_forced(slice_to_translate=forced_slice)

        if len(sequence.slice2_aa) / len(sequence.slice1_aa) >= params.ptc_threshold:
            sequence.stitch_to_canonical_aa(slice_to_stitch=2,
                                            slice_has_ptc=True)

        sequence.write_slices(out_dir=args.out,
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

        sequence.write_slices(out_dir=args.out,
                              suffix='T4',
                              )

        return fates.TIER4.format(forced_slice)

    #
    # if nothing works, write FAIL fate
    #
    else:
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
                        )
    parser.add_argument('genome',
                        help='path to genome file (optionally .gz)',
                        )

    # parser.add_argument('-n', '--num_threads', help='number of threads for concurrency [default: 6]',
    #                     default=6,
    #                     type=int)

    parser.add_argument('-o', '--out',
                        help='output directory, created if needed [default: out]',
                        default='out')

    parser.add_argument('--overwrite',
                        help='replace the output of an earlier run in the output directory; without it, a directory '
                             'holding JCAST output is refused [default: False]',
                        action='store_true',
                        default=False,
                        )

    parser.add_argument('--cache_dir',
                        help='directory for the GTF cache and the UniProt API cache [default: the user cache '
                             'directory, e.g. ~/Library/Caches/jcast or ~/.cache/jcast]',
                        default=None,
                        )

    parser.add_argument('-r', '--read',
                        help='the lowest read count for a junction to be translated: the smaller of its inclusion and '
                             'skipping junction counts, each summed over samples and replicates [default: 1, which '
                             'keeps every junction]',
                        default=1,
                        type=int,
                        )

    parser.add_argument('-m', '--model',
                        help='take the read count cutoff from a Gamma + Normal mixture model of the read counts, '
                             'instead of -r (used only if the fit fails) [default: False]',
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

    parser.add_argument('--uniprot',
                        help='UniProt FASTA (e.g. Swiss-Prot with isoforms for the species, optionally .gz) to take '
                             'canonical proteins from instead of the UniProt API; genes without an entry use the GTF. '
                             'Also writes jcast_uniprot_merged.fasta, psq_isoforms.tsv and uniprot_canonicals.tsv '
                             '[default: none, use the API]',
                        default=None,
                        )

    parser.add_argument('-p', '--processes',
                        help='number of processes to translate junctions in; the output is the same for any number. '
                             'More than 1 needs --uniprot (canonicals from the UniProt API are fetched and cached by '
                             'one process only), and macOS or Linux [default: 1]',
                        default=1,
                        type=int,
                        )

    parser.add_argument('--mask',
                        help='deprecated and ignored: slices containing any base other than A, C, G or T '
                             '(case-insensitive) are always skipped. Will be removed in the next release.',
                        default=None,
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

    if args.processes < 1:
        parser.error('-p/--processes must be at least 1')
    if args.processes > 1 and not args.uniprot:
        parser.error('-p/--processes > 1 needs --uniprot: without it, canonicals come from the UniProt API, whose '
                     'calls and local cache cannot be shared between processes')

    for path in [args.gtf_file, args.genome] + ([args.uniprot] if args.uniprot else []):
        if not os.path.isfile(path):
            parser.error(f'file not found: {path}')

    found = existing_outputs(args.out)
    if found and not args.overwrite:
        parser.error(f'output directory {args.out} already holds JCAST output ({", ".join(found)}); use --overwrite '
                     f'to replace it, or another -o')

    # Check if the rMATS directory contains the expected files
    for splice_type in args.splice_type:
        assert os.path.exists(os.path.join(args.rmats_folder, f'{splice_type}.MATS.JC.txt')), \
            'rMATS files not found, check directory.'

    # run the function in the argument
    args.func(args)
