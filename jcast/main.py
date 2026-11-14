# -*- coding: utf-8 -*-

"""jcast.main: Main function."""

import os
import collections
import datetime
import enum
import glob
import hashlib
import argparse
import json
import logging
import multiprocessing
import signal
import sys
import time
import traceback
from typing import Annotated, Optional

import tqdm
import typer

from jcast import params, fates, model, cache, outcomes
from jcast import helpers as h
from jcast.junctions import Junction, RmatsResults, rmats_path
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.sequences import Sequence
from jcast.uniprot import UniProtFasta, MERGED_FILE, ISOFORMS_FILE, CANONICALS_FILE
from jcast import __version__
from jcast.logger import get_logger, close_logger

# The files a run writes into its output directory. A directory holding any of them is used again only with
# --overwrite, which removes them first (reports/2026-11-12-thresholds-runs-canonical-audit.md §3.2)
OUTPUT_FILES = ['logfile.log', 'model.pdf', outcomes.JUNCTIONS_FILE, MERGED_FILE, ISOFORMS_FILE, CANONICALS_FILE]
OUTPUT_PATTERN = 'psq_*.fasta'


def existing_outputs(out_dir) -> list:
    """ JCAST output files already in a directory """
    if not os.path.isdir(out_dir):
        return []
    found = [f for f in OUTPUT_FILES if os.path.isfile(os.path.join(out_dir, f))]
    found += [os.path.basename(f) for f in glob.glob(os.path.join(glob.escape(out_dir), OUTPUT_PATTERN))]
    return sorted(found)


def run_jcast(args, progress=None) -> None:
    """
    Main function for jcast

    The output directory is created if needed. If it holds output of an earlier run, the run is refused unless
    args.overwrite is set, which removes that output first: appending to it mixed the two runs (REVIEW.md MED-12).

    :param args: parsed arguments
    :param progress: None, or a function that takes each step of the run as a dict (run --progress; _progress)
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
    logger = get_logger('jcast', args.out, verbose=getattr(args, 'verbose', False))
    try:
        return _run(args, logger, removed=found, progress=progress or _no_progress)
    except Stopped:
        logger.warning('Stopped by SIGTERM: the output is incomplete. Run again with --overwrite to replace it.')
        raise
    finally:
        close_logger(logger)


class Stopped(BaseException):
    """
    SIGTERM, during a run from the command line: the GUI's cancel. A BaseException, so that the junction it interrupts
    does not take it for its own failure (_translate_safely catches Exception)
    """


# the exit status of a run stopped by SIGTERM, as for a process the signal ends
STOPPED = 128 + signal.SIGTERM


def _stop(signum, frame):
    raise Stopped()


def _no_progress(step: dict) -> None:
    pass


def _progress(step: dict) -> None:
    """ run --progress: each step of the run as a line of JSON on stdout, for the GUI """
    print(json.dumps(step), flush=True)


def _describe(path, checksum: bool = True) -> str:
    """ a file's path and size, and its sha256 or its modification time, for the log """
    if checksum:
        digest = hashlib.sha256()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                digest.update(chunk)
        identity = f'sha256 {digest.hexdigest()}'
    else:
        identity = 'modified ' + datetime.datetime.fromtimestamp(os.path.getmtime(path)).isoformat(' ', 'seconds')
    return f'{os.path.abspath(path)} ({os.path.getsize(path):,} bytes, {identity})'


def _run(args, logger, removed=(), progress=_no_progress) -> bool:
    """
    The run, logging to logger; removed: the output of an earlier run that --overwrite removed; progress: takes each
    step: {'stage': <rmats, gtf, genome, uniprot, model, junctions, outputs or done>}, and for junctions also
    splice_type, done and total
    """
    started = time.monotonic()
    logger.info(f'JCAST {__version__}')
    logger.info(f'Arguments: {args}')
    if removed:
        logger.info(f'--overwrite: removed the output of an earlier run: {", ".join(removed)}.')

    # ---- Read in rMATS files: those of the -s splice types, and with -m, the others in the folder too, since the
    # model fits on every junction read (reports/2026-11-13-command-line-checks.md) ----
    progress({'stage': 'rmats'})
    rmats_results = RmatsResults(logger=logger,
                                 rmats_dir=args.rmats_folder,
                                 splice_types=args.splice_type,
                                 others_if_present=args.model,
                                 )

    # ---- The inputs. The GTF and the genome (GBs) get no checksum: it would read them again, about 20 s for GENCODE
    # and GRCh38 (reports/2026-11-12-log-audit.md §7) ----
    for splice_type in rmats_results.tables:
        unselected = '' if splice_type in args.splice_type else ' (not selected by -s; read for the -m model)'
        logger.info(f'rMATS {splice_type}{unselected}: {_describe(rmats_path(args.rmats_folder, splice_type))}')
    logger.info(f'GTF: {_describe(args.gtf_file, checksum=False)}')
    logger.info(f'Genome: {_describe(args.genome, checksum=False)}')

    # ---- Read the gtf file using  gtfpase then write as a pandas data frame. ----
    progress({'stage': 'gtf'})
    gtf = ReadAnnotations(logger=logger,
                          path=args.gtf_file,
                          )

    # ---- Read genome file into memory ----
    progress({'stage': 'genome'})
    genome = ReadGenome(logger=logger,
                        path=args.genome,
                        )

    # ---- Read the UniProt FASTA if supplied; otherwise canonicals come from the UniProt API ----
    if getattr(args, 'uniprot', None):
        progress({'stage': 'uniprot'})
    uniprot = UniProtFasta(logger=logger, path=args.uniprot) if getattr(args, 'uniprot', None) else None

    #
    # Model read count cutoff.
    # TODO: move this to a separate class
    #
    if args.model:

        logger.info('The -m flag is set. The modeled read count will override -r --read values.')
        progress({'stage': 'model'})

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

    table = outcomes.JunctionTable(args.out)
    translating = time.monotonic()
    try:
        failed = _translate_all(args, logger, rmats_results, gtf, genome, min_count, uniprot, pool, table, progress)
    except BaseException:
        if pool is not None:
            pool.terminate()
        raise
    finally:
        table.close()
    if pool is not None:
        pool.close()
        pool.join()
    translated = time.monotonic()

    progress({'stage': 'outputs'})
    if uniprot is not None:
        uniprot.write_outputs(out_dir=args.out)

    # ---- The end-of-run summary (reports/2026-11-12-log-audit.md §4) ----
    logger.info('\n'.join(['Summary.'] + table.summary(args.splice_type, args.out)))
    logger.info(f'Time: {_duration(time.monotonic() - started)} in all, of which {_duration(translated - translating)} '
                f'translating junctions.')

    if failed:
        logger.warning(f'{sum(failed.values())} junction(s) failed with an exception and wrote nothing (fate FAILURE '
                       f'(exception)). Their tracebacks are in {os.path.join(args.out, "logfile.log")}.')
        if failed['UniProtAPIError']:
            logger.warning(f'{failed["UniProtAPIError"]} of them failed on a UniProt API request, after retries. Rerun, '
                           f'or use --uniprot with a UniProt FASTA, which needs no network.')

    progress({'stage': 'done'})
    return True


def _duration(seconds: float) -> str:
    """ e.g. 4 min 05 s """
    minutes, seconds = divmod(round(seconds), 60)
    return f'{minutes} min {seconds:02d} s' if minutes else f'{seconds} s'


def _translate_all(args, logger, rmats_results, gtf, genome, min_count, uniprot, pool, table,
                   progress=_no_progress) -> int:
    """
    Translate every junction of the selected splice types, in file order, in this process or in the pool; write
    their records, and their rows of junctions.tsv. Progress is shown by tqdm, or with run --progress, reported to
    progress

    :return: the junctions that failed with an exception, counted by the exception's type
    """
    failed = collections.Counter()

    #
    # Main loop through every line of each of the five rMATS files to make junction object, then translate them
    #
    for splice_type in outcomes.SPLICE_TYPES:

        # ---- Skip if the splice type is not in the list of splice types to run ----
        if splice_type not in args.splice_type:
            logger.info(f'Skipping {splice_type} because it was not specified in the -s argument.')
            continue
        rma = rmats_results.tables[splice_type]

        # ---- Get a list of junctions from each splice type ----
        junctions = [Junction(**rma.iloc[i].to_dict()) for i in range(len(rma))]

        # ---- Worker processes translate; this process logs and writes each junction's result in file order, so
        # the output and the log are those of a single process (reports/2026-11-12-parallel-audit.md) ----
        if pool is not None:
            results = pool.imap(_translate_in_worker, junctions, chunksize=JUNCTIONS_PER_TASK)
            for jx, result in _each(zip(junctions, results), len(junctions), splice_type, progress):
                _log_start(jx, logger)
                fate, records, found = _write_worker_result(result, logger, uniprot)
                _finish(jx, fate, records, found, table, failed, logger)
            continue

        # ---- Single threaded for loop ----
        for jx in _each(junctions, len(junctions), splice_type, progress):
            _log_start(jx, logger)
            fate, records, found = _translate_safely(jx, gtf, genome, args, pred_bound=min_count, uniprot=uniprot)
            _finish(jx, fate, records, found, table, failed, logger)

    return failed


def _each(items, total: int, splice_type: str, progress):
    """
    The items, with a progress bar on the console; or with run --progress, about each 1 % of them reported to progress
    instead, as {'stage': 'junctions', 'splice_type', 'done', 'total'}
    """
    if progress is _no_progress:
        yield from tqdm.tqdm(items, total=total, desc='Processing {0} Junctions'.format(splice_type))
        return
    step = max(1, total // 100)
    progress({'stage': 'junctions', 'splice_type': splice_type, 'done': 0, 'total': total})
    for done, item in enumerate(items, 1):
        yield item
        if done % step == 0 or done == total:
            progress({'stage': 'junctions', 'splice_type': splice_type, 'done': done, 'total': total})


def _log_start(junction, logger) -> None:
    logger.debug('>>>>>> Doing {0} junction {1} for gene {2} {3}'.format(junction.junction_type,
                                                                          junction.name,
                                                                          junction.gene_symbol,
                                                                          junction.gene_id,
                                                                          ))


def _finish(junction, fate, records, found, table, failed, logger) -> None:
    """ In the main process: write the junction's records, log its fate, and add its row to junctions.tsv """
    new_records = 0
    for output, suffix, seqrecord in records:
        new_records += h.write_seqrecord_to_fasta(seqrecord, output, suffix)
    logger.debug(fate.rstrip())
    table.add(junction, fate, found, [suffix for _, suffix, _ in records], new_records)
    _count_failure(junction, fate, failed, logger)


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

    :return: the fate; the records to write as (output directory, file suffix, SeqRecord); and what translating
        found, for junctions.tsv (outcomes.outcome; {} after an exception)
    """
    h._collected = []
    try:
        fate = _translate_one(junction, gtf, genome, args, pred_bound=pred_bound, uniprot=uniprot)
        records, found = h._collected, outcomes.outcome(junction)
    except Exception as e:
        logging.getLogger('jcast').error('Junction {0} {1} of gene {2} {3} failed:\n{4}'.format(
            junction.junction_type, junction.name, junction.gene_symbol, junction.gene_id,
            traceback.format_exc().rstrip()), extra={'jcast_console': False})
        fate, records, found = fates.EXCEPTION.format(type(e).__name__, ' '.join(str(e).split())), [], {}
    finally:
        h._collected = None
    return fate, records, found


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
    """
    In each new worker: capture the jcast logger instead of writing to the log file and the console. SIGTERM, which
    pool.terminate sends, ends the worker, as before the command line handled it (Stopped)
    """
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    logger = logging.getLogger('jcast')
    _WORKER['capture'] = _CaptureLog()
    logger.handlers = [_WORKER['capture']]


def _translate_in_worker(junction: Junction):
    """
    Translate one junction in a worker process.

    :return: the fate; the records to write, as (output directory, file suffix, SeqRecord); what translating found,
        for junctions.tsv; the log records; and the canonical this worker chose for the junction's gene, if not yet
        reported
    """
    _WORKER['capture'].records = []
    fate, records, found = _translate_safely(junction, _WORKER['gtf'], _WORKER['genome'], _WORKER['args'],
                                             pred_bound=_WORKER['pred_bound'], uniprot=_WORKER['uniprot'])

    chosen = {}
    uniprot = _WORKER['uniprot']
    if junction.gene_id in uniprot.chosen and junction.gene_id not in _WORKER['reported']:
        chosen[junction.gene_id] = uniprot.chosen[junction.gene_id]
        _WORKER['reported'].add(junction.gene_id)

    return fate, records, found, _WORKER['capture'].records, chosen


def _write_worker_result(result, logger, uniprot) -> tuple:
    """
    In the main process: log one junction's result from a worker, as a single process would have.

    :return: the fate, the records to write, and what translating found
    """
    fate, records, found, log_records, chosen = result

    # a gene's canonical is logged when first chosen; another worker may have chosen it again
    for record in log_records:
        if getattr(record, 'jcast_canonical_for', None) not in uniprot.chosen:
            logger.handle(record)
    for gene_id, choice in chosen.items():
        uniprot.chosen.setdefault(gene_id, choice)

    return fate, records, found


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

    junction.canonical = sequence.canonical_source
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
    skipped = fates.filtered(junction.min_jc, junction.fdr, pred_bound, args.qvalue)
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




# ---- The command line: `jcast run` and `jcast gui`, on Typer. 2026-11-13 argparse, with one command; see
# reports/2026-11-13-gui-cli-scope.md ----

app = typer.Typer(add_completion=False,
                  pretty_exceptions_enable=False,
                  context_settings={'help_option_names': ['-h', '--help']},
                  help='JCAST retrieves transcript splice junctions and translates them into amino acid sequences. '
                       'jcast RMATS_FOLDER GTF_FILE GENOME [OPTIONS] is short for jcast run.',
                  )

COMMANDS = ['run', 'gui']

# -s/--splice_type's choices
SpliceType = enum.Enum('SpliceType', [(t, t) for t in outcomes.SPLICE_TYPES], type=str)


def _show_version(value: bool) -> None:
    if value:
        typer.echo(f'jcast {__version__}')
        raise typer.Exit()


@app.callback()
def _jcast(version: Annotated[bool, typer.Option('--version', help='show the version and exit', is_eager=True,
                                                 callback=_show_version)] = False):
    pass


@app.command(no_args_is_help=True)
def run(ctx: typer.Context,
        rmats_folder: Annotated[str, typer.Argument(metavar='RMATS_FOLDER', help='path to folder storing rMATS output')],
        gtf_file: Annotated[str, typer.Argument(metavar='GTF_FILE', help='path to Ensembl/Gencode GTF file')],
        genome: Annotated[str, typer.Argument(metavar='GENOME', help='path to genome FASTA (optionally .gz)')],
        out: Annotated[str, typer.Option('-o', '--out', rich_help_panel='Output', metavar='DIR',
                                         help='output directory, created if needed')] = 'out',
        overwrite: Annotated[bool, typer.Option('--overwrite', rich_help_panel='Output',
                                                help='replace the output of an earlier run in the output directory; '
                                                     'without it, a directory holding JCAST output is refused')] = False,
        cache_dir: Annotated[Optional[str], typer.Option('--cache_dir', rich_help_panel='Output', metavar='DIR',
                                                         help='directory for the GTF cache and the UniProt API cache',
                                                         show_default='the user cache directory, e.g. '
                                                                      '~/Library/Caches/jcast or ~/.cache/jcast'
                                                         )] = None,
        read: Annotated[int, typer.Option('-r', '--read', rich_help_panel='Filters',
                                          help='the lowest read count for a junction to be translated: the smaller of '
                                               'its inclusion and skipping junction counts, each summed over samples '
                                               'and replicates. 1 keeps every junction')] = 1,
        model: Annotated[bool, typer.Option('-m', '--model', rich_help_panel='Filters',
                                            help='take the read count cutoff from a Gamma + Normal mixture model of the '
                                                 'read counts of every rMATS file in the folder (also those -s leaves '
                                                 'out), instead of -r (used only if the fit fails)')] = False,
        qvalue: Annotated[tuple[float, float], typer.Option('-q', '--qvalue', rich_help_panel='Filters',
                                                            metavar='Q_LO Q_HI',
                                                            help="keep junctions whose rMATS FDR is in this range"
                                                            )] = (0, 1),
        splice_type: Annotated[Optional[list[SpliceType]], typer.Option('-s', '--splice_type', rich_help_panel='Filters',
                                                                        metavar='MXE|SE|RI|A5SS|A3SS',
                                                                        help='splice types to translate (-s SE MXE, or '
                                                                             '-s SE -s MXE); only their rMATS files '
                                                                             'need to be in the folder',
                                                                        show_default='MXE SE RI A5SS A3SS')] = None,
        canonical: Annotated[bool, typer.Option('-c', '--canonical', rich_help_panel='Canonicals',
                                                help="write the canonical protein of every gene with a junction to "
                                                     "psq_canonical.fasta, whatever the junction's fate")] = False,
        uniprot: Annotated[Optional[str], typer.Option('--uniprot', rich_help_panel='Canonicals', metavar='FASTA',
                                                       help='UniProt FASTA (e.g. Swiss-Prot with isoforms for the '
                                                            'species, optionally .gz) to take canonical proteins from '
                                                            'instead of the UniProt API; genes without an entry use '
                                                            'the GTF. Also writes jcast_uniprot_merged.fasta, '
                                                            'psq_isoforms.tsv and uniprot_canonicals.tsv',
                                                       show_default='none, use the API')] = None,
        processes: Annotated[int, typer.Option('-p', '--processes', rich_help_panel='The run',
                                               help='number of processes to translate junctions in; the output is the '
                                                    'same for any number. More than 1 needs --uniprot (canonicals from '
                                                    'the UniProt API are fetched and cached by one process only), and '
                                                    'macOS or Linux')] = 1,
        verbose: Annotated[bool, typer.Option('-v', '--verbose', rich_help_panel='The run',
                                              help="also log each junction's steps and its sequences, at DEBUG level "
                                                   "(large: 300 MB for a whole-genome -m run); without it, the log has "
                                                   "the run's settings, warnings, errors and a summary, and each "
                                                   "junction's outcome is in junctions.tsv")] = False,
        progress: Annotated[bool, typer.Option('--progress', hidden=True,
                                               help='write each step of the run as a line of JSON on stdout, for '
                                                    'the GUI')] = False,
        ) -> None:
    """ Translate the splice junctions in rMATS output into protein sequences. """
    # the arguments as argparse gave them before PR 18, in its order, so the log's Arguments line keeps its form
    args = argparse.Namespace(rmats_folder=rmats_folder, gtf_file=gtf_file, genome=genome, out=out, overwrite=overwrite,
                              cache_dir=cache_dir, read=read, model=model, canonical=canonical, uniprot=uniprot,
                              processes=processes, verbose=verbose,
                              splice_type=[t.value for t in splice_type] if splice_type else set(outcomes.SPLICE_TYPES),
                              qvalue=list(qvalue),
                              )
    problems = validate(args)
    if problems:
        ctx.fail('\n'.join(message for _, message in problems))

    # SIGTERM (the GUI's cancel) stops the run cleanly: the worker processes ended, the log says so
    previous = signal.signal(signal.SIGTERM, _stop)
    try:
        run_jcast(args, progress=_progress if progress else None)
    except Stopped:
        raise typer.Exit(STOPPED)
    finally:
        signal.signal(signal.SIGTERM, previous)


@app.command()
def gui() -> None:
    """ Open the graphical interface (a preview; it needs PySide6: pip install 'jcast\\[gui]'). """
    from jcast.gui import launch_gui

    raise typer.Exit(launch_gui())


def validate(args) -> list:
    """
    What keeps a run from starting, for the command line to report as a usage error and the GUI next to its fields

    :param args: the run's arguments, as run_jcast takes them
    :return: (argument, message) for each problem, in the order the command line checked them before PR 18; empty if
        the run can start
    """
    problems = []
    q_lo, q_hi = args.qvalue
    if not 0 <= q_lo <= 1:
        problems.append(('qvalue', 'qvalue lower bound must be between 0 and 1'))
    if not 0 <= q_hi <= 1:
        problems.append(('qvalue', 'qvalue upper bound must be between 0 and 1'))
    if q_lo >= q_hi:
        problems.append(('qvalue', 'qvalue lower bound must be smaller than upper bound'))

    if args.processes < 1:
        problems.append(('processes', '-p/--processes must be at least 1'))
    if args.processes > 1 and not args.uniprot:
        problems.append(('processes', '-p/--processes > 1 needs --uniprot: without it, canonicals come from the UniProt '
                                      'API, whose calls and local cache cannot be shared between processes'))

    for name in ['gtf_file', 'genome'] + (['uniprot'] if args.uniprot else []):
        if not os.path.isfile(getattr(args, name)):
            problems.append((name, f'file not found: {getattr(args, name)}'))

    # only the -s splice types' rMATS files must exist (REVIEW.md MED-8)
    if not os.path.isdir(args.rmats_folder):
        problems.append(('rmats_folder', f'rMATS folder not found: {args.rmats_folder}'))
    else:
        missing = [rmats_path(args.rmats_folder, t) for t in outcomes.SPLICE_TYPES
                   if t in args.splice_type and not os.path.isfile(rmats_path(args.rmats_folder, t))]
        if missing:
            problems.append(('splice_type', f'rMATS file not found: {", ".join(missing)}; select the splice types the '
                                            f'folder has with -s'))

    found = existing_outputs(args.out)
    if found and not args.overwrite:
        problems.append(('out', f'output directory {args.out} already holds JCAST output ({", ".join(found)}); use '
                                f'--overwrite to replace it, or another -o'))
    return problems


def _as_before(argv: list) -> list:
    """
    The command line as Typer takes it, from the forms accepted before PR 18: no command name (jcast RMATS GTF GENOME
    runs `jcast run`), --gui, no arguments (the help), and several splice types after one -s (-s SE MXE, which Click
    cannot parse: an option takes a fixed number of values), given each its own -s
    """
    if not argv:
        return ['--help']
    if '--gui' in argv:
        return ['gui']
    if argv[0] not in COMMANDS + ['-h', '--help', '--version']:
        argv = ['run'] + argv

    words, option = [], None    # option: -s or --splice_type, while the words after it are splice types
    for i, word in enumerate(argv):
        if word == '--':
            return words + argv[i:]
        if option and word in outcomes.SPLICE_TYPES:
            if words[-1] != option:
                words.append(option)
        else:
            option = None
        words.append(word)
        if word in ('-s', '--splice_type'):
            option = word
    return words


def main(argv: list = None) -> None:
    """
    The jcast command. It returns when the command succeeds, and exits with its status otherwise: 2 for a usage error
    """
    try:
        app(args=_as_before(sys.argv[1:] if argv is None else list(argv)), prog_name='jcast')
    except SystemExit as e:
        if e.code:
            raise
