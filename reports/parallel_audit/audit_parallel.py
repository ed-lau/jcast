# -*- coding: utf-8 -*-

"""Audit: where JCAST spends its time on a whole-genome run, and what parallelizing it would gain.

On a seeded random sample of junctions from a real rMATS run (whole-genome GTF and genome, --uniprot), measures:

    1. load       time to read the GTF, genome, rMATS files and UniProt FASTA
    2. profile    the current code, one junction at a time (cProfile), split into: whole-GTF queries, FASTA writes,
                  the rest
    3. variants   the same sample with candidate changes patched in at runtime (jcast code is not modified):
                    gene_index   each junction sees only its gene's GTF rows (a per-gene index built once)
                    memory_dedup FASTA de-duplication against an in-memory set instead of re-reading the file
                  each checked to give byte-identical output files and identical fates
    4. write cost write_seqrecord_to_fasta against output files of growing size (the current re-read)
    5. parallel   gene_index + memory_dedup with a process pool: workers translate and return (fate, records);
                  the main process writes records in junction order, so output is identical to the serial run.
                  Timed for several worker counts on a larger sample.

Writes, into reports/parallel_audit/<run>/: summary.txt, profile.txt, variants.tsv, parallel.tsv.

Usage (from the repository root):
    python reports/parallel_audit/audit_parallel.py --rmats <rmats dir> --gtf <gtf> --genome <fa> --uniprot <fasta> \\
        [--sample 300] [--parallel_sample 3000] [--workers 1 4 8 16] [--out_name <run>]
"""

import argparse
import collections
import cProfile
import hashlib
import io
import logging
import multiprocessing as mp
import os
import pstats
import random
import shutil
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1].parent
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))

from Bio import SeqIO  # noqa: E402

from jcast import helpers as h, main  # noqa: E402
from jcast.annots import ReadAnnotations, ReadGenome  # noqa: E402
from jcast.junctions import Junction, RmatsResults  # noqa: E402
from jcast.uniprot import UniProtFasta  # noqa: E402

TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']

# shared with forked workers
G = {}


class GeneGTF(object):
    """ A GTF holding one gene's rows: every GTF query in jcast first filters by the junction's gene id """
    def __init__(self, annot):
        self.annot = annot


def args_for(out_dir):
    return argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])


def sample(rmats, n, seed):
    rows = [(t, i) for t, rma in zip(TYPES, rmats.rmats_list) for i in range(len(rma))]
    random.Random(seed).shuffle(rows)
    return sorted(rows[:n], key=lambda r: (TYPES.index(r[0]), r[1]))


def junction(rmats, t, i):
    return Junction(**rmats.rmats_list[TYPES.index(t)].iloc[i].to_dict())


def digest(out_dir):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in sorted(Path(out_dir).glob('psq_*.fasta'))}


def run_serial(rows, out_dir, gene_index=False, memory_dedup=False, profile=None):
    """ the current loop, optionally with the variants patched in """
    out_dir.mkdir(parents=True)
    G['uniprot'].chosen = {}    # every run starts with no canonical chosen yet
    original = h.write_seqrecord_to_fasta
    if memory_dedup:
        seen = collections.defaultdict(set)

        def write(seqrecord, output, suffix):
            if str(seqrecord.seq) in seen[suffix]:
                return True
            seen[suffix].add(str(seqrecord.seq))
            with open(os.path.join(output, 'psq_' + suffix + '.fasta'), 'a') as fh:
                SeqIO.write(seqrecord, fh, 'fasta')
            return True
        h.write_seqrecord_to_fasta = write
    fates = []
    t0 = time.perf_counter()
    try:
        if profile:
            profile.enable()
        for t, i in rows:
            jx = junction(G['rmats'], t, i)
            gtf = GeneGTF(G['by_gene'].get(jx.gene_id, G['empty'])) if gene_index else G['gtf']
            fates.append(main._translate_one(jx, gtf, G['genome'], args_for(out_dir), pred_bound=1,
                                             uniprot=G['uniprot']))
        if profile:
            profile.disable()
    finally:
        h.write_seqrecord_to_fasta = original
    return time.perf_counter() - t0, fates


def _worker(task):
    """ translate a chunk of junctions; return each fate and the records it would have written, in order """
    out = []

    def collect(seqrecord, output, suffix):
        records.append((suffix, seqrecord))
        return True
    h.write_seqrecord_to_fasta = collect
    for t, i in task:
        records = []
        jx = junction(G['rmats'], t, i)
        fate = main._translate_one(jx, GeneGTF(G['by_gene'].get(jx.gene_id, G['empty'])), G['genome'],
                                   args_for('unused'), pred_bound=1, uniprot=G['uniprot'])
        out.append((fate, records))
    return out


def run_parallel(rows, out_dir, workers, chunk):
    out_dir.mkdir(parents=True)
    G['uniprot'].chosen = {}
    seen = collections.defaultdict(set)
    fates = []
    tasks = [rows[k:k + chunk] for k in range(0, len(rows), chunk)]
    t0 = time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        for result in pool.imap(_worker, tasks):
            for fate, records in result:
                fates.append(fate)
                for suffix, r in records:
                    if str(r.seq) not in seen[suffix]:
                        seen[suffix].add(str(r.seq))
                        with open(out_dir / f'psq_{suffix}.fasta', 'a') as fh:
                            SeqIO.write(r, fh, 'fasta')
    return time.perf_counter() - t0, fates


def rss_gb():
    import psutil
    p = psutil.Process()
    return (p.memory_info().rss + sum(c.memory_info().rss for c in p.children(recursive=True))) / 1e9


def main_(a):
    import psutil
    logger = logging.getLogger('jcast')
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    logger.setLevel(logging.WARNING)

    res = HERE / a.out_name
    res.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='parallel_audit_', dir=a.work))
    out = [f'inputs: rmats {a.rmats}; gtf {a.gtf}; genome {a.genome}; uniprot {a.uniprot}\n'
           f'machine: {os.cpu_count()} CPUs, {psutil.virtual_memory().total / 1e9:.0f} GB; '
           f'load average at start {os.getloadavg()[0]:.1f}\n']

    # 1. load
    t = {}
    t0 = time.perf_counter()
    G['gtf'] = ReadAnnotations(logger=logger, path=a.gtf)
    t['GTF (from .cached)'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    G['genome'] = ReadGenome(logger=logger, path=a.genome)
    t['genome'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    G['rmats'] = RmatsResults(logger=logger, rmats_dir=a.rmats)
    t['rMATS'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    G['uniprot'] = UniProtFasta(logger=logger, path=a.uniprot)
    t['UniProt FASTA'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    G['by_gene'] = {g: d for g, d in G['gtf'].annot.groupby('gene_id', sort=False)}
    G['empty'] = G['gtf'].annot.iloc[:0]
    t['per-gene GTF index (variant)'] = time.perf_counter() - t0
    n_junctions = sum(len(r) for r in G['rmats'].rmats_list)
    out.append('1. load (s)\n' + ''.join(f'  {k:32s}{v:8.1f}\n' for k, v in t.items()) +
               f'  GTF rows {len(G["gtf"].annot):,}; genes {len(G["by_gene"]):,}; junctions {n_junctions:,}; '
               f'memory after load {rss_gb():.1f} GB\n')

    # 2. profile, current code
    rows = sample(G['rmats'], a.sample, a.seed)
    prof = cProfile.Profile()
    base_t, base_fates = run_serial(rows, work / 'current', profile=prof)
    s = io.StringIO()
    st = pstats.Stats(prof, stream=s)
    st.sort_stats('cumulative').print_stats(40)
    (res / 'profile.txt').write_text(s.getvalue())
    stats = st.stats

    def cum(name_part, file_part=''):
        return sum(v[3] for k, v in stats.items() if k[2] == name_part and file_part in k[0])
    query_t = cum('query', 'pandas')
    write_t = cum('write_seqrecord_to_fasta')
    gp_t = cum('gene_proteins')
    prof_t = sum(v[2] for v in stats.values())
    out.append(f'2. current code, {len(rows)} sampled junctions (seed {a.seed}): {base_t:.1f} s, '
               f'{len(rows) / base_t:.2f} junctions/s (unprofiled rate below)\n'
               f'  under the profiler: total {prof_t:.1f} s; DataFrame.query {query_t:.1f} s '
               f'({100 * query_t / prof_t:.0f}%); write_seqrecord_to_fasta {write_t:.1f} s '
               f'({100 * write_t / prof_t:.0f}%); of which gene_proteins (--uniprot mapping) {gp_t:.1f} s\n')

    # 3. variants
    lines = ['variant\tseconds\tjunctions_per_s\tspeedup\tidentical_output\tidentical_fates']
    ref_t, ref_fates = run_serial(rows, work / 'current_unprofiled')
    ref = digest(work / 'current_unprofiled')
    assert ref_fates == base_fates
    for name, kw in [('current', {}), ('gene_index', {'gene_index': True}),
                     ('memory_dedup', {'memory_dedup': True}),
                     ('gene_index+memory_dedup', {'gene_index': True, 'memory_dedup': True})]:
        if name == 'current':
            secs, fates, dg = ref_t, ref_fates, ref
        else:
            secs, fates = run_serial(rows, work / name, **kw)
            dg = digest(work / name)
        lines.append(f'{name}\t{secs:.1f}\t{len(rows) / secs:.2f}\t{ref_t / secs:.1f}\t{dg == ref}\t{fates == ref_fates}')
    (res / 'variants.tsv').write_text('\n'.join(lines) + '\n')
    out.append('3. variants, same sample\n' + ''.join(f'  {x}\n' for x in lines))

    # 4. write cost against a growing file
    rec = next(SeqIO.parse(a.write_probe, 'fasta')) if a.write_probe else None
    if rec is not None:
        lines = []
        probe = work / 'probe'
        probe.mkdir()
        records = list(SeqIO.parse(a.write_probe, 'fasta'))
        for n in [1000, 10000, len(records)]:
            n = min(n, len(records))
            with open(probe / 'psq_X.fasta', 'w') as fh:
                SeqIO.write(records[:n], fh, 'fasta')
            size = (probe / 'psq_X.fasta').stat().st_size / 1e6
            t0 = time.perf_counter()
            for _ in range(5):
                h.write_seqrecord_to_fasta(records[0], str(probe), 'X')   # already present: read only
            lines.append(f'  {n:7,d} records ({size:5.1f} MB): {(time.perf_counter() - t0) / 5 * 1000:7.0f} ms per write')
        out.append(f'4. one write_seqrecord_to_fasta call (current code) against a file of\n' + '\n'.join(lines) + '\n')

    # 5. parallel
    rows = sample(G['rmats'], a.parallel_sample, a.seed + 1)
    lines = ['workers\tseconds\tjunctions_per_s\tspeedup_vs_1\tidentical_to_serial\tpeak_memory_gb']
    s_t, s_fates = run_serial(rows, work / 'p_serial', gene_index=True, memory_dedup=True)
    s_dg = digest(work / 'p_serial')
    one = None
    for w in a.workers:
        import threading
        peak = [rss_gb()]
        stop = threading.Event()

        def watch():
            while not stop.is_set():
                peak[0] = max(peak[0], rss_gb())
                time.sleep(0.5)
        th = threading.Thread(target=watch, daemon=True)
        th.start()
        secs, fates = run_parallel(rows, work / f'p_{w}', w, a.chunk)
        stop.set()
        th.join()
        one = one or secs
        ok = digest(work / f'p_{w}') == s_dg and fates == s_fates
        lines.append(f'{w}\t{secs:.1f}\t{len(rows) / secs:.1f}\t{one / secs:.1f}\t{ok}\t{peak[0]:.1f}')
    (res / 'parallel.tsv').write_text('\n'.join(lines) + '\n')
    out.append(f'5. parallel (gene_index + memory_dedup; fork; chunks of {a.chunk}), {len(rows)} sampled junctions '
               f'(seed {a.seed + 1}); serial with the same variants: {s_t:.1f} s ({len(rows) / s_t:.1f}/s)\n' +
               ''.join(f'  {x}\n' for x in lines))
    out.append(f'fates in the parallel sample: {dict(collections.Counter(f.split(".")[0] for f in s_fates))}\n')
    (res / 'summary.txt').write_text('\n'.join(out))
    print((res / 'summary.txt').read_text())
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--rmats', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--uniprot', required=True)
    p.add_argument('--write_probe', help='an existing large psq_T1.fasta, to time the current write against')
    p.add_argument('--sample', type=int, default=300)
    p.add_argument('--parallel_sample', type=int, default=3000)
    p.add_argument('--workers', type=int, nargs='+', default=[1, 4, 8, 16])
    p.add_argument('--chunk', type=int, default=20)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--work', default=None, help='scratch directory for outputs')
    p.add_argument('--out_name', default='genome')
    main_(p.parse_args())
