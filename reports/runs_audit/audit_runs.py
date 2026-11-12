# -*- coding: utf-8 -*-

"""Audit of PRs 6, 7 and 8 (reports/2026-11-12-thresholds-runs-canonical-audit.md): read-count thresholds, the
robustness of real runs, and canonicals under -c.

Runs the working tree in this process; to reproduce the evidence exactly, check out the commit named in the report.
Three parts, each writing into --out (default reports/runs_audit/audit):

    thresholds  (PR 6)  thresholds.tsv: per rMATS set (chr15, the 12 bundled model sets, and optionally a
                        whole-genome set), the read counts and how many junctions each cutoff rule keeps.
                        threshold_fates.tsv: chr15 fates and output records under each rule.
    runs        (PR 7)  runs.txt: chr15 runs through run_jcast with --uniprot showing each failure mode: a missing
                        output directory, a rerun into the same directory, an exception in one junction (with and
                        without -p), logger handlers across runs, and where the two caches are written.
    canonical   (PR 8)  canonical.tsv: genes with a canonical but no record under -c (chr15 GTF-only and --uniprot,
                        and optionally a whole-genome run's output); canonical.txt: the masked path, the UniProt API
                        calls made per path, and what a network failure or a refused request does to a run.

summary.txt collects the numbers quoted in the report. The chr15 runs need no network: the UniProt API is replaced
by a stub that serves the local API cache (cache/uniprot-cache.db, keys without the gene version) and an empty
answer otherwise.

Usage (from the repository root):
    python reports/runs_audit/audit_runs.py [--rmats <whole-genome rMATS dir> --jcast_out <its jcast -c --uniprot
        output dir> --uniprot <the UniProt FASTA it used>]
"""

import argparse
import collections
import contextlib
import hashlib
import logging
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import traceback
import warnings
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
DATA = REPO / 'tests' / 'data'
RMATS = DATA / 'rmats'
GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'
API_CACHE = REPO / 'cache' / 'uniprot-cache.db'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']
FAULT_AT = 500    # the junction made to raise in the exception runs, by its number in run order

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import main, model, params
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.junctions import Junction, RmatsResults
from jcast.sequences import Sequence
from jcast.uniprot import UniProtFasta

LOG = logging.getLogger('jcast-audit')


@contextlib.contextmanager
def chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


@contextlib.contextmanager
def patched(obj, name, value):
    old = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, old)


def label(fate):
    """ fate class, keeping the reason for skips """
    head = fate.split('.')
    return f'{head[0]}: {head[1].strip()}' if head[0] == 'SKIPPED' else head[0]


def records(out_dir):
    """ psq file -> number of records """
    return {p.name: sum(1 for _ in SeqIO.parse(p, 'fasta')) for p in sorted(Path(out_dir).glob('psq_*.fasta'))}


def digests(out_dir):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
            for p in sorted(Path(out_dir).iterdir()) if p.is_file() and p.name != 'logfile.log'}


# ---- read counts (PR 6) ------------------------------------------------------------------------------------------

def _total(a, b):
    """ summed counts of two rMATS count fields (comma-separated replicates), as in Junction.sum_sjc """
    return sum(int(x) for x in str(a).split(',')) + sum(int(x) for x in str(b).split(','))


def read_counts(rmats_dir):
    """ per junction: (min of summed IJC and SJC, which is Junction.sum_sjc; summed SJC) """
    import pandas as pd
    out = []
    for t in TYPES:
        df = pd.read_table(Path(rmats_dir) / f'{t}.MATS.JC.txt', low_memory=False)
        for s1, s2, i1, i2 in zip(df.SJC_SAMPLE_1, df.SJC_SAMPLE_2, df.IJC_SAMPLE_1, df.IJC_SAMPLE_2):
            try:
                sjc = _total(s1, s2)
            except ValueError:
                sjc = 0
            try:
                ijc = _total(i1, i2)
            except ValueError:
                ijc = 0
            out.append((min(sjc, ijc), sjc))
    return out


def threshold_row(name, counts):
    import numpy as np
    m = np.array([c[0] for c in counts])
    sjc = np.array([c[1] for c in counts])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, _, mc = model.general_mixture_model(np.array([[x + 1] for x in m]))
    row = dict(dataset=name, junctions=len(m), min_count_0=int((m == 0).sum()), min_count_1=int((m == 1).sum()),
               kept_r1_now=int((m > 1).sum()), kept_r1_as_documented=int((m >= 1).sum()),
               kept_sjc_only_r1_now=int((sjc > 1).sum()), sjc_ge2_but_min_lt2=int(((sjc > 1) & (m <= 1)).sum()),
               model_min_count=mc if mc is not None else 'degenerate')
    if mc is not None:
        row.update(kept_model_now=int((m > mc).sum()), kept_model_boundary=int((m >= mc - 1).sum()),
                   lost_to_offset=int(((m >= mc - 1) & (m <= mc)).sum()))
    return row


def part_thresholds(a, ctx, out, summary):
    rows = [threshold_row('chr15 (tests/data/rmats)', read_counts(RMATS))]
    for d in sorted((DATA / 'model').iterdir()):
        try:
            counts = read_counts(d)
        except Exception as e:
            summary.append(f'  {d.name}: not read ({type(e).__name__})')
            continue
        if counts:
            rows.append(threshold_row(d.name, counts))
        else:
            summary.append(f'  {d.name}: header only, skipped')
    if a.rmats:
        rows.append(threshold_row(f'whole genome ({a.rmats})', read_counts(a.rmats)))
    cols = list(dict.fromkeys(k for r in rows for k in r))
    (out / 'thresholds.tsv').write_text('\t'.join(cols) + '\n' + ''.join(
        '\t'.join(str(r.get(c, '')) for c in cols) + '\n' for r in rows))
    summary.append('thresholds.tsv (counts are min(summed IJC, summed SJC) unless named sjc):')
    for r in rows:
        summary.append(f"  {r['dataset']}: {r['junctions']} junctions; count 0: {r['min_count_0']}, count 1: "
                       f"{r['min_count_1']}; -r 1 keeps {r['kept_r1_now']} now, {r['kept_r1_as_documented']} as "
                       f"documented; SJC alone keeps {r['kept_sjc_only_r1_now']} ({r['sjc_ge2_but_min_lt2']} with "
                       f"SJC >= 2 but min < 2); model min_count {r['model_min_count']}" +
                       (f"; -m keeps {r['kept_model_now']} now (count >= {r['model_min_count'] + 1}), "
                        f"{r['kept_model_boundary']} at the model's boundary (count >= {r['model_min_count'] - 1}); "
                        f"{r['lost_to_offset']} lost to the offset" if 'kept_model_now' in r else ''))

    # chr15 fates and records under each rule (GTF canonicals, UniProt stubbed, as the golden test)
    chr15_mc = rows[0]['model_min_count']
    variants = [
        ('now: -r 1 (count > 1)', dict(read=1)),
        ('-r 1 as documented (count >= 1)', dict(read=0)),
        ('-r 1 on summed SJC (before 2022-05)', dict(read=1, sjc=True)),
        (f'now: -m (count > {chr15_mc})', dict(read=1, model=True, bound=chr15_mc)),
        (f'-m at the boundary (count >= {chr15_mc - 1})', dict(read=1, model=True, bound=chr15_mc - 2)),
    ]
    fates, files, per_junction = {}, {}, {}
    for name, v in variants:
        d = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='thr_'))
        f = golden_run(ctx, d, read=v['read'], model_flag=v.get('model', False), bound=v.get('bound', v['read']),
                       sjc=v.get('sjc', False))
        per_junction[name] = f
        fates[name] = collections.Counter(f.values())
        files[name] = records(d)
    labels = sorted({k for c in fates.values() for k in c})
    names = [n for n, _ in variants]
    lines = ['\t'.join(['fate or file'] + names)]
    lines += ['\t'.join([k] + [str(fates[n].get(k, 0)) for n in names]) for k in labels]
    lines += ['\t'.join([k] + [str(files[n].get(k, 0)) for n in names]) for k in sorted({k for c in files.values()
                                                                                          for k in c})]
    (out / 'threshold_fates.tsv').write_text('\n'.join(lines) + '\n')
    summary.append('threshold_fates.tsv (chr15, GTF canonicals; changed fates against the first rule of each pair):')
    for base, alt in [(names[0], names[1]), (names[0], names[2]), (names[3], names[4])]:
        changed = collections.Counter(f'{per_junction[base][j]} -> {per_junction[alt][j]}'
                                      for j in per_junction[base] if per_junction[base][j] != per_junction[alt][j])
        summary.append(f'  {alt}: {sum(changed.values())} junctions change: ' +
                       '; '.join(f'{k}: {n}' for k, n in changed.most_common()))
        summary.append('    records: ' + ', '.join(f'{k} {files[base].get(k, 0)} -> {files[alt].get(k, 0)}'
                                                   for k in sorted(set(files[base]) | set(files[alt]))
                                                   if files[base].get(k, 0) != files[alt].get(k, 0)))


def golden_run(ctx, out_dir, read, model_flag=False, bound=None, sjc=False, canonical=True, uniprot=None,
               genome=None, junctions=None, api=None):
    """ chr15 through main._translate_one with GTF canonicals and the UniProt API stubbed; returns fate per junction """
    args = argparse.Namespace(out=str(out_dir), read=read, model=model_flag, canonical=canonical, mask=None,
                              qvalue=[0, 1])
    fates = {}
    sum_sjc = Junction.sum_sjc
    if sjc:
        sum_sjc = property(lambda j: _total(j.sjc_s1, j.sjc_s2))
    with patched(params, 'use_gtf_only', True), patched(Junction, 'sum_sjc', sum_sjc), \
            patched(Sequence, 'get_canonical_aa_uniprot', api or (lambda self, reviewed='true': SeqRecord(Seq('')))):
        for jx in junctions or ctx['junctions']():
            fates[(jx.junction_type, jx.name)] = label(main._translate_one(
                jx, ctx['gtf'], genome or ctx['genome'], args, pred_bound=bound if bound is not None else read,
                uniprot=uniprot))
    return fates


# ---- real runs (PR 7) --------------------------------------------------------------------------------------------

def run(ctx, out_dir, clear_handlers=True, **overrides):
    """ chr15 through run_jcast with --uniprot, as the CLI does after parsing """
    args = argparse.Namespace(rmats_folder=str(RMATS), gtf_file=str(ctx['gtf_copy']), genome=str(GENOME),
                              out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1],
                              splice_type=set(TYPES), uniprot=str(UNIPROT), processes=1)
    for k, v in overrides.items():
        setattr(args, k, v)
    if clear_handlers:
        logging.getLogger('jcast').handlers = []
    try:
        main.run_jcast(args)
    finally:
        if clear_handlers:
            logging.getLogger('jcast').handlers = []


def log_runs(path):
    """ number of runs in a log file, and of junctions started in it """
    text = Path(path).read_text()
    return text.count('~INFO~Namespace('), text.count('>>>>>> Doing')


def part_runs(ctx, out, summary):
    lines = []
    w = ctx['work']

    # 1. the output directory must exist
    missing = w / 'no' / 'such' / 'dir'
    try:
        run(ctx, missing)
        lines.append('1. missing -o directory: ran')
    except Exception as e:
        lines.append(f'1. missing -o directory: {type(e).__name__}: {e}'.replace(str(w), '<work>'))

    # 2. a rerun into the same directory
    first, alone = w / 'rerun', w / 'r30_alone'
    first.mkdir()
    alone.mkdir()
    run(ctx, first, read=1)
    before = (records(first), digests(first))
    run(ctx, alone, read=30)
    run(ctx, first, read=30)
    after = (records(first), digests(first))
    lines.append('2. -r 1 into a directory, then -r 30 into the same directory:')
    lines.append(f'   records after -r 1:            {before[0]}')
    lines.append(f'   records of -r 30 on its own:   {records(alone)}')
    lines.append(f'   records after both:            {after[0]}')
    lines.append(f'   output files after both identical to the -r 1 run alone: {after[1] == before[1]}; '
                 f'to -r 30 alone: {after[1] == digests(alone)}')
    lines.append(f'   log: (runs, junctions started) {log_runs(first / "logfile.log")}')

    # 3. an exception in one junction, in one process and in three
    translate_one = main._translate_one

    def faulty(junction, *args, **kwargs):
        if (junction.junction_type, junction.name) == fault:
            raise RuntimeError('injected fault')
        return translate_one(junction, *args, **kwargs)

    fault = ctx['order'][FAULT_AT - 1]
    lines.append(f'3. an exception in junction {fault[0]} {fault[1]} (number {FAULT_AT} of {len(ctx["order"])} in run '
                 f'order):')
    for processes in (1, 3):
        d = w / f'fault_p{processes}'
        d.mkdir()
        with patched(main, '_translate_one', faulty):
            try:
                run(ctx, d, processes=processes)
                outcome = 'run completed'
            except Exception as e:
                outcome = f'run aborted: {type(e).__name__}: {e}'
        recs = records(d)
        lines.append(f'   -p {processes}: {outcome}; junctions logged {log_runs(d / "logfile.log")[1]}; '
                     f'records written {sum(recs.values())} of {sum(before[0].values())}; '
                     f'UniProt outputs written: {(d / "uniprot_canonicals.tsv").exists()}')

    # 4. logger handlers across runs in one process
    logging.getLogger('jcast').handlers = []
    a_dir, b_dir = w / 'log_a', w / 'log_b'
    a_dir.mkdir()
    b_dir.mkdir()
    run(ctx, a_dir, clear_handlers=False)
    h1 = [type(h).__name__ for h in logging.getLogger('jcast').handlers]
    a_lines = len((a_dir / 'logfile.log').read_text().splitlines())
    run(ctx, b_dir, clear_handlers=False)
    h2 = [type(h).__name__ for h in logging.getLogger('jcast').handlers]
    lines.append('4. two runs in one process (as the GUI or a notebook would):')
    lines.append(f'   handlers after run 1: {h1}; after run 2: {h2}')
    lines.append(f'   run 1 log lines after run 1: {a_lines}; after run 2: '
                 f'{len((a_dir / "logfile.log").read_text().splitlines())} (run 2 also writes into run 1\'s log)')
    lines.append(f'   run 2 log lines: {len((b_dir / "logfile.log").read_text().splitlines())}')
    for h in logging.getLogger('jcast').handlers:
        h.close()
    logging.getLogger('jcast').handlers = []

    # 5. the caches
    lines.append('5. caches:')
    cwd = w / 'cwd'
    cwd.mkdir()
    with chdir(cwd), patched(requests.Session, 'get', lambda self, url, **kw: response(200, '')):
        jx = next(iter(ctx['junctions']()))
        Sequence(jx).get_canonical_aa_uniprot()
    lines.append(f'   UniProt API cache written to the working directory: {(cwd / "cache" / "uniprot-cache.db").exists()}')
    ro = w / 'readonly_cwd'
    ro.mkdir()
    ro.chmod(0o555)
    try:
        with chdir(ro), patched(requests.Session, 'get', lambda self, url, **kw: response(200, '')):
            Sequence(jx).get_canonical_aa_uniprot()
        lines.append('   API path from a read-only working directory: no error')
    except Exception as e:
        lines.append(f'   API path from a read-only working directory: {type(e).__name__}: {e}'.replace(str(w), '<work>'))
    finally:
        ro.chmod(0o755)

    g = w / 'gtf_cache'
    g.mkdir()
    gtf_path = g / 'annotation.gtf'
    shutil.copy(GTF, gtf_path)
    full = ReadAnnotations(LOG, gtf_path)
    lines.append(f'   GTF cache written next to the GTF: {(g / "annotation.gtf.cached").exists()}; '
                 f'extra index column: {[c for c in full.annot.columns if c.startswith("Unnamed")]}')
    genes_full = full.annot.gene_id.nunique()
    text = GTF.read_text().splitlines(keepends=True)
    gtf_path.write_text(''.join(text[:len(text) // 2]))
    stale = ReadAnnotations(LOG, gtf_path)
    lines.append(f'   GTF replaced by its first half under the same name: genes read {stale.annot.gene_id.nunique()} '
                 f'(the replacement has {len({l.split("gene_id ")[1].split(";")[0] for l in text[:len(text) // 2] if "gene_id " in l})}, '
                 f'the old file {genes_full}): the stale cache is used')
    cached = g / 'annotation.gtf.cached'
    shutil.copy(GTF, gtf_path)
    cached.unlink()
    ReadAnnotations(LOG, gtf_path)
    data = cached.read_bytes()
    cached.write_bytes(data[:len(data) // 2])
    try:
        part = ReadAnnotations(LOG, gtf_path)
        lines.append(f'   cache half written (as a second run sees while the first writes it): read without error, '
                     f'{len(part.annot)} of {len(full.annot)} rows, {part.annot.gene_id.nunique()} of {genes_full} genes')
    except Exception as e:
        lines.append(f'   cache half written: {type(e).__name__}: {e}')
    r = w / 'gtf_readonly'
    r.mkdir()
    shutil.copy(GTF, r / 'annotation.gtf')
    r.chmod(0o555)
    try:
        ReadAnnotations(LOG, r / 'annotation.gtf')
        lines.append('   GTF in a read-only directory, no cache yet: read without error')
    except Exception as e:
        lines.append(f'   GTF in a read-only directory, no cache yet: {type(e).__name__}: {e}'.replace(str(w), '<work>'))
    finally:
        r.chmod(0o755)

    (out / 'runs.txt').write_text('\n'.join(lines) + '\n')
    summary.append('runs.txt:')
    summary.extend('  ' + line for line in lines)


# ---- canonicals under -c and the UniProt API (PR 8) ----------------------------------------------------------------

def response(status, text):
    r = requests.models.Response()
    r.status_code = status
    r._content = text.encode()
    r.encoding = 'utf-8'
    r.url = 'stub'
    return r


def api_stub(calls, cache, fail=None):
    """ a Session.get serving the local API cache by gene id without version, else an empty answer, counting calls """
    def get(self, url, **kwargs):
        caller = sys._getframe(2).f_code.co_name
        calls[caller] += 1
        if fail == 'connection':
            raise requests.exceptions.ConnectionError('stub: connection refused')
        if fail is not None:
            return response(fail, '')
        return response(200, cache.get(url.split('Ensembl:')[1].split('?')[0], ''))
    return get


def canonical_coverage(ctx, out_dir, uniprot):
    """ genes with a canonical and whether a record of it is written, from a chr15 -c run """
    found = {}
    get = Sequence.get_canonical_aa

    def spy(self, gtf, genome_index):
        get(self, gtf, genome_index)
        found[self.j.gene_id] = str(self.canonical_aa.seq)
    genes = collections.defaultdict(list)
    with patched(Sequence, 'get_canonical_aa', spy):
        fates = golden_run(ctx, out_dir, read=1, uniprot=uniprot)
    for (t, name), f in fates.items():
        genes[ctx['gene_of'][(t, name)]].append(f)
    written = {str(r.seq) for r in SeqIO.parse(out_dir / 'psq_canonical.fasta', 'fasta')}
    with_canonical = {g for g in genes if found.get(g)}
    missing = sorted(g for g in with_canonical if found[g] not in written)
    return genes, with_canonical, len(written), missing


def part_canonical(a, ctx, out, summary):
    rows, lines = [], []
    for name, uniprot in [('GTF only', None), ('--uniprot', ctx['uniprot'])]:
        d = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='canon_'))
        genes, has, n, missing = canonical_coverage(ctx, d, uniprot)
        summary.append(f'canonical.tsv: chr15 -c {name}: {len(genes)} genes with junctions, {len(has)} with a '
                       f'canonical, {n} canonical records, {len(missing)} genes with a canonical but no record: ' +
                       ', '.join(f'{g} ({dict(collections.Counter(genes[g]))})' for g in missing))
        rows += [(f'chr15 {name}', g, ctx['symbol'][g], ', '.join(f'{k} {v}' for k, v in
                                                                 collections.Counter(genes[g]).items()))
                 for g in missing]

    # the masked path: N in one exon of a gene's only junction
    by_gene = collections.defaultdict(list)
    for jx in ctx['junctions']():
        by_gene[jx.gene_id].append((jx.junction_type, jx.name))
    seq = ctx['genome'].genome['15']

    def fresh(key):
        return next(j for j in ctx['junctions']() if (j.junction_type, j.name) == key)
    for key in (j[0] for g, j in sorted(by_gene.items()) if len(j) == 1 and j[0][0] == 'SE'):
        jx = fresh(key)
        masked = argparse.Namespace(genome={'15': seq[:jx.alt1_es - 1] + 'N' * (jx.alt1_ee - jx.alt1_es + 1)
                                            + seq[jx.alt1_ee:]})
        d = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='mask_'))
        if list(golden_run(ctx, d, read=1, genome=masked, junctions=[fresh(key)]).values())[0].endswith('ambiguous bases'):
            break
    for name, uniprot in [('GTF only', None), ('--uniprot', ctx['uniprot'])]:
        d = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='mask_'))
        calls = collections.Counter()
        fates = golden_run(ctx, d, read=1, uniprot=uniprot, genome=masked, junctions=[fresh(key)],
                           api=lambda self, reviewed='true': calls.update([sys._getframe(1).f_code.co_name])
                           or SeqRecord(Seq('')))
        lines.append(f'masked path, {name}: {jx.gene_symbol} {jx.junction_type} {jx.name} (its gene\'s only junction) '
                     f'with N in the skipped exon: {list(fates.values())[0]}; canonical records written: '
                     f'{sum(1 for _ in SeqIO.parse(d / "psq_canonical.fasta", "fasta")) if (d / "psq_canonical.fasta").exists() else 0}; '
                     f'UniProt API lookups by caller: {dict(calls)}')

    # UniProt API calls by caller, chr15 through the API path with -c
    cache = {}
    if API_CACHE.exists():
        con = sqlite3.connect(API_CACHE)
        cache = {i.split('.')[0]: f for i, f in con.execute('SELECT id, seq FROM sequences')}
        con.close()
    for fail in [None, 429, 'connection']:
        d = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='api_'))
        calls = collections.Counter()
        cwd = Path(tempfile.mkdtemp(dir=ctx['work'], prefix='cwd_'))
        args = argparse.Namespace(out=str(d), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
        fates, error = {}, None
        with chdir(cwd), patched(params, 'use_gtf_only', False), \
                patched(requests.Session, 'get', api_stub(calls, cache, fail)):
            for jx in ctx['junctions']():
                try:
                    fates[(jx.junction_type, jx.name)] = label(main._translate_one(jx, ctx['gtf'], ctx['genome'],
                                                                                   args, pred_bound=1))
                except Exception as e:
                    error = f'{type(e).__name__}: {e}'
                    break
        name = {None: 'API answers (cache, else empty)', 429: 'API refuses every request (HTTP 429)',
                'connection': 'network down'}[fail]
        lines.append(f'API path, -c, {name}: junctions done {len(fates)} of {len(ctx["order"])}' +
                     (f', then the run aborts: {error}' if error else '') +
                     f'; requests by caller {dict(calls)}; records {records(d)}')
        if fail is None:
            ctx['api_ok'] = (fates, records(d))
        elif fail == 429:
            changed = sum(1 for k in fates if fates[k] != ctx['api_ok'][0].get(k))
            lines.append(f'   HTTP 429 against a working API: {changed} fates change; records '
                         + ', '.join(f'{k} {ctx["api_ok"][1].get(k, 0)} -> {records(d).get(k, 0)}'
                                     for k in sorted(set(records(d)) | set(ctx['api_ok'][1]))
                                     if records(d).get(k, 0) != ctx['api_ok'][1].get(k, 0)))

    if a.jcast_out:
        miss = whole_genome_missing(a)
        summary.append(f'canonical.tsv: whole genome ({a.jcast_out}): {miss[0]} genes with a canonical, '
                       f'{len(miss[1])} without a record: {dict(collections.Counter(s for _, s in miss[1]))}')
        rows += [('whole genome', g, '', s) for g, s in miss[1]]

    (out / 'canonical.tsv').write_text('run\tgene_id\tgene_symbol\tfates or canonical source\n' +
                                       ''.join('\t'.join(r) + '\n' for r in rows))
    (out / 'canonical.txt').write_text('\n'.join(lines) + '\n')
    summary.append('canonical.txt:')
    summary.extend('  ' + line for line in lines)


def whole_genome_missing(a):
    """ from a -c --uniprot run's uniprot_canonicals.tsv: genes with a canonical whose sequence is not written """
    import csv
    accession = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(a.uniprot, 'fasta')}
    recs = list(SeqIO.parse(Path(a.jcast_out) / 'psq_canonical.fasta', 'fasta'))
    seqs, ids = {str(r.seq) for r in recs}, {r.id for r in recs}
    n, missing = 0, []
    for r in csv.DictReader(open(Path(a.jcast_out) / 'uniprot_canonicals.tsv'), delimiter='\t'):
        if r['source'] == 'none':
            continue
        n += 1
        ok = (any(i.startswith(r['gene_id'] + '|') for i in ids) if r['source'] == 'gtf'
              else accession[r['accession']] in seqs)
        if not ok:
            missing.append((r['gene_id'], r['source']))
    return n, missing


def main_(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='runs_audit_'))
    os.chdir(work)
    logging.basicConfig(level=logging.WARNING)

    gtf_copy = work / GTF.name
    shutil.copy(GTF, gtf_copy)
    gtf = ReadAnnotations(LOG, gtf_copy)
    genome = ReadGenome(LOG, GENOME)
    rmats = RmatsResults(LOG, RMATS)

    def junctions():
        return [Junction(**rma.iloc[i].to_dict()) for rma in rmats.rmats_list for i in range(len(rma))]
    js = junctions()
    ctx = dict(work=work, gtf=gtf, gtf_copy=gtf_copy, genome=genome, junctions=junctions,
               order=[(j.junction_type, j.name) for j in js],
               gene_of={(j.junction_type, j.name): j.gene_id for j in js},
               symbol={j.gene_id: j.gene_symbol for j in js},
               uniprot=UniProtFasta(LOG, UNIPROT))

    commit = subprocess.run(['git', '-C', str(REPO), 'rev-parse', '--short', 'HEAD'], capture_output=True,
                            text=True).stdout.strip()
    dirty = subprocess.run(['git', '-C', str(REPO), 'status', '--porcelain', '--', 'jcast'], capture_output=True,
                           text=True).stdout.strip()
    summary = [f'code: {commit}{" + working tree" if dirty else ""}; python {sys.version.split()[0]}', '']
    if 'thresholds' in a.parts:
        part_thresholds(a, ctx, out, summary)
        summary.append('')
    if 'runs' in a.parts:
        part_runs(ctx, out, summary)
        summary.append('')
    if 'canonical' in a.parts:
        part_canonical(a, ctx, out, summary)
    (out / 'summary.txt').write_text('\n'.join(summary) + '\n')
    print('\n'.join(summary))
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--rmats', help='a whole-genome rMATS directory (read counts only)')
    p.add_argument('--jcast_out', help='output directory of a jcast -c --uniprot run')
    p.add_argument('--uniprot', help='the UniProt FASTA that run used')
    p.add_argument('--parts', nargs='+', default=['thresholds', 'runs', 'canonical'],
                   choices=['thresholds', 'runs', 'canonical'])
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'audit'))
    main_(p.parse_args())
