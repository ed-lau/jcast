# -*- coding: utf-8 -*-

"""Verification of PR 8 (reports/2026-11-12-thresholds-runs-canonical-audit.md §10).

Each part runs in its own process with PYTHONPATH set to a checkout, so the code under test can be frozen in a git
worktree:

    files       the chr15 test set through the CLI with -c --uniprot (-r 1 and -r 2; -p 1 and -p 3), before and after:
                every output file but psq_canonical.fasta byte for byte; psq_canonical.fasta as sets of records
    coverage    for each after run, the genes with a canonical (uniprot_canonicals.tsv) that have no record in
                psq_canonical.fasta
    api         the UniProt API path on chr15 (-c, no --uniprot), with the API replaced by a stub serving the local
                API cache (cache/uniprot-cache.db, keys without the gene version), before and after: the requests
                by caller, and what a refusal (HTTP 429) or a network failure does to the run's fates and records
    genome      with --genome_before/--genome_after (output directories of whole-genome runs, e.g. from verify_6.py)
                and --uniprot_full (the FASTA they used): the same file comparison and coverage, and the genes whose
                canonical source changed

Writes verify_8.txt into --out.

Usage (from the repository root):
    python reports/runs_audit/verify_8.py --before <checkout> --after <checkout> \\
        [--genome_before <dir> --genome_after <dir> --uniprot_full <fasta>]
"""

import argparse
import collections
import csv
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'
RMATS = DATA / 'rmats'
API_CACHE = REPO / 'cache' / 'uniprot-cache.db'


# ---- worker: the API path, inside one code version ------------------------------------------------------------------

def worker_api(work, fail):
    """ chr15 through _translate_safely (or _translate_one) junction by junction, as a run does; prints JSON """
    import argparse as ap
    import logging
    import requests
    from jcast import cache, main, params
    from jcast.annots import ReadAnnotations, ReadGenome
    from jcast.junctions import Junction, RmatsResults

    work = Path(work)
    log = logging.getLogger('jcast-verify')
    if hasattr(cache, 'set_directory'):
        cache.set_directory(work / 'cache')
    os.chdir(work)
    gtf, genome, rmats = ReadAnnotations(log, GTF), ReadGenome(log, GENOME), RmatsResults(log, RMATS)
    con = sqlite3.connect(API_CACHE)
    answers = {i.split('.')[0]: f for i, f in con.execute('SELECT id, seq FROM sequences')}
    con.close()
    calls = collections.Counter()

    def get(self, url, **kwargs):
        calls[sys._getframe(2).f_code.co_name] += 1
        if fail == 'connection':
            raise requests.exceptions.ConnectionError('stub: connection refused')
        r = requests.models.Response()
        r.status_code, r.encoding, r.url = (int(fail) if fail else 200), 'utf-8', 'stub'
        r._content = b'' if fail else answers.get(url.split('Ensembl:')[1].split('?')[0], '').encode()
        return r
    requests.Session.get = get
    params.use_gtf_only = False

    out = work / 'out'
    out.mkdir()
    args = ap.Namespace(out=str(out), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])
    fates, error = collections.Counter(), None
    from jcast import helpers as h
    for rma in rmats.rmats_list:
        for i in range(len(rma)):
            jx = Junction(**rma.iloc[i].to_dict())
            try:
                if hasattr(main, '_translate_safely'):
                    fate, records = main._translate_safely(jx, gtf, genome, args, pred_bound=1, uniprot=None)
                    for output, suffix, record in records:
                        h.write_seqrecord_to_fasta(record, output, suffix)
                else:
                    fate = main._translate_one(jx, gtf, genome, args, pred_bound=1)
            except Exception as e:
                error = f'{type(e).__name__}: {e}'
                break
            fates[fate.split('.')[0] if not fate.startswith('SKIPPED') else 'SKIPPED'] += 1
    records = {p.name: sum(1 for _ in SeqIO.parse(p, 'fasta')) for p in sorted(out.glob('psq_*.fasta'))}
    print(json.dumps(dict(fates=fates, calls=calls, records=records, error=error)))


# ---- driver ------------------------------------------------------------------------------------------------------

def _in(repo, args, cwd):
    out = subprocess.run([sys.executable] + args, cwd=cwd, capture_output=True, text=True,
                         env=dict(os.environ, PYTHONPATH=str(repo), MPLBACKEND='Agg'))
    if out.returncode:
        raise RuntimeError(out.stderr[-3000:])
    return out.stdout


def _commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + working tree' if d.stdout.strip() else '')


def _records(path):
    return [(r.description, str(r.seq)) for r in SeqIO.parse(path, 'fasta')] if Path(path).exists() else []


def compare(before, after, uniprot):
    """ lines comparing two output directories of -c --uniprot runs """
    files = sorted(p.name for p in Path(after).iterdir() if p.name not in ('logfile.log', 'model.pdf'))
    other = [f for f in files if f != 'psq_canonical.fasta']
    differ = [f for f in other if not (Path(before) / f).exists() or
              (Path(before) / f).read_bytes() != (Path(after) / f).read_bytes()]
    b, a = _records(Path(before) / 'psq_canonical.fasta'), _records(Path(after) / 'psq_canonical.fasta')
    lines = [f'    other output files: {len(other)}; differing: {differ or "none"}',
             f'    psq_canonical.fasta: {len(b)} -> {len(a)} records; every record before is still there: '
             f'{set(b) <= set(a)}; added {len(set(a) - set(b))}; removed {len(set(b) - set(a))}']
    lines.append('    ' + coverage(after, uniprot))
    return lines, differ


def coverage(out, uniprot):
    """ genes with a canonical (uniprot_canonicals.tsv) whose canonical sequence is not in psq_canonical.fasta """
    import gzip
    opener = gzip.open if str(uniprot).endswith('.gz') else open
    with opener(uniprot, 'rt') as fh:
        accession = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(fh, 'fasta')}
    recs = _records(Path(out) / 'psq_canonical.fasta')
    seqs = {s for _, s in recs}
    ids = {d.split()[0] for d, _ in recs}
    n, missing = 0, []
    for r in csv.DictReader(open(Path(out) / 'uniprot_canonicals.tsv'), delimiter='\t'):
        if r['source'] == 'none':
            continue
        n += 1
        ok = (any(i.startswith(r['gene_id'] + '|') for i in ids) if r['source'] == 'gtf'
              else accession[r['accession']] in seqs)
        if not ok:
            missing.append(r['gene_symbol'])
    return f'-c coverage: {n} genes with a canonical, {len(missing)} without a record {missing[:10]}'


def source_changes(before, after):
    rows = {}
    for label, d in (('before', before), ('after', after)):
        for r in csv.DictReader(open(Path(d) / 'uniprot_canonicals.tsv'), delimiter='\t'):
            rows.setdefault(r['gene_id'], {})[label] = (r['source'], r['accession'], r['gene_symbol'])
    changed = {g: v for g, v in rows.items() if v.get('before') != v.get('after')}
    return [f'    canonical source changed for {len(changed)} genes:'] + [
        f'      {v["after"][2]} {g}: {v["before"][0]} {v["before"][1]} -> {v["after"][0]} {v["after"][1]}'
        for g, v in sorted(changed.items())]


def main_(a):
    work = Path(tempfile.mkdtemp(prefix='verify_8_'))
    lines = [f'before: {_commit(a.before)}; after: {_commit(a.after)}', '', 'files (chr15, CLI, -c --uniprot):']
    for extra in (['-r', '1'], ['-r', '2'], ['-r', '1', '-p', '3']):
        outs = {}
        for label, repo in (('before', a.before), ('after', a.after)):
            d = work / f'{label}_{"_".join(extra)}'
            _in(repo, ['-m', 'jcast', str(RMATS), str(GTF), str(GENOME), '-c', '--uniprot', str(UNIPROT), '-o', str(d),
                       '--cache_dir', str(work / 'cache')] + extra, cwd=work)
            outs[label] = d
        lines.append(f'  {" ".join(extra)}:')
        lines += compare(outs['before'], outs['after'], UNIPROT)[0]
        lines.append('    before: ' + coverage(outs['before'], UNIPROT))

    lines += ['', 'api (chr15, -c, the UniProt API path with a stub serving the local API cache):']
    for fail in ('', '429', 'connection'):
        name = {'': 'answering', '429': 'refusing every request (HTTP 429)', 'connection': 'network down'}[fail]
        lines.append(f'  {name}:')
        for label, repo in (('before', a.before), ('after', a.after)):
            d = work / f'api_{label}_{fail or "ok"}'
            d.mkdir()
            r = json.loads(_in(repo, [__file__, '--worker_api', str(d), '--fail', fail], cwd=d).strip().splitlines()[-1])
            lines.append(f'    {label}: fates {r["fates"]}; requests by caller {r["calls"]}; records {r["records"]}' +
                         (f'; the run aborts: {r["error"]}' if r['error'] else ''))

    if a.genome_before:
        lines += ['', f'genome (-m -c --uniprot -p 16; {a.genome_before} -> {a.genome_after}):']
        lines += compare(a.genome_before, a.genome_after, a.uniprot_full)[0]
        lines.append('    before: ' + coverage(a.genome_before, a.uniprot_full))
        lines += source_changes(a.genome_before, a.genome_after)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'verify_8.txt').write_text('\n'.join(l.replace(str(work), '<work>') for l in lines) + '\n')
    print((out / 'verify_8.txt').read_text())
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--before')
    p.add_argument('--after')
    p.add_argument('--genome_before')
    p.add_argument('--genome_after')
    p.add_argument('--uniprot_full')
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'verify_8'))
    p.add_argument('--worker_api')
    p.add_argument('--fail', default='')
    a = p.parse_args()
    if a.worker_api:
        worker_api(a.worker_api, a.fail)
    else:
        main_(a)
