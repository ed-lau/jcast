# -*- coding: utf-8 -*-

"""Verification of PR 7 (reports/2026-11-12-thresholds-runs-canonical-audit.md §8) on chr15.

Each part runs in its own process with PYTHONPATH set to a checkout, so the code under test can be frozen in a git
worktree:

    identical       the chr15 test set through run_jcast with --uniprot, with -p 1 and -p 3, in the code before
                    (--before) and after (--after): every output file byte for byte, and the log from the first
                    junction on
    failure modes   the failure modes of the audit's runs.txt (§3.1-3.5), in the code after: a missing output
                    directory, a rerun into the same directory, an exception in one junction, two runs in one process,
                    and the caches

Writes verify_7.txt into --out.

Usage (from the repository root):
    python reports/runs_audit/verify_7.py --before <checkout of 82a309b> --after <checkout with PR 7>
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'
RMATS = DATA / 'rmats'
TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']
FAULT_AT = 500      # the junction made to raise, by its number in run order, as in the audit


# ---- workers: run inside one code version ----------------------------------------------------------------------------

def _args(out, gtf, **extra):
    a = dict(rmats_folder=str(RMATS), gtf_file=str(gtf), genome=str(GENOME), out=str(out), read=1, model=False,
             canonical=True, mask=None, qvalue=[0, 1], splice_type=set(TYPES), uniprot=str(UNIPROT), processes=1)
    a.update(extra)
    return argparse.Namespace(**a)


def _stub_api():
    from jcast.sequences import Sequence

    def no_api(self, reviewed='true'):
        raise AssertionError('UniProt API called with --uniprot')
    Sequence.get_canonical_aa_uniprot = no_api


def worker_run(out, processes):
    """ one chr15 run, as the code version does it; the old code needs the directory to exist """
    import logging
    from jcast import main
    _stub_api()
    gtf = Path(out).parent / f'gtf_{processes}' / GTF.name     # the old code writes <gtf>.cached next to the GTF
    gtf.parent.mkdir()
    shutil.copy(GTF, gtf)
    Path(out).mkdir()
    extra = {'cache_dir': str(Path(out).parent / 'cache')} if hasattr(main, 'existing_outputs') else {}
    main.run_jcast(_args(out, gtf, processes=processes, **extra))
    logging.getLogger('jcast').handlers = []


def worker_modes(work):
    """ the failure modes of the audit, in this code version; prints one JSON object """
    import logging
    import requests
    from Bio import SeqIO
    from jcast import cache, main
    from jcast.annots import ReadAnnotations
    from jcast.junctions import Junction, RmatsResults
    from jcast.sequences import Sequence
    real_api = Sequence.get_canonical_aa_uniprot
    _stub_api()
    work = Path(work)
    cache_dir = work / 'cache'
    log = logging.getLogger('jcast-verify')
    results = {}

    def recs(d):
        return sum(1 for p in Path(d).glob('psq_*.fasta') for _ in SeqIO.parse(p, 'fasta'))

    def run(out, **extra):
        main.run_jcast(_args(out, GTF, cache_dir=str(cache_dir), **extra))

    # 1. a missing output directory
    run(work / 'no' / 'such' / 'dir')
    results['1. missing -o directory'] = f'created; {recs(work / "no" / "such" / "dir")} records written'

    # 2. a rerun into the same directory
    run(work / 'r1')
    run(work / 'r30_alone', read=30)
    try:
        run(work / 'r1', read=30)
        refused = 'not refused'
    except FileExistsError as e:
        refused = f'refused: {str(e)[:80]}...'
    run(work / 'r1', read=30, overwrite=True)
    same = all((work / 'r1' / p.name).read_bytes() == p.read_bytes()
               for p in (work / 'r30_alone').glob('psq_*.fasta')) and \
        sorted(p.name for p in (work / 'r1').glob('psq_*.fasta')) == \
        sorted(p.name for p in (work / 'r30_alone').glob('psq_*.fasta'))
    results['2. -r 30 into a directory holding a -r 1 run'] = (
        f'{refused}; with --overwrite: psq files identical to -r 30 alone: {same}; '
        f'runs in the log: {(work / "r1" / "logfile.log").read_text().count("~INFO~Namespace(")}')

    # 3. an exception in one junction
    run(work / 'reference')
    rmats = RmatsResults(log, RMATS)
    order = [(t, str(rma.iloc[i]['id'])) for t, rma in zip(TYPES, rmats.rmats_list) for i in range(len(rma))]
    fault = order[FAULT_AT - 1]
    translate_one = main._translate_one

    def faulty(junction, *a, **k):
        if (junction.junction_type, junction.name) == fault:
            raise RuntimeError('injected fault')
        return translate_one(junction, *a, **k)
    main._translate_one = faulty
    for processes in (1, 3):
        d = work / f'fault_p{processes}'
        run(d, processes=processes)
        text = (d / 'logfile.log').read_text()
        results[f'3. exception in junction {fault[0]} {fault[1]} (number {FAULT_AT}), -p {processes}'] = (
            f'run completed; junctions logged {text.count(">>>>>> Doing")} of {len(order)}; exception fates '
            f'{text.count("FAILURE (exception). RuntimeError: injected fault.")}; records {recs(d)} of '
            f'{recs(work / "reference")}; UniProt outputs written: {(d / "uniprot_canonicals.tsv").exists()}')
    main._translate_one = translate_one
    p1, p3 = work / 'fault_p1', work / 'fault_p3'
    results['3. -p 1 and -p 3 psq files identical'] = all(
        (p1 / p.name).read_bytes() == p.read_bytes() for p in p3.glob('psq_*.fasta'))

    # 4. two runs in one process
    run(work / 'log_a')
    a_log = (work / 'log_a' / 'logfile.log').read_text()
    run(work / 'log_b')
    results['4. two runs in one process'] = (
        f'run 1 log unchanged by run 2: {(work / "log_a" / "logfile.log").read_text() == a_log}; handlers after the '
        f'runs: {len(logging.getLogger("jcast").handlers)}; log lines: {len(a_log.splitlines())} and '
        f'{len((work / "log_b" / "logfile.log").read_text().splitlines())}')

    # 5. the caches
    cache.set_directory(cache_dir)
    cwd = work / 'cwd'
    cwd.mkdir()
    os.chdir(cwd)
    ok = requests.models.Response()
    ok.status_code, ok._content, ok.encoding = 200, b'', 'utf-8'
    requests.Session.get = lambda self, url, **kw: ok
    Sequence.get_canonical_aa_uniprot = real_api
    jx = Junction(**rmats.rmats_se.iloc[0].to_dict())
    Sequence(jx).get_canonical_aa_uniprot()
    results['5. UniProt API cache'] = (f'in the cache directory: {(cache_dir / "uniprot-cache.db").exists()}; '
                                       f'under the working directory: {(cwd / "cache").exists()}')
    ro = work / 'readonly_cwd'
    ro.mkdir()
    ro.chmod(0o555)
    os.chdir(ro)
    try:
        Sequence(jx).get_canonical_aa_uniprot()
        results['5. API path from a read-only working directory'] = 'no error'
    except Exception as e:
        results['5. API path from a read-only working directory'] = f'{type(e).__name__}: {e}'
    finally:
        os.chdir(work)
        ro.chmod(0o755)

    g = work / 'gtf_dir'
    g.mkdir()
    path = g / 'annotation.gtf'
    shutil.copy(GTF, path)
    full = ReadAnnotations(log, path)
    genes_full = full.annot.gene_id.nunique()
    results['5. GTF cache'] = (f'files next to the GTF: {sorted(p.name for p in g.iterdir())}; caches of it: '
                               f'{len(list((cache_dir / "gtf").glob("annotation.gtf.*.pkl")))}; index column: '
                               f'{[c for c in full.annot.columns if c.startswith("Unnamed")]}')
    text = GTF.read_text().splitlines(keepends=True)
    path.write_text(''.join(text[:len(text) // 2]))
    half = ReadAnnotations(log, path).annot.gene_id.nunique()
    results['5. GTF replaced by its first half under the same name'] = (
        f'genes read {half} (the old file has {genes_full}); caches of it: '
        f'{len(list((cache_dir / "gtf").glob("annotation.gtf.*.pkl")))}')
    shutil.copy(GTF, path)
    ReadAnnotations(log, path)
    cached = next((cache_dir / 'gtf').glob('annotation.gtf.*.pkl'))
    data = cached.read_bytes()
    Path(f'{cached}.999.tmp').write_bytes(data[:len(data) // 2])     # another run's partial cache
    part = ReadAnnotations(log, path).annot
    results['5. another run\'s partial cache present'] = f'not read: {len(part)} of {len(full.annot)} rows'
    cached.write_bytes(data[:len(data) // 2])                         # a damaged cache
    part = ReadAnnotations(log, path).annot
    results['5. damaged cache'] = f'read again: {len(part)} of {len(full.annot)} rows'
    r = work / 'gtf_readonly'
    r.mkdir()
    shutil.copy(GTF, r / 'annotation.gtf')
    r.chmod(0o555)
    try:
        n = ReadAnnotations(log, r / 'annotation.gtf').annot.gene_id.nunique()
        results['5. GTF in a read-only directory'] = f'read without error, {n} genes'
    except Exception as e:
        results['5. GTF in a read-only directory'] = f'{type(e).__name__}: {e}'
    finally:
        r.chmod(0o755)
    print(json.dumps(results))


# ---- driver ------------------------------------------------------------------------------------------------------

def _in(repo, *args, cwd):
    out = subprocess.run([sys.executable, __file__, *args], cwd=cwd, capture_output=True, text=True,
                         env=dict(os.environ, PYTHONPATH=str(repo), MPLBACKEND='Agg'))
    if out.returncode:
        raise RuntimeError(out.stderr[-3000:])
    return out.stdout


def _messages(path):
    lines = Path(path).read_text().split('\n')
    messages = [line.split('~', 1)[1] if line[:4].isdigit() else line for line in lines]
    return messages[next(i for i, m in enumerate(messages) if '>>>>>> Doing' in m):]


def _commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + working tree' if d.stdout.strip() else '')


def main_(a):
    work = Path(tempfile.mkdtemp(prefix='verify_7_'))
    lines = [f'before: {_commit(a.before)}; after: {_commit(a.after)}', '', 'identical (chr15, --uniprot, -c):']
    for processes in (1, 3):
        outs = {}
        for label, repo in [('before', a.before), ('after', a.after)]:
            d = work / f'{label}_p{processes}'
            d.mkdir()
            _in(repo, '--worker_run', str(d / 'out'), '--processes', str(processes), cwd=d)
            outs[label] = d / 'out'
        files = {k: {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                     for p in sorted(v.iterdir()) if p.name != 'logfile.log'} for k, v in outs.items()}
        same_log = _messages(outs['before'] / 'logfile.log') == _messages(outs['after'] / 'logfile.log')
        lines.append(f'  -p {processes}: {len(files["after"])} output files, identical: '
                     f'{files["before"] == files["after"]}; log from the first junction identical: {same_log}')
    modes = work / 'modes'
    modes.mkdir()
    lines += ['', 'failure modes (after):']
    lines += [f'  {k}: {v}' for k, v in json.loads(_in(a.after, '--worker_modes', str(modes), cwd=modes)
                                                    .strip().splitlines()[-1]).items()]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'verify_7.txt').write_text('\n'.join(l.replace(str(work), '<work>') for l in lines) + '\n')
    print((out / 'verify_7.txt').read_text())
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--before')
    p.add_argument('--after')
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'verify_7'))
    p.add_argument('--worker_run')
    p.add_argument('--processes', type=int, default=1)
    p.add_argument('--worker_modes')
    a = p.parse_args()
    if a.worker_run:
        worker_run(a.worker_run, a.processes)
    elif a.worker_modes:
        worker_modes(a.worker_modes)
    else:
        main_(a)
