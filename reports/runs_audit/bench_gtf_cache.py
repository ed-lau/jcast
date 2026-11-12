# -*- coding: utf-8 -*-

"""GTF cache formats (reports/2026-11-12-thresholds-runs-canonical-audit.md §3.5): reading the current tab-separated
cache against a pickle of the same parsed table.

Reads an existing `<gtf>.cached` file as ReadAnnotations does, writes it as a pickle, reads that back and checks that
the table and its column types are identical. Writes gtf_cache.txt into --out.

Usage (from the repository root):
    python reports/runs_audit/bench_gtf_cache.py --cached <gtf>.cached --scratch <dir on a local disk>
"""

import argparse
import os
import time
from pathlib import Path

import pandas as pd


def main_(a):
    pickled = Path(a.scratch) / 'gtf_cache_bench.pkl'
    t = time.perf_counter()
    df = pd.read_table(a.cached, sep='\t', low_memory=False)
    tsv_read = time.perf_counter() - t
    t = time.perf_counter()
    df.to_pickle(pickled, protocol=5)
    pickle_write = time.perf_counter() - t
    t = time.perf_counter()
    back = pd.read_pickle(pickled)
    pickle_read = time.perf_counter() - t
    text = (f'{a.cached}: {len(df)} rows, {df.shape[1]} columns; pandas {pd.__version__}\n'
            f'tab-separated: {os.path.getsize(a.cached) / 1e9:.2f} GB, read {tsv_read:.1f} s\n'
            f'pickle: {os.path.getsize(pickled) / 1e9:.2f} GB, written {pickle_write:.1f} s, read {pickle_read:.1f} s\n'
            f'identical table and column types: {df.equals(back) and list(df.dtypes) == list(back.dtypes)}\n')
    pickled.unlink()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / 'gtf_cache.txt').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--cached', required=True)
    p.add_argument('--scratch', required=True)
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'audit'))
    main_(p.parse_args())
