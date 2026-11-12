# -*- coding: utf-8 -*-

"""Time spent on junctions skipped for low coverage (reports/2026-11-12-log-audit.md §3).

Every junction is trimmed, framed, sliced, translated and given a canonical before the read-count filter skips it.
Runs every chr15 junction through main._translate_one (GTF canonicals, UniProt API stubbed out) and times each one;
writes timing.txt into --out with the mean time per junction by fate, and the share of the junction loop spent on
junctions skipped for low coverage, here and at the whole-genome run's share of such junctions (--skipped_share).

Usage (from the repository root):
    python reports/log_audit/time_filter_first.py [--skipped_share 0.787]
"""

import argparse
import collections
import logging
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import main, params
from jcast.annots import ReadAnnotations, ReadGenome
from jcast.junctions import Junction, RmatsResults
from jcast.sequences import Sequence

DATA = REPO / 'tests' / 'data'


def main_(a):
    work = Path(tempfile.mkdtemp(prefix='time_filter_'))
    os.chdir(work)
    log = logging.getLogger('jcast-audit')
    gtf_copy = work / 'chr15.gtf'
    shutil.copy(DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf', gtf_copy)
    gtf = ReadAnnotations(log, gtf_copy)
    genome = ReadGenome(log, DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    rmats = RmatsResults(log, DATA / 'rmats')
    params.use_gtf_only = True
    Sequence.get_canonical_aa_uniprot = lambda self, reviewed='true': SeqRecord(Seq(''))
    out_dir = work / 'out'
    out_dir.mkdir()
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1])

    times = collections.defaultdict(list)
    for repeat in range(a.repeats):
        for rma in rmats.rmats_list:
            for i in range(len(rma)):
                jx = Junction(**rma.iloc[i].to_dict())
                t = time.perf_counter()
                fate = main._translate_one(jx, gtf, genome, args, pred_bound=1)
                if repeat == a.repeats - 1:
                    times['low coverage' if 'low coverage' in fate else 'other'].append(time.perf_counter() - t)
    low, other = times['low coverage'], times['other']
    mean_low, mean_other = sum(low) / len(low), sum(other) / len(other)
    share_here = sum(low) / (sum(low) + sum(other))
    s = a.skipped_share
    share_genome = s * mean_low / (s * mean_low + (1 - s) * mean_other)
    text = (f'chr15, GTF canonicals, timed on pass {a.repeats} of {a.repeats} (earlier passes warm the per-gene caches)\n'
            f'junctions skipped for low coverage: {len(low)}, mean {mean_low * 1000:.1f} ms\n'
            f'other junctions: {len(other)}, mean {mean_other * 1000:.1f} ms\n'
            f'share of the junction loop spent on low-coverage skips: {share_here:.0%} here; '
            f'{share_genome:.0%} at a {s:.1%} share of such junctions (the whole-genome run)\n')
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / 'timing.txt').write_text(text)
    print(text, end='')
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--skipped_share', type=float, default=471821 / 599644)
    p.add_argument('--repeats', type=int, default=2)
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'genome'))
    main_(p.parse_args())
