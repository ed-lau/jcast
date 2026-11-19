# -*- coding: utf-8 -*-

"""Verification of PR 27 (reports/2026-11-15-read-count-model-audit.md §18): the Model page's two charts of what a
cutoff keeps (§16.2): the junctions a run would translate at each cutoff, annotated and novel, and by read count the
events whose alternative region is a multiple of 3 nt.

    chr15       the bundled rMATS folder, default and -m, run with the code before and after (verify_24.py's runs):
                output files byte-identical, logs identical but for timestamps, the output directory and the time taken
    events      (--folder) each rMATS folder read as the Model page reads it (model.preview, the code after), against
                the audit's own reading of the rMATS files (cutoff_recovery.frame_kept): for every event, whether its
                alternative region is a multiple of 3, and whether rMATS calls it novel
    bins        (--bins) the frame chart's numbers for the junctions a -r 1 run translated (T1, T1_alt_terminus, T2,
                T4), against the audit's bins.tsv: §4.4's table, which the chart shows for every event
    runs        (--run) the junctions translated at a run's cutoff, as the chart counts them, against the run's
                junctions.tsv: the junctions neither skipped_low_count nor skipped_fdr
    time        model.preview on each --folder with the code before and after, each twice, alternately
    screenshot  (--screenshot) the Model page's two charts for a folder, at its -m cutoff (the code after)

Writes <out>/verify_27.txt, and with --screenshot <out>/kept.png and <out>/frame.png.

Usage (from the repository root):
    python reports/model_audit/verify_27.py --before <checkout> --after <checkout> [--folder <label>=<rMATS> ...]
        [--bins <label>=<rMATS>,<-r 1 output dir>,<bins.tsv> ...] [--run <label>=<rMATS>,<output dir> ...]
        [--screenshot <rMATS>] [--work <dir>] [--out <dir>]
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cutoff_recovery as cr  # noqa: E402
import verify_24 as v24  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
KEY = ['junction_type', 'junction_id']

TIME = """
import sys, time
from jcast import model
t = time.time()
model.preview(sys.argv[1])
print(round(time.time() - t, 1))
"""

SCREENSHOT = """
import sys
from PySide6.QtCore import QEventLoop, QSettings, QTimer
from PySide6.QtWidgets import QApplication
app = QApplication([])
from jcast.gui import main_window
main_window.user_settings = lambda: QSettings('settings.ini', QSettings.IniFormat)


def spin(ms):
    loop = QEventLoop()
    QTimer.singleShot(ms, loop.quit)
    loop.exec()


w = main_window.MainWindow()
w.resize(980, 800)
w.show()
w.setup.rmats_folder.setText(sys.argv[1])
w.setup.refresh()
w.tabs.setCurrentWidget(w.model)
while w.model._fit is None and 'could not' not in w.model.status.text():
    spin(100)
w.model.use_model.setChecked(True)
spin(500)
for i, name in [(1, 'kept'), (2, 'frame')]:
    w.model.charts.setCurrentIndex(i)
    spin(300)
    w.grab().save(f'{sys.argv[2]}/{name}.png')
print(w.model.status.text())
"""


def events(folder) -> list:
    """ the Model page's reading of each event against the audit's """
    from jcast import model
    counts = model.preview(folder)[0]
    audit = cr.frame_kept(folder).assign(junction_id=lambda d: d.junction_id.astype(int))
    both = counts.merge(audit, on=KEY, how='left', suffixes=('', '_audit'), indicator=True)
    found = (both['_merge'] == 'both').all()
    frame = (both.frame_kept == both.frame_kept_audit).all()
    novel = (both.novel.fillna(False).to_numpy(dtype=bool) == both.novel_audit.to_numpy(dtype=bool)).all()
    called = 'no novel-event files' if counts.novel.isna().all() else \
        f'{int(counts.novel.sum()):,} novel ({counts.novel.mean():.1%})'
    return [f'{len(counts):,} events; each in the audit\'s reading: {found}; the frame as the audit\'s: {frame}; '
            f'rMATS\'s call as the audit\'s: {novel} ({called}); alternative region a multiple of 3: '
            f'{counts.frame_kept.mean():.1%}']


def bins(folder, run, table) -> list:
    """ the frame chart for the junctions a -r 1 run translated, against the audit's bins.tsv ('all') """
    from jcast import model
    counts = model.preview(folder)[0]
    fates = pd.read_csv(Path(run) / 'junctions.tsv', sep='\t', usecols=KEY + ['fate'])
    translated = counts.merge(fates, on=KEY)
    translated = translated[translated.fate.isin(cr.TRANSLATED)]
    frame = model.frame_by_count(translated)
    audit = pd.read_csv(table, sep='\t')
    audit = audit[audit.splice_type == 'all'].set_index('bin')
    lines, differ = [], 0
    for group in ('annotated', 'novel'):
        mine = frame[f'{group} frame kept'].round(4).rename(index=int)
        theirs = audit[f'frame kept, {group} events'].reindex(mine.index)
        differ += int((mine.fillna(-1) != theirs.fillna(-1)).sum())
        lines.append(f'  {group}: ' + ', '.join(f'{b} {100 * v:.1f}%' for b, v in mine.dropna().items()))
    events_equal = (frame['annotated events'] + frame['novel events']).rename(index=int).equals(
        audit.junctions.reindex(frame.index.astype(int)).fillna(0).astype(int))
    return [f'{len(translated):,} junctions translated; per bin as bins.tsv: junctions {events_equal}, the shares '
            f'(to 4 decimals) {"all equal" if not differ else f"{differ} differ"}'] + lines


def runs(folder, out) -> list:
    """ the junctions translated at the run's cutoff, as the chart counts them, against its junctions.tsv """
    from jcast import model
    log = (Path(out) / 'logfile.log').read_text()
    cutoff = int(re.search(r'with a read count of at least (\d+)', log).group(1))
    q = tuple(float(v) for v in re.search(r'qvalue=\[([\d.e-]+), ([\d.e-]+)\]', log).groups())
    types = re.findall(r"'(\w+)'", re.search(r'splice_type=\{([^}]*)\}', log).group(1))
    counts = model.preview(folder)[0]
    at = model.kept_at_cutoffs(counts, q, types, cutoffs=[cutoff]).iloc[0]
    fates = pd.read_csv(Path(out) / 'junctions.tsv', sep='\t', usecols=['fate']).fate
    ran = int((~fates.isin(['skipped_low_count', 'skipped_fdr'])).sum())
    return [f'cutoff {cutoff}, -q {q[0]:g} {q[1]:g}: the chart {int(at.sum()):,} (' +
            ', '.join(f'{n:,} {g}' for g, n in at.items()) + f'); the run {ran:,}; equal: {int(at.sum()) == ran}']


def timed(repo, folder, work) -> float:
    done = subprocess.run([sys.executable, '-c', TIME, str(folder)], capture_output=True, text=True, cwd=work,
                          env=v24.env(repo), check=True)
    return float(done.stdout.strip().splitlines()[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--folder', action='append', default=[], help='<label>=<rMATS folder>')
    parser.add_argument('--bins', action='append', default=[], help='<label>=<rMATS>,<-r 1 output dir>,<bins.tsv>')
    parser.add_argument('--run', action='append', default=[], help='<label>=<rMATS>,<output dir>')
    parser.add_argument('--screenshot', help='an rMATS folder')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_27')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_27')
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.after).resolve()))     # the code after, for model.preview here

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    args.out.mkdir(parents=True, exist_ok=True)
    cache_dir = args.work / 'cache'
    v24.run(args.after, v24.DATA / 'rmats', ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)

    lines = [f'{side}: {v24.imported(getattr(args, side), args.work).replace(str(args.work), "<work>")}'
             for side in ('before', 'after')]
    lines += ['', '== chr15 runs: before and after']
    for label, extra in v24.RUNS.items():
        (rb, b), (ra, a) = (v24.run(getattr(args, side), v24.DATA / 'rmats', extra,
                                    args.work / re.sub(r'[^A-Za-z0-9]+', '_', f'{side} {label}'), cache_dir)
                            for side in ('before', 'after'))
        fb, fa = v24.files(b), v24.files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lb, la = v24.log_lines(b), v24.log_lines(a)
        lines.append(f'{label}: exit {rb} / {ra}; {len(fa)} output files, identical: '
                     f'{"yes" if not differ else "no, " + ", ".join(differ)}; logs identical: {lb == la}')

    folders = [f.split('=', 1) for f in args.folder]
    lines += ['', '== events: the Model page\'s reading against the audit\'s (cutoff_recovery.frame_kept)']
    for label, folder in folders:
        lines += [f'{label}: ' + events(folder)[0]]

    lines += ['', '== the frame chart for the junctions a -r 1 run translated, against the audit\'s bins.tsv']
    for spec in args.bins:
        label, paths = spec.split('=', 1)
        found = bins(*paths.split(','))
        lines += [f'{label}: {found[0]}'] + found[1:]

    lines += ['', '== the junctions translated at a run\'s cutoff: the chart against the run']
    for spec in args.run:
        label, paths = spec.split('=', 1)
        lines += [f'{label}: ' + runs(*paths.split(','))[0]]

    lines += ['', '== model.preview, seconds: before, after, before, after']
    for label, folder in folders:
        times = [timed(getattr(args, side), folder, args.work) for _ in range(2) for side in ('before', 'after')]
        lines.append(f'{label}: ' + ', '.join(f'{t:g}' for t in times))

    if args.screenshot:
        env = dict(v24.env(args.after), QT_QPA_PLATFORM='offscreen')
        done = subprocess.run([sys.executable, '-c', SCREENSHOT, args.screenshot, str(args.out.resolve())],
                              capture_output=True, text=True, cwd=args.work, env=env, check=True)
        lines += ['', f'== screenshots (kept.png, frame.png): {args.screenshot}', done.stdout.strip()]

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    (args.out / 'verify_27.txt').write_text(text)


if __name__ == '__main__':
    main()
