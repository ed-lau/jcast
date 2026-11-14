# -*- coding: utf-8 -*-

"""Verification of PR 21 (reports/2026-11-13-gui-cli-scope.md §13): the GUI after its first run on a Mac.

    pipeline      run --progress names the last junction done, and reports at least every PROGRESS_SECONDS. chr15
                  runs (<inputs> -c --uniprot) with the code before and after, with -p 1 and -p 2, without and with
                  --progress: output files byte-identical, logs identical but for timestamps, the output directory and
                  the time taken. The progress lines: how many, the longest time between two about junctions, and
                  (after) whether each names the junction at its place in junctions.tsv.
    packages      the wheel and the sdist built from the code after hold jcast/gui/icon.png (--build-python: a Python
                  with the build package).
    screenshots   with the code after: the Setup page with the help line, and the Run page during a run, saved in
                  <out> (paths relative to the scratch directory, where data links to the test data).

Each process runs from a scratch directory, with PYTHONPATH set to its checkout.

Writes <out>/verify_21.txt and <out>/<page>.png.

Usage (from the repository root):
    python reports/gui_scope/verify_21.py --before <checkout> --after <checkout> [--build-python <python>]
        [--work <dir>] [--out <dir>]
"""

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
INPUTS = [DATA / 'rmats', DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', '-c', '--uniprot',
          DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz']
CASES = {'-p 1': [], '-p 2': ['-p', '2']}
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
SET = re.compile(r"splice_type=\{([^}]*)\}")


def env(repo):
    return dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg', COLUMNS='1000',
                QT_QPA_PLATFORM='offscreen')


def run(repo, extra, out, cache_dir, progress):
    """ python -P -m jcast run, as the GUI runs it; with --progress, each line of stdout and when it came """
    out.mkdir(parents=True)
    cmd = [sys.executable, '-P', '-m', 'jcast', 'run'] + INPUTS + ['-o', out, '--cache_dir', cache_dir] + extra
    cmd += ['--progress'] if progress else []
    process = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                               cwd=out, env=env(repo))
    lines = [(time.monotonic(), line) for line in process.stdout]
    return process.wait(), [(t, json.loads(line)) for t, line in lines if line.startswith('{')]


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(out).iterdir())
            if p.is_file() and p.name not in ('logfile.log', 'model.pdf')}


def log_lines(out):
    text = (Path(out) / 'logfile.log').read_text().replace(str(out), '<out>')
    text = SET.sub(lambda m: 'splice_type={' + ', '.join(sorted(m.group(1).split(', '))) + '}', text)
    return [STAMP.sub('', line) for line in text.splitlines() if not STAMP.sub('', line).startswith('INFO~Time: ')]


def compare(a, b) -> str:
    fa, fb = files(a), files(b)
    differ = sorted(n for n in set(fa) | set(fb) if fa.get(n) != fb.get(n))
    return (f'{len(fa)} output files, identical: {"yes" if not differ else "NO, " + ", ".join(differ)}; log identical '
            f'but for timestamps, the output directory and the time taken: '
            f'{"yes" if log_lines(a) == log_lines(b) else "NO"}')


def progress_summary(steps, out) -> str:
    junctions = [(t, s) for t, s in steps if s['stage'] == 'junctions']
    gaps = [b[0] - a[0] for a, b in zip(junctions, junctions[1:])]
    named = [s for _, s in junctions if 'gene' in s]
    text = f'{len(steps)} lines, {len(junctions)} about junctions, at most {max(gaps):.2f} s apart'
    if not named:
        return text + '; none names a junction'
    with open(Path(out) / 'junctions.tsv') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    by_type = {t: [r for r in rows if r['junction_type'] == t] for t in {r['junction_type'] for r in rows}}
    right = sum(1 for s in named if (s['junction'], s['gene'], s['gene_id']) ==
                (lambda r: (r['junction_id'], r['gene_symbol'], r['gene_id']))(by_type[s['splice_type']][s['done'] - 1]))
    return text + f'; {len(named)} name the last junction done, {right} of them the junction at that place in ' \
                  f'junctions.tsv'


# ---- the screenshots, in the checkout's process ----

def gui_probe(spec: dict) -> None:
    from PySide6.QtCore import QEventLoop, QPointF, QTimer
    from PySide6.QtGui import QEnterEvent
    from PySide6.QtWidgets import QApplication

    app = QApplication([])
    from jcast.gui.main_window import MainWindow

    def spin(ms):
        loop = QEventLoop()
        QTimer.singleShot(ms, loop.quit)
        loop.exec()

    os.symlink(spec['data'], 'data')
    w = MainWindow()
    w.resize(980, 800)
    w.show()
    s = w.setup
    s.rmats_folder.setText('data/rmats')
    s.gtf_file.setText('data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf')
    s.genome.setText('data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
    s.uniprot.setText('data/uniprot/uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz')
    s.out.setText('screenshots out')
    s.cache_dir.setText('cache')
    s.canonical.setChecked(True)
    s.refresh()
    QApplication.sendEvent(s.uniprot, QEnterEvent(QPointF(1, 1), QPointF(1, 1), QPointF(1, 1)))
    spin(300)
    w.grab().save(os.path.join(spec['out'], 'setup_help.png'))

    w.tabs.setCurrentWidget(w.run)
    w.run.start.click()
    while w.run.running() and not (w.run.bars and w.run.bars['SE'].value() > 200):
        spin(20)
    QApplication.sendEvent(w.run.last, QEnterEvent(QPointF(1, 1), QPointF(1, 1), QPointF(1, 1)))
    spin(100)
    w.grab().save(os.path.join(spec['out'], 'run_progress.png'))
    while w.run.running():
        spin(50)
    app.quit()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before')
    parser.add_argument('--after')
    parser.add_argument('--build-python', default=sys.executable)
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_21')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_21')
    parser.add_argument('--probe', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe:
        gui_probe(json.loads(Path(args.probe).read_text()))
        return

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    args.out.mkdir(parents=True, exist_ok=True)
    cache_dir = args.work / 'cache'
    run(args.after, ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir, False)    # the GTF cache, once

    lines = []
    for side in ('before', 'after'):
        c = subprocess.run(['git', '-C', getattr(args, side), 'rev-parse', '--short', 'HEAD'], capture_output=True,
                           text=True)
        lines.append(f'{side}: {c.stdout.strip()}')

    lines += ['', '== chr15 runs (<inputs> -c --uniprot <args>): before and after']
    for label, extra in CASES.items():
        for progress in (False, True):
            name = f'{label}{" --progress" if progress else ""}'
            outs, steps = {}, {}
            for side in ('before', 'after'):
                outs[side] = args.work / f'{side} {name}'
                code, steps[side] = run(getattr(args, side), extra, outs[side], cache_dir, progress)
                assert code == 0, (side, name)
            lines.append(f'{name}: ' + compare(outs['before'], outs['after']))
            if progress:
                for side in ('before', 'after'):
                    lines.append(f'  progress, {side}: ' + progress_summary(steps[side], outs[side]))

    lines += ['', '== The packages built from the code after']
    dist = args.work / 'dist'
    subprocess.run([args.build_python, '-m', 'build', '--outdir', str(dist), str(Path(args.after).resolve())],
                   check=True, capture_output=True)
    for package in sorted(dist.iterdir()):
        if package.suffix == '.whl':
            members = {i.filename: i.file_size for i in zipfile.ZipFile(package).infolist()}
        else:
            members = {m.name.split('/', 1)[1]: m.size for m in tarfile.open(package).getmembers() if '/' in m.name}
        icon = members.get('jcast/gui/icon.png')
        lines.append(f'{package.name}: jcast/gui/icon.png {f"{icon:,} bytes" if icon else "MISSING"}')

    spec = args.work / 'spec.json'
    spec.write_text(json.dumps(dict(data=str(DATA), out=str(args.out.resolve()))))
    subprocess.run([sys.executable, __file__, '--probe', str(spec)], check=True, cwd=args.work, env=env(args.after),
                   capture_output=True)
    lines += ['', '== Screenshots: setup_help.png, run_progress.png']

    text = '\n'.join(lines).rstrip() + '\n'
    for path, name in ((args.work, '<work>'), (Path(args.before).resolve(), '<before>'),
                       (Path(args.after).resolve(), '<after>'), (REPO, '<repo>')):
        text = text.replace(str(path), name)
    print(text)
    (args.out / 'verify_21.txt').write_text(text)


if __name__ == '__main__':
    main()
