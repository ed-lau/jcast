# -*- coding: utf-8 -*-

"""Verification of PR 13 (reports/2026-11-13-gui-cli-scope.md §12): the GUI.

    pipeline      what PR 13 changes in the pipeline (outcomes.fate_table; model.draw_model without a fit, returning the
                  cutoff's line and label; model.move_cutoff), with the code before and after: chr15 runs, output files
                  byte-identical, logs identical but for timestamps, the output directory and the time taken; and
                  model.pdf of each bundled rMATS folder with a fit, drawn as PNG, identical.
    gui runs      with the code after: runs started from the GUI, without a display, and their command line, as the
                  Run page shows it, run with another -o: the same output files and logs.
    stop          a GUI run stopped while it translates: its exit status, the log's last line, the worker processes,
                  and the files left after Delete.
    screenshots   the four pages after a -m run, saved in <out> (paths relative to the scratch directory, where data
                  links to the test data).
    model page    each bundled rMATS folder on the Model page: what it shows, and the time to fit and draw.
    whole genome  (--rmats, --output) the Model page's times on a whole genome: the fit, the plot, a cutoff moved,
                  the table; and the Results page's time to read a whole-genome output folder.

Each process runs from a scratch directory, with PYTHONPATH set to its checkout: python -m and -c put the current
directory on sys.path first. The GUI runs offscreen (QT_QPA_PLATFORM=offscreen).

Writes <out>/verify_13.txt and <out>/<page>.png.

Usage (from the repository root):
    python reports/gui_scope/verify_13.py --before <checkout> --after <checkout> [--rmats <dir>] [--output <dir>]
        [--work <dir>] [--out <dir>]
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
RMATS = DATA / 'rmats'
GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'
FOLDERS = [RMATS] + sorted(p for p in (DATA / 'model').iterdir() if p.is_dir())
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
SET = re.compile(r"splice_type=\{([^}]*)\}")

# pipeline: chr15, before and after
PIPELINE = {'-c --uniprot': ['-c', '--uniprot', UNIPROT], '-m --uniprot': ['-m', '--uniprot', UNIPROT]}

# gui runs: the Setup page's fields, and the Model page's cutoff ('model', or a click on the plot at a read count)
GUI_RUNS = {
    'all splice types, --uniprot': dict(uniprot=True),
    '-m (the Model page), --uniprot -c': dict(uniprot=True, canonical=True, cutoff='model'),
    '-r 10 (a click on the plot), -s RI A5SS, -q 0 0.05, --uniprot': dict(uniprot=True, cutoff=10,
                                                                         splice_types=['RI', 'A5SS'], q_hi=0.05),
    '--uniprot -c -p 2 -v': dict(uniprot=True, canonical=True, processes=2, verbose=True),
}


def env(repo, **extra):
    return dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg', COLUMNS='1000',
                QT_QPA_PLATFORM='offscreen', **extra)


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(out).iterdir())
            if p.is_file() and p.name not in NOT_COMPARED}


def log_lines(out):
    """ the log's lines, without timestamps or the time taken; the output directory written <out>; -s's default set,
    whose order changes with each process, sorted """
    text = (Path(out) / 'logfile.log').read_text().replace(str(out), '<out>')
    text = SET.sub(lambda m: 'splice_type={' + ', '.join(sorted(m.group(1).split(', '))) + '}', text)
    return [STAMP.sub('', line) for line in text.splitlines() if not STAMP.sub('', line).startswith('INFO~Time: ')]


def compare(a, b) -> str:
    fa, fb = files(a), files(b)
    differ = sorted(n for n in set(fa) | set(fb) if fa.get(n) != fb.get(n))
    la, lb = log_lines(a), log_lines(b)
    return (f'{len(fa)} output files, identical: {"yes" if not differ else "NO, " + ", ".join(differ)}; '
            f'log identical but for timestamps, the output directory and the time taken: '
            f'{"yes" if la == lb else "NO"}')


def cli(repo, extra, out, cache_dir):
    """ python -P -m jcast run, as the GUI runs it """
    out.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, '-P', '-m', 'jcast', 'run', RMATS, GTF, GENOME, '-o', out, '--cache_dir', cache_dir] + extra
    started = time.monotonic()
    done = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=out, env=env(repo))
    return done.returncode, round(time.monotonic() - started, 1)


# ---- the probes, run in a checkout's process ----

PLOT_PROBE = r'''
import hashlib, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.figure
from jcast import model

png_dir, folders = sys.argv[1], sys.argv[2:]
saved = []
def savefig(self, fname, **kwargs):
    path = f'{png_dir}/{len(saved)}.png'
    self.canvas.print_figure(path, dpi=100, format='png', metadata={'Software': None})
    saved.append(path)
matplotlib.figure.Figure.savefig = savefig

import jcast
found = {'jcast': jcast.__file__}
for folder in folders:
    _, ln_sjc, params, cutoff, _ = model.preview(folder)
    if params is not None:
        model.plot_general_mixture_model(ln_sjc, params, cutoff, write_dir=png_dir, filename='model')
        found[folder] = hashlib.sha256(open(saved[-1], 'rb').read()).hexdigest()
print(json.dumps(found))
'''


def gui_probe(spec: dict) -> dict:
    """ In the checkout's process: drive the GUI offscreen as spec says """
    import math
    import subprocess as sp
    from types import SimpleNamespace

    from PySide6.QtCore import QEventLoop, QTimer
    from PySide6.QtWidgets import QApplication, QMessageBox

    app = QApplication([])
    import jcast
    from jcast import main
    from jcast.gui import command
    from jcast.gui.main_window import MainWindow
    from jcast.gui.results_page import read_output
    from jcast.gui.run_page import RunPage

    QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Yes)
    stages = []
    shown = RunPage._show_stage
    RunPage._show_stage = lambda self: stages.append(self._stage) or shown(self)

    def spin(ms):
        loop = QEventLoop()
        QTimer.singleShot(ms, loop.quit)
        loop.exec()

    def until(condition, timeout=600):
        started = time.monotonic()
        while not condition():
            if time.monotonic() - started > timeout:
                raise TimeoutError
            spin(20)

    def window(rmats, out, case):
        w = MainWindow()
        w.show()
        s = w.setup
        s.rmats_folder.setText(rmats)
        s.gtf_file.setText(spec['gtf'])
        s.genome.setText(spec['genome'])
        s.out.setText(out)
        s.cache_dir.setText(spec['cache_dir'])
        s.uniprot.setText(spec['uniprot'] if case.get('uniprot') else '')
        s.refresh()
        for t, box in s.splice_types.items():
            box.setChecked(t in case.get('splice_types', ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']))
        s.processes.setValue(case.get('processes', 1))
        s.q_hi.setValue(case.get('q_hi', 1))
        s.canonical.setChecked(case.get('canonical', False))
        s.verbose.setChecked(case.get('verbose', False))
        s.refresh()
        return w

    def fitted(w):
        w.tabs.setCurrentWidget(w.model)
        until(lambda: w.model._fit is not None or 'could not' in w.model.status.text())

    found = {'jcast': jcast.__file__, 'child PYTHONPATH': None, 'runs': {}, 'model': {}}
    for label, case in spec['runs'].items():
        out = os.path.join(spec['work'], 'gui ' + label)
        w = window(spec['rmats'], out, case)
        if 'cutoff' in case:
            fitted(w)
            if case['cutoff'] == 'model':
                w.model.use_model.setChecked(True)
            else:
                w.model._clicked(SimpleNamespace(name='button_press_event', xdata=math.log(case['cutoff'] + 1),
                                                 button=1))
        found['child PYTHONPATH'] = command.process(w.arguments())[2].split(os.pathsep)[0]
        del stages[:]
        started = time.monotonic()
        w.tabs.setCurrentWidget(w.run)
        w.run.start.click()
        until(lambda: not w.run.running())
        took = round(time.monotonic() - started, 1)
        until(lambda: w.results.shown is not None)
        found['runs'][label] = dict(out=out, command=w.run.command.toPlainText(), arguments=command.arguments(
            w.run.args), exit=w.run.process.exitCode(), took=took, status=w.run.status.text(),
            stages=list(dict.fromkeys(s for s in stages if s != 'Starting')),
            bars={t: [b.value(), b.maximum()] for t, b in w.run.bars.items()}, results=w.results.status.text())
        w.close()

    if spec.get('stop'):
        out = os.path.join(spec['work'], 'gui stopped')
        w = window(spec['rmats'], out, dict(uniprot=True, processes=2))
        w.run.start.click()
        until(lambda: any(b.value() > 0 for b in w.run.bars.values()))
        pid = w.run.process.processId()
        workers = sp.run(['pgrep', '-P', str(pid)], capture_output=True, text=True).stdout.split()
        at = {t: b.value() for t, b in w.run.bars.items() if b.value()}
        w.run.stop_run(ask=False)
        until(lambda: not w.run.running())
        left = [p for p in workers if sp.run(['kill', '-0', p], capture_output=True).returncode == 0]
        log = open(os.path.join(out, 'logfile.log')).read().splitlines()
        written = main.existing_outputs(out)
        delete = not w.run.delete.isHidden()
        w.run._delete_output()
        found['stop'] = dict(at=at, exit=w.run.process.exitCode(), status=w.run.status.text(), workers=len(workers),
                             left=len(left), last=STAMP.sub('', log[-1]), written=written, delete_offered=delete,
                             after_delete=sorted(os.listdir(out)))
        w.close()

    if spec.get('screenshots'):
        os.symlink(spec['data'], 'data')
        w = MainWindow()
        w.resize(980, 760)
        w.show()
        s = w.setup
        s.rmats_folder.setText('data/rmats')
        s.gtf_file.setText('data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf')
        s.genome.setText('data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz')
        s.uniprot.setText('data/uniprot/uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz')
        s.out.setText('out')
        s.cache_dir.setText('cache')
        s.canonical.setChecked(True)
        s.refresh()
        s.processes.setValue(2)
        s.refresh()
        fitted(w)
        w.model.use_model.setChecked(True)
        w.tabs.setCurrentWidget(w.run)
        w.run.start.click()
        until(lambda: not w.run.running())
        until(lambda: w.results.shown is not None)
        for name, page in (('setup', w.setup), ('model', w.model), ('run', w.run), ('results', w.results)):
            w.tabs.setCurrentWidget(page)
            spin(300)
            w.grab().save(os.path.join(spec['screenshots'], f'{name}.png'))
        w.close()

    for folder in spec['folders']:
        w = window(folder, os.path.join(spec['work'], 'model'), {})
        started = time.monotonic()
        fitted(w)
        found['model'][folder] = dict(status=w.model.status.text(), took=round(time.monotonic() - started, 2),
                                      use_model=w.model.use_model.isEnabled())
        w.close()

    if spec.get('whole_rmats'):
        w = window(spec['whole_rmats'], os.path.join(spec['work'], 'model'), {})
        started = time.monotonic()
        fitted(w)
        fit = time.monotonic() - started
        page = w.model
        started = time.monotonic()
        page._draw()
        page.canvas.draw()
        draw = time.monotonic() - started
        moves = []
        for cutoff in (10, 20, 40, 80, 160):
            started = time.monotonic()
            page.read.setValue(cutoff)
            page.canvas.draw()
            moves.append(time.monotonic() - started)
        started = time.monotonic()
        page.update_table()
        table = time.monotonic() - started
        found['whole genome'] = dict(status=page.status.text(), fit=round(fit, 1), draw=round(draw, 2),
                                     move=round(max(moves), 3), table=round(table, 2))
        w.close()
    if spec.get('whole_output'):
        started = time.monotonic()
        result = read_output(spec['whole_output'])
        read = time.monotonic() - started
        w = MainWindow()
        started = time.monotonic()
        w.results._show(result)
        fill = time.monotonic() - started
        found['whole output'] = dict(read=round(read, 1), fill=round(fill, 3), status=w.results.status.text(),
                                     junctions=result['fates'][-1][-1], records=sum(result['records'].values()))
    app.quit()
    return found


def run_probe(repo, spec: dict, work: Path) -> dict:
    path = work / 'spec.json'
    path.write_text(json.dumps(spec))
    done = subprocess.run([sys.executable, __file__, '--probe', str(path)], capture_output=True, text=True, cwd=work,
                          env=env(repo))
    if done.returncode:
        raise RuntimeError(done.stderr)
    return json.loads(done.stdout.splitlines()[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before')
    parser.add_argument('--after')
    parser.add_argument('--rmats', type=Path, help='a whole-genome rMATS folder, for the Model page')
    parser.add_argument('--output', type=Path, help='a whole-genome JCAST output folder, for the Results page')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_13')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_13')
    parser.add_argument('--probe', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe:
        print(json.dumps(gui_probe(json.loads(Path(args.probe).read_text()))))
        return

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    cache_dir = args.work / 'cache'
    cli(args.after, ['-s', 'RI', '-r', '100000'], args.work / 'warm', cache_dir)     # the GTF cache, written once

    lines = []
    for side in ('before', 'after'):
        c = subprocess.run(['git', '-C', getattr(args, side), 'rev-parse', '--short', 'HEAD'], capture_output=True,
                           text=True)
        lines.append(f'{side}: {c.stdout.strip()}')

    # ---- pipeline ----
    lines += ['', '== The pipeline, chr15 (<inputs> <args>): before and after']
    for label, extra in PIPELINE.items():
        (rb, _), (ra, _) = [cli(getattr(args, side), extra, args.work / f'{side} {label}', cache_dir)
                            for side in ('before', 'after')]
        lines.append(f'{label}: exit {rb} / {ra}; '
                     + compare(args.work / f'before {label}', args.work / f'after {label}'))
    plots = {}
    for side in ('before', 'after'):
        png = args.work / f'png {side}'
        png.mkdir()
        done = subprocess.run([sys.executable, '-c', PLOT_PROBE, str(png)] + [str(f) for f in FOLDERS], check=True,
                              capture_output=True, text=True, cwd=png, env=env(getattr(args, side)))
        plots[side] = json.loads(done.stdout.splitlines()[-1])
    lines.append(f'model.pdf as PNG, the bundled folders with a fit (jcast from {plots["before"].pop("jcast")} / '
                 f'{plots["after"].pop("jcast")}):')
    for folder, digest in plots['before'].items():
        lines.append(f'  {Path(folder).name}: {"identical" if plots["after"].get(folder) == digest else "DIFFERENT"}')

    # ---- the GUI, with the code after ----
    args.out.mkdir(parents=True, exist_ok=True)
    spec = dict(work=str(args.work), rmats=str(RMATS), gtf=str(GTF), genome=str(GENOME), uniprot=str(UNIPROT),
                cache_dir=str(cache_dir), runs=GUI_RUNS, stop=True, folders=[str(f) for f in FOLDERS],
                screenshots=str(args.out.resolve()), data=str(DATA),
                whole_rmats=str(args.rmats) if args.rmats else None,
                whole_output=str(args.output) if args.output else None)
    found = run_probe(args.after, spec, args.work)
    lines += ['', f'== GUI runs (jcast from {found["jcast"]}; the child process PYTHONPATH starts with '
                  f'{found["child PYTHONPATH"]}), and their command line run with another -o']
    for label, run in found['runs'].items():
        out = args.work / f'cli {label}'
        arguments = run['arguments']
        arguments[arguments.index('-o') + 1] = str(out)
        out.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, '-P', '-m', 'jcast', 'run'] + arguments
        started = time.monotonic()
        code = subprocess.run(cmd, capture_output=True, cwd=out, env=env(args.after)).returncode
        took = round(time.monotonic() - started, 1)
        lines += [f'{label}:',
                  f'  command: {run["command"]}',
                  f'  GUI exit {run["exit"]} in {run["took"]} s, command line exit {code} in {took} s; '
                  + compare(run['out'], out),
                  f'  stages shown: {", ".join(run["stages"])}',
                  f'  progress: {", ".join(f"{t} {v}/{m}" for t, (v, m) in run["bars"].items())}',
                  f'  Run page: {run["status"]}',
                  f'  Results page: {run["results"]}']

    stop = found['stop']
    lines += ['', '== A GUI run stopped while translating (every chr15 junction, --uniprot -p 2)',
              f'stopped at {stop["at"]}; exit status {stop["exit"]}; {stop["workers"]} worker processes, '
              f'{stop["left"]} left after the stop',
              f'Run page: {stop["status"]}',
              f'the log\'s last line: {stop["last"]}',
              f'files written: {", ".join(stop["written"])}; Delete offered: {stop["delete_offered"]}; '
              f'left after Delete: {stop["after_delete"] or "none"}']

    lines += ['', '== The Model page, each bundled rMATS folder: what it shows, and the time to fit and draw']
    for folder, model in found['model'].items():
        lines.append(f'{Path(folder).name} ({model["took"]} s; -m {"enabled" if model["use_model"] else "disabled"}): '
                     f'{model["status"]}')

    if 'whole genome' in found:
        w = found['whole genome']
        lines += ['', f'== The Model page on a whole genome ({args.rmats})', w['status'],
                  f'fit, in the background: {w["fit"]} s; the plot drawn: {w["draw"]} s; a cutoff moved and drawn: '
                  f'{w["move"]} s at most; the table: {w["table"]} s']
    if 'whole output' in found:
        w = found['whole output']
        lines += ['', f'== The Results page on a whole-genome output folder ({args.output})', w['status'],
                  f'{w["junctions"]:,} junctions, {w["records"]:,} records; read, in the background: {w["read"]} s; '
                  f'the page filled: {w["fill"]} s']

    text = '\n'.join(lines).rstrip() + '\n'
    for path, name in ((args.work, '<work>'), (Path(args.before).resolve(), '<before>'),
                       (Path(args.after).resolve(), '<after>'), (REPO, '<repo>')):
        text = text.replace(str(path), name)
    print(text)
    (args.out / 'verify_13.txt').write_text(text)


if __name__ == '__main__':
    main()
