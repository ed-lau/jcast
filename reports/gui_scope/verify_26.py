# -*- coding: utf-8 -*-

"""Verification of PR 26 (reports/2026-11-13-gui-cli-scope.md §14): the GUI keeps the Setup page's settings between
sessions.

    code          the files the PR changes, against --before: none outside jcast/gui/ and tests/, so no run changes
    sessions      with the code after (--after), in Qt's native format (on macOS a .plist) in a file of the scratch
                  directory, never the user's: the first session's settings (the defaults); settings filled in and the
                  window closed; the second session's settings, the run's arguments, and the file's content; Restore
                  the defaults, then the third session's settings
    screenshot    the second session's Setup page, with the help line over Restore the defaults (paths relative to
                  the scratch directory, where data links to the test data)

The probe runs in a process of its own, from the scratch directory, with PYTHONPATH set to the checkout after.

Writes <out>/verify_26.txt and <out>/setup_kept.png.

Usage (from the repository root):
    python reports/gui_scope/verify_26.py --before <commit> --after <checkout> [--work <dir>] [--out <dir>]
"""

import argparse
import json
import os
import plistlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
FILLED = {'rmats_folder': 'data/rmats', 'gtf_file': 'data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf',
          'genome': 'data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz', 'out': 'screenshots out',
          'cache_dir': 'cache', 'uniprot': 'data/uniprot/uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz',
          'splice_types': 'RI,A5SS', 'q_lo': 0.0, 'q_hi': 0.05, 'canonical': True, 'processes': 2, 'verbose': True}


def probe(spec: dict) -> None:
    from PySide6.QtCore import QEventLoop, QPointF, QSettings, QTimer
    from PySide6.QtGui import QEnterEvent
    from PySide6.QtWidgets import QApplication, QMessageBox

    app = QApplication([])
    from jcast.gui import main_window
    from jcast.gui.setup_page import DEFAULTS

    def spin(ms):
        loop = QEventLoop()
        QTimer.singleShot(ms, loop.quit)
        loop.exec()

    found = {'user settings (not written)': main_window.user_settings().fileName()}
    main_window.user_settings = lambda: QSettings(spec['settings'], QSettings.NativeFormat)
    os.symlink(spec['data'], 'data')

    first = main_window.MainWindow()
    found['first session: the defaults'] = first.setup.values() == DEFAULTS
    first.setup.set_values(FILLED)
    first.close()

    second = main_window.MainWindow()
    second.resize(980, 800)
    second.show()
    found['second session: the settings filled in'] = second.setup.values() == FILLED
    found['second session: problems'] = second.setup.problems
    found['second session: the cutoff'] = list(second.model.options())
    found['second session: the run'] = vars(second.arguments())
    QApplication.sendEvent(second.setup.restore, QEnterEvent(QPointF(1, 1), QPointF(1, 1), QPointF(1, 1)))
    spin(300)
    second.grab().save(os.path.join(spec['out'], 'setup_kept.png'))

    QMessageBox.question = lambda *a, **k: QMessageBox.Yes
    second.setup.restore.click()
    second.close()
    third = main_window.MainWindow()
    found['third session, after Restore the defaults: the defaults'] = third.setup.values() == DEFAULTS
    third.close()
    print(json.dumps(found, default=sorted))
    app.quit()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before')
    parser.add_argument('--after')
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_26')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_26')
    parser.add_argument('--probe', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe:
        probe(json.loads(Path(args.probe).read_text()))
        return

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    args.out.mkdir(parents=True, exist_ok=True)

    changed = subprocess.run(['git', '-C', args.after, 'diff', '--name-only', args.before], capture_output=True,
                             text=True, check=True).stdout.split()
    outside = [f for f in changed if not f.startswith(('jcast/gui/', 'tests/', 'reports/', 'README', 'CHANGELOG',
                                                      'REVIEW', 'HANDOFF', 'jcast_website/'))]
    lines = [f'before: {args.before}', f'after: {args.after} (the working tree)', '',
             f'== code: {len(changed)} files changed; outside jcast/gui/, the tests and the documents: '
             f'{outside or "none"}', '']

    spec = {'settings': str(args.work / 'settings.plist'), 'data': str(DATA), 'out': str(args.out.resolve())}
    (args.work / 'spec.json').write_text(json.dumps(spec))
    env = dict(os.environ, PYTHONPATH=str(Path(args.after).resolve()), QT_QPA_PLATFORM='offscreen')
    probe_args = ['--probe', str(args.work / 'spec.json')]
    done = subprocess.run([sys.executable, '-P', str(Path(__file__).resolve()), *probe_args], capture_output=True,
                          text=True, cwd=args.work, env=env)
    if done.returncode:
        sys.exit(done.stderr)
    found = json.loads(done.stdout.strip().splitlines()[-1])
    lines.append('== sessions, the settings in Qt\'s native format, in <work>/settings.plist')
    lines += [f'{k}: {v}' for k, v in found.items()]
    lines += ['', '== <work>/settings.plist after the third session']
    lines += [f'{k}: {v!r}' for k, v in sorted(plistlib.loads((args.work / 'settings.plist').read_bytes()).items())]

    text = '\n'.join(lines).replace(str(args.work), '<work>').replace(str(REPO), '<repo>') + '\n'
    text = text.replace(str(Path.home()), '~')
    print(text)
    (args.out / 'verify_26.txt').write_text(text)


if __name__ == '__main__':
    main()
