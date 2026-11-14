# -*- coding: utf-8 -*-

"""Verification of PR 18 (reports/2026-11-13-gui-cli-scope.md §9): the command line on Typer.

Runs `python -m jcast <args>` with the code before and after (PYTHONPATH set to each checkout, so the code under test
can be frozen in a git worktree), on the chr15 test set:

    same        runs that worked before, as written before (no command name, -s SE MXE): every output file but the
                log and model.pdf is byte-identical, and so are the cutoff, the log's summary and its other lines; the
                Arguments line holds the same arguments, less --mask and func.
    forms       after only: `jcast run …`, and -s SE -s MXE, give the output of the forms above.
    errors      command lines that fail, or that changed: the exit status before and after, and the message.

Writes <out>/verify_18.txt.

Usage (from the repository root):
    python reports/cli_checks/verify_18.py --before <checkout> --after <checkout> [--work <dir>] [--out <dir>]
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / 'tests' / 'data'
RMATS = DATA / 'rmats'
GTF = DATA / 'genome' / 'Homo_sapiens.GRCh38.89.chromosome.15.gtf'
GENOME = DATA / 'genome' / 'Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz'
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'
INPUTS = [RMATS, GTF, GENOME, '-c', '--uniprot', UNIPROT]
NOT_COMPARED = ('logfile.log', 'model.pdf')
STAMP = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~')
CUTOFF = re.compile(r'Translating junctions with a read count of at least (\d+)')
FUNC = re.compile(r', func=<function run_jcast at 0x[0-9a-f]+>')

# runs that worked before, as written before: <inputs> <args>
SAME = {'default': [],
        '-m': ['-m'],
        '-s SE MXE -m': ['-s', 'SE', 'MXE', '-m'],
        '-r 30 -q 0 0.05 -s RI A5SS': ['-r', '30', '-q', '0', '0.05', '-s', 'RI', 'A5SS'],
        '-p 2': ['-p', '2']}

# after only: (label, the command line, the SAME run it should equal)
FORMS = [('run <inputs> -m', ['run'] + INPUTS + ['-m'], '-m'),
         ('<inputs> -s SE -s MXE -m', INPUTS + ['-s', 'SE', '-s', 'MXE', '-m'], '-s SE MXE -m')]

# a quick run: one splice type, every junction skipped by its read count
QUICK = [RMATS, GTF, GENOME, '-s', 'RI', '-r', '100000']
ERRORS = {'no arguments': [],
          '--version': ['--version'],
          'missing GTF': [RMATS, 'missing.gtf', GENOME],
          'missing rMATS folder': ['missing_rmats', GTF, GENOME],
          'missing rMATS file of -s': ['<se_only>', GTF, GENOME, '-s', 'SE', 'MXE'],
          '-p 2 without --uniprot': [RMATS, GTF, GENOME, '-p', '2'],
          '-p 0': [RMATS, GTF, GENOME, '-p', '0'],
          '-q 0.5 0.2': [RMATS, GTF, GENOME, '-q', '0.5', '0.2'],
          '-q 0 1.5': [RMATS, GTF, GENOME, '-q', '0', '1.5'],
          '-s XX': [RMATS, GTF, GENOME, '-s', 'XX'],
          '--g_or_ln 1': [RMATS, GTF, GENOME, '--g_or_ln', '1'],
          'output directory holding JCAST output': [RMATS, GTF, GENOME, '-o', '<used>'],
          'missing GTF and genome, -p 0': [RMATS, 'missing.gtf', 'missing.fa', '-p', '0'],
          '--mask N (quick run)': QUICK + ['--mask', 'N'],
          '--overw, abbreviated (quick run)': QUICK + ['--overw'],
          '-s SE -s MXE (quick run)': [RMATS, GTF, GENOME, '-s', 'SE', '-s', 'MXE', '-r', '100000']}


def run(repo, argv, out, cache_dir, work):
    """ the run's exit status, its stdout and stderr, and its output directory """
    argv = [str(a).replace('<se_only>', str(work / 'se_only')).replace('<used>', str(work / 'used')) for a in argv]
    if argv and not argv[0].startswith('-'):
        argv += ['-o', str(out)] if '-o' not in argv else []
        argv += ['--cache_dir', str(cache_dir)]
    out.mkdir(parents=True, exist_ok=True)
    done = subprocess.run([sys.executable, '-m', 'jcast'] + argv, cwd=out, capture_output=True, text=True,
                          env=dict(os.environ, PYTHONPATH=str(Path(repo).resolve()), MPLBACKEND='Agg', COLUMNS='1000',
                                   NO_COLOR='1'))
    return done.returncode, done.stdout, done.stderr, out


def message(stdout, stderr) -> str:
    """ the error, from argparse (<prog>: error: …; run as python -m jcast, its prog is __main__.py) or from Typer's
    box; or the warnings; or the first line printed """
    lines = [line.strip(' │╭╮╰╯─') for line in stderr.splitlines()]
    errors = [line.split(': error: ', 1)[1] for line in lines if re.match(r'^\S+: error: ', line)]
    if not errors:
        start = next((i for i, line in enumerate(lines) if line == 'Error'), None)
        errors = [line for line in lines[start + 1:] if line] if start is not None else []
    if errors:
        return ' / '.join(errors)
    warnings = [line for line in lines if line.startswith('WARNING -- ')]
    first = next((line.strip() for line in stdout.splitlines() if line.strip()), '')
    return ' / '.join(warnings) or (f'stdout: {first}' if first else '(nothing printed)')


def files(out):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
            if p.is_file() and p.name not in NOT_COMPARED}


def log_lines(out):
    """ the log's lines, without their timestamps, with the output directory written <out> """
    text = (out / 'logfile.log').read_text().replace(str(out), '<out>')
    return [STAMP.sub('', line) for line in text.splitlines()]


def arguments(out) -> dict:
    """ the run's arguments, from the log's Arguments line """
    line = next(line for line in log_lines(out) if line.startswith('INFO~Arguments: '))
    return vars(eval(FUNC.sub('', line.split('Arguments: ', 1)[1]), {'Namespace': argparse.Namespace}))


def shown(value) -> str:
    """ the value as the log prints it; a set (-s's default) in sorted order, since its order changes with each run """
    return '{' + repr(sorted(value))[1:-1] + '}' if isinstance(value, set) else repr(value)


def summary(out):
    lines = log_lines(out)
    start = next(i for i, l in enumerate(lines) if l.endswith('Summary.'))
    end = next(i for i, l in enumerate(lines) if l.startswith('INFO~Time: '))
    return lines[start:end]


def cutoff(out):
    found = CUTOFF.search((out / 'logfile.log').read_text())
    return int(found.group(1)) if found else None


def commit(repo):
    c = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    d = subprocess.run(['git', '-C', str(repo), 'status', '--porcelain', '--', 'jcast'], capture_output=True, text=True)
    return c.stdout.strip() + (' + changes to jcast/' if d.stdout.strip() else '')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--work', type=Path, default=Path(tempfile.gettempdir()) / 'jcast_verify_18')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'verify_18')
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    cache_dir = args.work / 'cache'
    (args.work / 'se_only').mkdir(parents=True)
    shutil.copy(RMATS / 'SE.MATS.JC.txt', args.work / 'se_only')
    (args.work / 'used').mkdir()
    (args.work / 'used' / 'logfile.log').write_text('')

    # warm the GTF cache once, so that the parallel runs do not all build it
    run(args.after, QUICK, args.work / 'warm', cache_dir, args.work)

    jobs = {}
    for label, extra in SAME.items():
        for side in ('before', 'after'):
            jobs[('same', label, side)] = (getattr(args, side), INPUTS + extra)
    for label, argv, _ in FORMS:
        jobs[('forms', label, 'after')] = (args.after, argv)
    for label, argv in ERRORS.items():
        for side in ('before', 'after'):
            jobs[('errors', label, side)] = (getattr(args, side), argv)

    def one(key):
        repo, argv = jobs[key]
        out = args.work / re.sub(r'[^A-Za-z0-9]+', '_', '-'.join(key))
        return key, run(repo, argv, out, cache_dir, args.work)

    with ThreadPoolExecutor(4) as pool:
        results = dict(pool.map(one, jobs))

    lines = [f'before: {commit(args.before)}; after: {commit(args.after)}',
             '<inputs>: <chr15 rMATS> <chr15 GTF> <chr15 genome> -c --uniprot <chr15 Swiss-Prot>', '',
             '== Runs that worked before, as written before: before and after']
    for label in SAME:
        (rb, ob, eb, b), (ra, oa, ea, a) = results[('same', label, 'before')], results[('same', label, 'after')]
        if rb or ra:
            lines.append(f'<inputs> {label}: exit {rb} / {ra}: {message(ob, eb)} / {message(oa, ea)}')
            continue
        fb, fa = files(b), files(a)
        differ = sorted(n for n in set(fb) | set(fa) if fb.get(n) != fa.get(n))
        lines.append(f'<inputs> {label}: {len(fa)} output files, identical: '
                     f'{"yes" if not differ else "no, " + ", ".join(differ)}; cutoff {cutoff(b)} / {cutoff(a)}; '
                     f'summary identical: {"yes" if summary(b) == summary(a) else "no"}')
        before, after = arguments(b), arguments(a)
        gone = {k: before.pop(k) for k in ('mask', 'func') if k in before}
        lines.append(f'  Arguments: the same{" values" if before == after else ": no"}, in the same order: '
                     f'{"yes" if list(before) == list(after) else "no"}; before only: {gone}; printed differently: '
                     f'{[f"{k}: {shown(before[k])} / {shown(after[k])}" for k in after if shown(before.get(k)) != shown(after[k])]}')
        skip = ('INFO~Arguments: ', 'INFO~Time: ')
        lb, la = log_lines(b), log_lines(a)
        lines += [f'  log, before only: {l}' for l in lb if l not in la and not l.startswith(skip)]
        lines += [f'  log, after only:  {l}' for l in la if l not in lb and not l.startswith(skip)]

    lines += ['', '== Forms new in PR 18 (after only), against the run above']
    for label, _, same_as in FORMS:
        r, o, e, a = results[('forms', label, 'after')]
        _, _, _, other = results[('same', same_as, 'after')]
        identical = not r and files(a) == files(other) and summary(a) == summary(other) and \
            arguments(a) == arguments(other)
        lines.append(f'{label}: exit {r}; output, summary and arguments identical to <inputs> {same_as}: '
                     f'{"yes" if identical else "no"}')

    lines += ['', '== Command lines that fail, or changed: exit status before / after, and the message']
    for label in ERRORS:
        (rb, ob, eb, b), (ra, oa, ea, a) = results[('errors', label, 'before')], results[('errors', label, 'after')]
        lines.append(f'{label}: exit {rb} / {ra}')
        lines.append(f'  before: {message(ob, eb)}')
        lines.append(f'  after:  {message(oa, ea)}')
        if (b / 'logfile.log').exists() and (a / 'logfile.log').exists():
            sb, sa = arguments(b)['splice_type'], arguments(a)['splice_type']
            lines.append(f'  -s: {sb} / {sa}')

    text = ('\n'.join(lines).rstrip() + '\n').replace(str(args.work), '<work>').replace(str(REPO), '<repo>')
    print(text)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'verify_18.txt').write_text(text)


if __name__ == '__main__':
    main()
