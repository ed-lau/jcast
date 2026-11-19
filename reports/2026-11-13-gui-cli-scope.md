# GUI and command-line scope for v0.5.0 (PRs 18, 19, 13)

| | |
|---|---|
| **Date** | 2026-11-13; §14 on 2026-11-19 |
| **Status** | Decisions taken 2026-11-13 (§7). To be implemented in PRs 18, 19 and 13, in that order (§8). **PR 18 implemented in `5479ee3`** (§9); **PR 19 implemented in `ec1068a`** (§10). No output change. **PR 20 in `2f0e54c`** (§11): `-m` without a fit falls back to `-r`. **PR 13 implemented in `12e0576`** (§12): the GUI runs JCAST; no output change. **PR 21 in `4cc72c2`** (§13): a help line, progress counts and the last junction done, the logo. CI's first run (2026-11-14) failed 7 command-line tests: fixed in the tests (§9.4). **2026-11-14: v0.5.0 widened** before its release: §4's junction browser moves into it, as a view of each junction and its peptides, after a feasibility analysis; and the `-m` model is assessed ([REVIEW.md §9](../REVIEW.md#v050-the-gui-and-the-command-line), rows 22 and 23). 2026-11-19: the view moves to v0.6.0. **PR 26 in `4fef85b`** (§14): the GUI keeps the Setup page's settings between sessions; no output change. |
| **Question** | What should the v0.5.0 GUI do: start a run only, or also show results, and which? How is it launched? What does moving the command line to Typer change for users? |
| **Reproduce** | `python reports/gui_scope/scope_evidence.py --work <scratch> [--python <interpreter for the venvs>] [--rmats <whole-genome rMATS folder>]` (≈ 2 min; needs PyPI) |
| **Evidence** (code at `f37faf9`) | [gui_scope/scope_evidence/scope_evidence.txt](gui_scope/scope_evidence/scope_evidence.txt) |
| **PR 18 verification** | `python reports/cli_checks/verify_18.py --before <checkout> --after <checkout>` (≈ 30 s); [cli_checks/verify_18/verify_18.txt](cli_checks/verify_18/verify_18.txt) |
| **PR 19 verification** | `python reports/gui_scope/verify_19.py --before <checkout> --after <checkout> [--rmats <whole-genome rMATS folder>]` (≈ 1.5 min); [gui_scope/verify_19/verify_19.txt](gui_scope/verify_19/verify_19.txt) |
| **PR 20 verification** | `python reports/gui_scope/verify_20.py --before <checkout> --after <checkout> [--earlier <checkout before PR 19>]` (≈ 1 min); [gui_scope/verify_20/verify_20.txt](gui_scope/verify_20/verify_20.txt) |
| **PR 13 verification** | `python reports/gui_scope/verify_13.py --before <checkout> --after <checkout> [--rmats <whole-genome rMATS folder>] [--output <whole-genome output folder>]` (≈ 3 min; needs PySide6); [gui_scope/verify_13/verify_13.txt](gui_scope/verify_13/verify_13.txt) |
| **PR 21 verification** | `python reports/gui_scope/verify_21.py --before <checkout> --after <checkout> [--build-python <Python with build>]` (≈ 1.5 min; needs PySide6); [gui_scope/verify_21/verify_21.txt](gui_scope/verify_21/verify_21.txt) |
| **PR 26 verification** | `python reports/gui_scope/verify_26.py --before <commit> --after <checkout>` (≈ 5 s; needs PySide6); [gui_scope/verify_26/verify_26.txt](gui_scope/verify_26/verify_26.txt) |

## 1. Summary

- **The GUI released in 0.4.0 cannot run JCAST.** It has three placeholder
  pages ([jcast/gui/](../jcast/gui/), 157 lines).
- **Decided (§7):**
  - `jcast gui` opens the GUI. `jcast` becomes a group of commands, `run`
    and `gui`, and `jcast RMATS GTF GENOME …` still runs as before.
  - The GUI runs the command line in a child process,
    `python -m jcast run …`, and shows that command. A GUI run is a
    command-line run.
  - The GUI shows results, within limits:
    - before the run, the `-m` model, where the cutoff can be set;
    - after it, the run's summary, for this run or any earlier one.

    A junction browser waits for v0.5.x.
  - The command line moves to Typer (`typer>=0.26,<1`). `-s SE MXE` keeps
    working.
  - No standalone app bundles in v0.5.0.
- **Two findings:**
  - Typer cannot parse `-s SE MXE`: it exits with a usage error (§5).
  - The model preview takes about 27 s on the whole genome, not seconds.
    Most of it is computing the read counts junction by junction (§6).

## 2. Launching: `jcast gui`

Today the GUI opens with `python -m jcast --gui` (checked in
`__main__.py` before argparse) or the `jcast-gui` script. The `jcast`
script itself rejects `--gui`.

With Typer, `jcast` becomes a group of two commands. A group takes a
command name first, so `jcast A B C` fails with "No such command 'A'"
(§5). `main()` therefore puts `run` in front when the first argument is
not a command name or a help option:
- `jcast gui`, `jcast run RMATS GTF GENOME …` and
  `jcast RMATS GTF GENOME …` all work.
- `jcast -h` lists the commands; `jcast run -h` gives the options.
- The planned UniProt download helper can become a third command.
- An rMATS folder named `gui` or `run` must be written `./gui`. Otherwise
  the command fails with a usage error; it does not run something else.

`jcast-gui` and `python -m jcast --gui` stay: the README documents both.

## 3. How the GUI runs JCAST

| | In the GUI's process (`run_jcast` in a `QThread`) | A child process (`python -m jcast run …` in a `QProcess`) |
|---|---|---|
| `-p` | Forking a process that runs Qt threads is unsafe ([REVIEW.md §10](../REVIEW.md#10-gui-plan--notes)). `spawn` or `forkserver` would copy the genome and GTF into every worker instead of sharing them. In practice, one process. | As on the command line: the child runs no Qt threads. |
| Cancel | Needs a flag checked at each junction in `run_jcast`. | Stop the process. A signal handler lets the run end its workers and log the cancel. |
| The same as the command line | The same function, with arguments the GUI builds. | The same command, shown in the GUI, and the same log. |
| A crash, or running out of memory | Takes the GUI down. | Ends the run only. |
| Progress | A callback. | A machine-readable progress output (a hidden option). |
| Genome and GTF | Could stay loaded between runs. | Read again at each run: about 45 s with GENCODE v50 and GRCh38 (the GTF cache 8 s, the genome 35 s). |
| Standalone bundle | No change. | Needs its own way to start itself. |

The child process keeps the `-p` speed-up (about 11× with 16 processes on
the whole genome; [parallel audit §9](2026-11-12-parallel-audit.md)),
and a GUI run is reproducible from the command it shows. Of the three GUI
constraints in the v0.4.0 handoff, it settles two: no fork from the Qt
process, and cancellation. The model plot drawn on an `Axes` the caller
supplies is still needed.

## 4. What the GUI shows

Four pages, as tabs rather than a one-way wizard: a rerun with another
cutoff means going back.

1. **Setup.**
   - The rMATS folder, GTF, genome, `--uniprot` FASTA and output folder.
   - The splice types, as checkboxes; those whose file is not in the
     folder are disabled.
   - `-q`, `-c`, `-p` (only with `--uniprot`), `-v` and `--cache_dir`.
   - The checks of the command line, from the same function, shown next
     to the fields. An output folder holding JCAST output asks before
     replacing it (`--overwrite`).
2. **Model, before the run.** The command line shows `model.pdf` only
   after the run.
   - The model needs only the rMATS tables, not the GTF or the genome
     (§6).
   - It shows the histogram, the fit and the cutoff, and how many
     junctions pass, per splice type, within the `-q` range.
   - Accepting the model's cutoff runs with `-m`. Setting another runs
     with `-r N`. A number field sets the cutoff; dragging the line on
     the plot is a convenience.
   - The preview runs the code the run does: with `-m`, every rMATS file
     in the folder, also those `-s` leaves out
     ([command-line checks §3](2026-11-13-command-line-checks.md#3-decisions-2026-11-13)).
     A test checks that its cutoff equals the command line's: 26 on the
     chr15 folder.
3. **Run.**
   - The command, with a copy button; start and cancel.
   - Progress per splice type, warnings as they come, and the time taken.
   - A cancelled run is marked incomplete, with an offer to delete its
     output.
4. **Results.**
   - Junctions by fate and splice type, and records per file: the log's
     summary, computed from `junctions.tsv` by the code that writes it.
   - Junctions that failed, and the warnings.
   - The file to search: `jcast_uniprot_merged.fasta` with `--uniprot`,
     otherwise as the README's FAQ says.
   - Buttons to open the files and the folder.
   - It opens any earlier run's output folder too, including runs made on
     the command line.

**Not in v0.5.0:**
- **A junction browser** (a gene's junctions and their sequences). A
  junction's records cannot be found by their headers alone: a record
  whose sequence is already in its file is not written again
  (`new_records` in `junctions.tsv`), so a later junction's sequence can
  sit under an earlier junction's header. It needs its own design.
  `junctions.tsv` opens in R or a spreadsheet meanwhile. v0.5.x.
  *2026-11-14: moved into v0.5.0, after a feasibility analysis
  ([REVIEW.md §9](../REVIEW.md#v050-the-gui-and-the-command-line), row
  22).*
- Exon or sequence diagrams; the thresholds in `params.py` (not options
  on the command line either); the UniProt download helper (command line
  first).
- **Standalone app bundles** (PyInstaller, Briefcase). On macOS, an app
  opens without Gatekeeper's warning only if signed and notarized, which
  needs an Apple developer account; and the child process needs a
  bundle-specific way to start itself. `pip install 'jcast[gui]'` first.

## 5. The command line on Typer

`typer_cli.py` is a command line shaped like JCAST's on Typer. The
evidence script ran it on Typer 0.26.0 and 0.29.1, with the same results:

| Command line | Exit status | Result |
|---|---|---|
| `run a b c -s SE` | 0 | `['SE']` |
| `run a b c -s SE MXE` | 2 | `Got unexpected extra argument(s) (MXE)` |
| `run a b c -s SE -s MXE` | 0 | `['SE', 'MXE']` |
| `run a b c -q 0 0.05` | 0 | `(0.0, 0.05)` |
| `run a b c --g_or_ln 1` | 2 | `No such option: --g_or_ln` |
| `a b c` | 2 | `No such command 'a'.` |
| `gui` | 0 | the `gui` command |

- **`-s`.** Click, which Typer is built on, gives an option a fixed number
  of values, so `-s SE MXE` fails. `main()` will rewrite `-s A B` as
  `-s A -s B` before parsing. The values come from a fixed list, so the
  rewrite cannot take a path. `-q 0 0.05` (two values) works as before.
- **Dependencies.** From each release's wheel metadata:

  | Typer | Click | Other dependencies |
  |---|---|---|
  | 0.24.0, 0.25.1 | a dependency, `click>=8.2.1` | rich, shellingham, annotated-doc |
  | 0.26.0, 0.28.0 | vendored | rich, shellingham, annotated-doc |
  | 0.29.0, 0.29.1 | vendored | also pydantic (compiled) |

  The CI job for the lowest versions resolves the direct dependencies
  only to their lowest (`--resolution lowest-direct`). With a floor below
  0.26 it would pair an old Typer with the newest Click. The floor is
  therefore 0.26. Typer also brings rich, markdown-it-py, pygments and
  shellingham.
- **Help text.** Rich, which formats Typer's help, read the hand-written
  `[default: 1, which keeps every junction]` as markup and dropped it;
  the `-r` line showed `[default: 1]`. Typer prints the defaults itself.
  Explanations go outside square brackets.
- **Unchanged:**
  - the options' spellings (`--cache_dir`, `--splice_type`); Typer would
    otherwise derive `--cache-dir`;
  - exit status 2 for a usage error;
  - `run_jcast` takes an `argparse.Namespace`, so the log's
    `Arguments: Namespace(…)` line keeps its form, less `func=` and
    `mask=`.
- **Added:**
  - `-h`, since Click knows only `--help`;
  - `--version`;
  - options grouped in the help (output; filters; canonicals; the run).
- **Removed:** `--mask`, deprecated in 0.4.0 to be removed in this
  release.
- **One `validate()`** holds the checks now in `main()`: the files, the
  rMATS files of `-s`, `-p` needs `--uniprot`, the `-q` range, an output
  folder holding JCAST output. The command line and the GUI both call it,
  so they report the same messages. Typer's own checks that a path exists
  are not used.
- **What else changes:** the error messages Click writes; six test files
  that call `main.main()`; and the help shown in the README and on the
  website (Rich draws boxes even when the output is piped).

## 6. The model preview's cost

The `-m` model fitted from an rMATS folder alone, as the preview will:

| rMATS folder | Junctions | Cutoff | Read the tables | Read counts | Fit | In all |
|---|---|---|---|---|---|---|
| chr15 (bundled) | 1,105 | 26 | 0.0 s | 0.0 s | 0.0 s | 0.1 s |
| whole genome | 599,644 | 40 | 0.7 s | 20.7 s | 5.2 s | 26.7 s |

- **Reading the tables.** The whole-genome files are on a network volume.
  The first read in this session took 13.5 s; the evidence run read them
  from the cache.
- **The read counts** take most of the time. `get_junction_count_array`
  builds a `Junction` for each row (`iloc` per row) to read its `min_jc`.
  Computing the counts from the tables' columns at once should take a
  fraction of that, and a `-m` run would save the same time. Its result
  must equal the current array (PR 19).
- One timing each, on a shared machine.

The preview runs in the background, with a busy indicator, either way.

## 7. Decisions (2026-11-13)

1. **`jcast` is a group of commands, `run` and `gui`.** A command line not
   starting with a command name runs `run`. Rejected: one command, with
   `gui` handled as a special case.
2. **`-s SE MXE` keeps working,** rewritten as `-s SE -s MXE` before
   parsing; the repeated form works too. Rejected: making it a breaking
   change.
3. **`typer>=0.26,<1`,** accepting pydantic from Typer 0.29. Rejected:
   capping below 0.29.
4. **The GUI runs JCAST in a child process.** Rejected: `run_jcast` in a
   `QThread`, without `-p`.
5. **Results in v0.5.0:** the model before the run, and the run's summary
   after it. A junction browser in v0.5.x.
6. **No standalone app bundles in v0.5.0.**
7. **No hyphenated aliases** (`--cache-dir`); the options keep their
   spellings.
8. *(After PR 18, §9.3.)* **Abbreviated long options are not restored**:
   argparse took them, Click does not.

Also agreed, from §4 and §5:
- `-h` and `--version`;
- `--mask` removed;
- one `validate()` for the command line and the GUI;
- the `gui` extra on `PySide6-Essentials` (QtCore, QtGui and QtWidgets,
  all the GUI uses) instead of the whole of PySide6;
- the GUI tested with pytest-qt, without a display, in one CI job.

## 8. Plan

| PR | What | Output |
|---|---|---|
| 18 | The command line on Typer: the `run` and `gui` commands, `validate()`, `-s` as before, `-h`, `--version`, grouped help, `--mask` removed | No change, but for the log's `Arguments` line |
| 19 | What the GUI needs from the pipeline: the model plot on a supplied `Axes` (`model.pdf` written without pyplot); a preview function on the run's model code; the read counts computed from the tables' columns; a hidden progress option; a clean end on the stop signal; the run summary from `junctions.tsv` | No change |
| 13 | The GUI: the four pages of §4, its tests, a CI job, the README and the website | No change |

Each PR is verified like PR 12: the same command lines run with the code
before and after, comparing output files, cutoffs, logs and exit statuses
([command-line checks §5](2026-11-13-command-line-checks.md#5-verification)).

Small items for v0.5.0: PyPI's trusted publishing; deleting the merged
`v0.4.0-hardening` branch.

## 9. Implementation and verification: PR 18 (`5479ee3`)

### 9.1 What PR 18 changes

- **[jcast/main.py](../jcast/main.py):**
  - A Typer app with two commands, `run` and `gui`, and `--version`.
  - `main(argv=None)` passes the command line through `_as_before`
    first, which:
    - puts `run` in front when the first argument is not a command name
      or `-h`, `--help`, `--version`;
    - turns `--gui` into `gui`;
    - shows the help when there are no arguments;
    - gives each splice type after one `-s` its own `-s`.

    `main()` returns when the command succeeds, as before, so the tests
    and the `jcast` script call it unchanged.
  - `validate(args)` holds the checks `main()` made, with their
    messages unchanged. It returns every problem, as (argument, message)
    for the GUI. The command line reports them together as a usage error;
    before, it stopped at the first.
  - `run` builds the `argparse.Namespace` that argparse gave `run_jcast`,
    in its order, less `mask` and `func`.
  - The help: options grouped as Output, Filters, Canonicals and The run;
    `-h` and `--help`; defaults printed by Typer. The `-c` and `-q` texts
    are now the README's, the two typos are gone, and the description is
    no longer misspelled. Typer's pretty tracebacks are off, so an error
    prints Python's usual traceback.
  - `--mask` is removed, with the warning `run_jcast` logged for it.
- [jcast/__main__.py](../jcast/__main__.py) calls `main()`;
  `python -m jcast --gui` goes through `_as_before`.
- `pyproject.toml`: `typer>=0.26,<1`.
- **Tests:**
  - [tests/test_command_line.py](../tests/test_command_line.py): 12 new
    tests and one changed (`--g_or_ln` and `--mask`, now one
    parametrized test). All 13 fail on the code before (`fff9981`).
  - `tests/conftest.py`: a wide console for every test, since Typer
    prints a usage error in a box as wide as the console.
  - The arguments the tests build for `run_jcast` no longer have `mask`.
- **Docs:** the README's and the website's options (the help as printed
  by Typer 0.29.1 at 88 columns), and `jcast gui`.

### 9.2 Verification

`verify_18.py` ran the chr15 test set with the code before (`fff9981`)
and after (`72dfab5`, a snapshot of the working tree whose `jcast/` and
`tests/` are those of `5479ee3`), each frozen in a worktree
([verify_18.txt](cli_checks/verify_18/verify_18.txt)):

| Run (`<inputs>` = rMATS, GTF, genome, `-c --uniprot`) | Output files | Cutoff | Log |
|---|---|---|---|
| default | 14, byte-identical | 1 / 1 | identical, but the `Arguments` line |
| `-m` | 11, byte-identical | 26 / 26 | the same |
| `-s SE MXE -m` | 11, byte-identical | 26 / 26 | the same |
| `-r 30 -q 0 0.05 -s RI A5SS` | 6, byte-identical | 30 / 30 | the same |
| `-p 2` | 14, byte-identical | 1 / 1 | the same |
| after only: `run <inputs> -m`; `<inputs> -s SE -s MXE -m` | identical to `-m` and to `-s SE MXE -m` | | |

- **The `Arguments` line** holds the same values in the same order, less
  `mask=None` and `func=<function run_jcast at 0x…>`. The default `-q`
  prints as `[0.0, 1.0]` (`[0, 1]` before).
- **JCAST's own usage errors** have the same exit status (2) and message:
  - a missing GTF, rMATS folder or rMATS file of `-s`;
  - `-p 2` without `--uniprot`, and `-p 0`;
  - `-q` out of range;
  - an output directory holding JCAST output.
- **Changed:**

  | Command line | Before | After |
  |---|---|---|
  | `jcast --version` | usage error | `jcast 0.4.0` |
  | `jcast` | the options (`usage: __main__.py …` with `python -m jcast`) | the commands |
  | several problems | the first | all of them |
  | `--mask N` | runs, with a deprecation warning | usage error: `No such option: --mask` |
  | `-s SE -s MXE` | MXE only: the second `-s` replaced the first | SE and MXE |
  | `-s XX`, `--g_or_ln 1` | argparse's message | Click's (`Invalid value for '-s' …`, `No such option: --g_or_ln`) |
  | `--overw` (an abbreviation) | `--overwrite` | usage error: `No such option: --overw (Possible options: --overwrite, --verbose)` |

- **Typer 0.26.0, the floor:** the 157 fast tests pass with it in place
  of 0.29.1. Its help names the types `TEXT` and `INTEGER`, where 0.29.1
  writes `<str>` and `<int>`.
- **The whole suite** passes: 162 tests offline, 165 with the three that
  query the UniProt API (153 before PR 18).

### 9.3 Abbreviated long options

Not among the decisions of §7: argparse took unambiguous abbreviations of
long options (`--overw`, `--uni`); Click takes only the full name, and
suggests it. `_as_before` could have expanded them as argparse did.

**Decided (2026-11-13): not restored.** Long options are written in full.

### 9.4 CI's first run (2026-11-14)

The branch's first push ran PR 18's tests on CI for the first time (run
43662548893). The GUI's job passed. The other twelve failed the same 7
tests: Linux and macOS, every Python version, the lowest dependency
versions, and the sdist's tests:
- `test_removed_options` (two);
- `test_every_problem_is_reported`;
- `test_help`;
- `test_cli_rejects_processes` (two);
- `test_cli_refuses_directory_with_output`.

**The cause was in the tests, not jcast.** Typer treats the console as a
terminal when `GITHUB_ACTIONS`, `FORCE_COLOR` or `PY_COLORS` is set. It
reads them once, at import, into `typer.rich_utils.FORCE_TERMINAL`, in
0.26.0 as in 0.29.1. Rich then styles the text, and escape codes split
the words the tests look for: `-p/--processes`, and words of the help.
The tests expect the plain text printed to a console that is not a
terminal. `GITHUB_ACTIONS=true` reproduces all 7 locally.

**Fixed in `tests/conftest.py`.** Every test sets `FORCE_TERMINAL` to
`False`, which Rich checks before any environment variable:
- the 7 pass with each of `GITHUB_ACTIONS`, `FORCE_COLOR`, `PY_COLORS`
  and `TTY_COMPATIBLE` set;
- the whole suite passes with `GITHUB_ACTIONS=true`.

jcast itself is unchanged: on CI its messages stay coloured, as Typer
intends.

## 10. Implementation and verification: PR 19 (`ec1068a`)

What the GUI needs from the pipeline (§8). No output file changes.

### 10.1 What PR 19 changes

- **The filters, shared by the run and the preview.**
  `fates.filtered(min_jc, fdr, cutoff, qvalue)` gives the fate of a
  junction that the read count or FDR filter skips, or `None`.
  `_translate` calls it for each junction, with the comparisons it made
  before.
- **The read counts from the tables' columns.**
  - `junctions.read_count(ijc_s1, sjc_s1, ijc_s2, sjc_s2)` is
    `Junction.min_jc`'s code, which now calls it.
  - `RmatsResults.junction_counts()` gives every junction's splice type,
    read count and FDR from the columns, without making a `Junction` per
    row.
  - `get_junction_count_array` builds the model's array from it, and
    keeps the table in `counts`.
- **[jcast/model.py](../jcast/model.py):**
  - `preview(rmats_dir)`: the `-m` model as a run fits it, on every rMATS
    file in the folder. It returns the junctions (type, read count, FDR),
    ln(count + 1), the parameters and the cutoff.
  - `filter_counts(counts, cutoff, qvalue)`: per splice type, the
    junctions, those the two filters skip, and those kept, by
    `fates.filtered`.
  - `draw_model(ax, …)` draws the model on a caller's `Axes`.
    `plot_general_mixture_model` draws on a `matplotlib.figure.Figure`
    with it, without pyplot, and saves `model.pdf` as before.
- **[jcast/outcomes.py](../jcast/outcomes.py):**
  - `read_counts(junctions.tsv)` counts a run's junctions by fate and
    splice type, as `JunctionTable` does.
  - `summary(counts, splice_types, out_dir)`, which writes the log's
    summary, and `records_per_file(out_dir)` are now module functions, so
    the GUI shows any run's summary from its output.
- **[jcast/main.py](../jcast/main.py):**
  - `run --progress` (hidden from the help): each step of the run as a
    line of JSON on stdout:
    - `{"stage": …}` for `rmats`, `gtf`, `genome`, `uniprot`, `model`,
      `outputs` and `done`;
    - `{"stage": "junctions", "splice_type", "done", "total"}` about every
      1 % of each splice type's junctions, instead of the progress bar.

    `run_jcast(args, progress=None)` takes the function that receives
    them.
  - **SIGTERM, the GUI's cancel,** stops a run from the command line
    cleanly. It raises `Stopped`, a `BaseException`, so the junction it
    interrupts does not record it as its own failure. The worker
    processes are terminated (they reset SIGTERM to its default), the log
    and the console say `Stopped by SIGTERM: the output is incomplete…`,
    `junctions.tsv` is closed, and the exit status is 143. The command
    restores the previous handler when the run ends.
- **Tests:** 25 new, all failing on the code before (`8975d02`):
  - in [tests/test_model.py](../tests/test_model.py):
    - the read counts equal those of a `Junction` per row, in the 14
      bundled rMATS folders;
    - `read_count`'s rules;
    - the preview's cutoff: 26 on chr15, 6 with SE only;
    - `filter_counts` against a run's `junctions.tsv`;
    - drawing on a caller's `Axes`;
  - the summary from `junctions.tsv`, line for line the log's;
  - `--progress`;
  - SIGTERM with 1 and 2 processes.
  - Three stubs of `run_jcast` in the tests take the new `progress`
    argument.

### 10.2 Verification

`verify_19.py` compared the code before (`8975d02`) and after (`dbf0d58`,
a snapshot of the working tree, whose code differs from `ec1068a` only in
the docstring of `read_count`), each frozen in a worktree
([verify_19.txt](gui_scope/verify_19/verify_19.txt)):

| Check | Result |
|---|---|
| chr15, `-c --uniprot` with default options, `-m`, `-r 30 -q 0 0.05`, `-p 2` | Output files byte-identical; logs identical but for timestamps, the output directory, the time taken and the order in which `-s`'s default set prints |
| The read count array of the 14 bundled rMATS folders | Identical values, dtype (int64) and shape |
| The cutoff | Unchanged in each folder: 26 on chr15; 4 to 29 in the others; no fit in `encode_human_stomach` |
| `model.preview`'s cutoff | The run's, in each folder |
| `model.pdf`'s figure, rendered as PNG | Byte-identical in the 12 folders with a fit |
| The whole-genome rMATS folder (599,644 junctions) | Identical read count array; cutoff 40 before, after and in the preview; the counts in 20.6 s before and 1.7 s after (back to back; the machine's load average was about 25) |

- **The verification's first run was wrong, and is not in the evidence.**
  It started its processes from the repository root. `python -m` and
  `-c` put the current directory on `sys.path` before `PYTHONPATH`, so
  both sides imported the working tree's code. `verify_19.py` now runs
  each process in a scratch directory, and records which `jcast` each
  side imported and how it computes the counts.
- **SIGTERM before PR 19** (a `-p 2` run of the code before):
  - the process ended with the signal (exit status −15), and the log
    ended at `Translating junctions in 2 worker processes.`;
  - nothing in the log or in `junctions.tsv` said the run was incomplete;
  - the workers exited too.
- **The whole suite** passes: 187 tests offline (162 before PR 19, and the
  25 new), 190 with the three that query the UniProt API.

### 10.3 Found while verifying

**`-m` on rMATS files without junctions fails.** `encode_mouse_fetalheart`
holds only the files' headers. On it, `general_mixture_model` raises
`ZeroDivisionError` (`numpy.average`, weights summing to zero), before
and after PR 19. A run with `-m` on such a folder therefore fails,
instead of falling back to `-r` as for a degenerate fit, and the preview
fails the same way. Not changed here: it changes what a run does, so it
needs a decision. **Decided 2026-11-13, done in PR 20 (§11).**

## 11. `-m` without a fit falls back to `-r` (PR 20)

### 11.1 Decision (2026-11-13)

**When the `-m` model cannot be fitted, for any reason, the run warns,
uses `-r` (1 by default: every junction kept) and continues.** Before,
only a degenerate fit fell back to `-r`. A fit that raised an exception,
such as `ZeroDivisionError` on rMATS files holding no junctions (§10.3),
stopped the run.

Which inputs fail: of the edge cases tried (no junctions; 1, 2 or 3
junctions; all counts equal; all counts 0; two values; one large count
among zeros), only no junctions raises. The others give a degenerate
fit, which already fell back to `-r`.

### 11.2 What PR 20 changes

- **`model.fit(min_jc_array)`** fits the model for the run and the
  preview alike. It returns `general_mixture_model`'s results and, when
  there is no fit, why:
  - `failed: the rMATS files read hold no junctions`, checked before
    fitting;
  - `failed (<exception>: <message>)` for any exception of the fit;
  - `was degenerate`, as before.
- **The run** logs `Mixture model fit <why>; falling back to -r <N>.` as
  a warning, on the console too. The message for a degenerate fit is
  unchanged.
- **The log line with the cutoff names `-r` when `-m` has no fit:**
  `Translating junctions with a read count of at least N (-r)`. It said
  `(-m)` whenever `-m` was given, so a degenerate fit's run said `(-m)`
  for `-r`'s value.
- **`model.preview`** returns the reason too, as its fifth value, for
  the GUI to show.
- **Tests:**
  - 10 new in [tests/test_model.py](../tests/test_model.py). 8 fail on
    the code before (`b3c0852`); the other 2 pin cases where the log line
    does not change (`-m` fitted, and `-r`).
  - 2 changed: the preview's fifth value.

### 11.3 Verification

`verify_20.py`, with the code before (`b3c0852`), after (`a2e5f6b`, a
snapshot of the working tree whose code is `2f0e54c`'s) and before PR 19 (`8975d02`), each run
from a scratch directory ([verify_20.txt](gui_scope/verify_20/verify_20.txt)):

| Run | Before | After |
|---|---|---|
| chr15, default and `-m` | | Output files byte-identical; logs identical (cutoffs 1 and 26) |
| `encode_human_stomach` (degenerate), `-m -s RI -r 100000` | `… at least 100000 (-m)` | Output files byte-identical; `… at least 100000 (-r)` |
| header-only rMATS files, `-m` | exit 1, `ZeroDivisionError` | exit 0; the warning; cutoff 1 (`-r`); `junctions.tsv` without rows |
| the same, `-m -r 3` | exit 1, `ZeroDivisionError` | exit 0; cutoff 3 (`-r`) |
| the same, without `-m` | exit 0 | exit 0 |

**PR 19 changed a run without saying so.** Before PR 19 (`8975d02`),
header-only rMATS files failed a run even without `-m`: `KeyError: 0`,
from the progress bar's label (`rma.jxn_type[0]` of an empty table).
PR 19's progress code labels the bar with the splice type instead, so
since PR 19 such a run completes with an empty `junctions.tsv`. The
report's "no output change" (§10) holds for runs that worked before, but
this change of behaviour went unnoticed until this verification.

The whole suite passes: 197 tests offline (187 before PR 20, and the 10
new), 200 with the three that query the UniProt API.

## 12. Implementation and verification: PR 13, the GUI (`12e0576`)

### 12.1 What PR 13 changes

**The GUI runs JCAST.** [jcast/gui/](../jcast/gui/) has 1,220 lines
instead of the 157 of the placeholder pages it replaces. Four pages, as
tabs, as decided in §4:

1. **Setup** ([setup_page.py](../jcast/gui/setup_page.py)). The paths,
   the splice types, `-q`, `-c`, `--uniprot`, `-p`, `-v` and
   `--cache_dir`.
   - The page builds the `argparse.Namespace` that `run_jcast` takes. It
     checks it with `main.validate`, the command line's checks, and shows
     each message under its field. Each field's tooltip is its help from
     `jcast run -h` (`command.option_help`).
   - Splice types whose rMATS file is not in the folder cannot be
     selected.
2. **Model** ([model_page.py](../jcast/gui/model_page.py)).
   - `model.preview` runs in a background thread when the page opens:
     the run's own fit, on every rMATS file in the folder.
   - `model.draw_model` draws the result on the page. Without a fit, the
     page says why and draws the histogram only.
   - The cutoff: the model's (`-m`), or another (`-r`) set in a number
     field or by clicking or dragging on the plot.
   - A table shows, for the selected splice types, the junctions, those
     below the cutoff, those outside `-q`, and those kept
     (`model.filter_counts`).
3. **Run** ([run_page.py](../jcast/gui/run_page.py)).
   - The page shows the command line (`jcast run …`, with a Copy
     button). It runs that command in a child process (`QProcess`):
     `python -P -m jcast run … --progress`, with the GUI's Python
     (`command.process`).
   - It shows the stage of the run, a progress bar per splice type (the
     JSON lines of `--progress`), and the console's warnings and errors
     as they come.
   - Stop sends SIGTERM; a run that has not ended 15 s later is killed.
     A stopped or failed run is marked incomplete, and Delete removes the
     files it wrote.
   - An output folder holding JCAST output is replaced only after a
     question (`--overwrite`).
4. **Results** ([results_page.py](../jcast/gui/results_page.py)).
   - The summary of any output folder: this run's, an earlier one's, or
     one made on the command line. It is read in the background.
   - Junctions by fate and splice type (`outcomes.read_counts` and
     `outcomes.fate_table`, the code of the log's summary), records per
     file, the junctions that failed, and the log's warnings and errors.
   - The version, cutoff and time from the log, and whether the run
     finished (its log has the summary).
   - The file to search (the README's FAQ), and buttons that open the
     folder, the log, `junctions.tsv` and `model.pdf`.

**`-P` and `PYTHONPATH`: the child runs the GUI's `jcast`.**
`python -m` puts the working directory first on `sys.path`. A child
started in a folder holding another `jcast` package (a checkout; a test's
temporary folder) would import that one, or fail if `jcast` is not
installed. `-P` (Python 3.11) leaves the working directory out. If the
GUI's `jcast` is not in a site-packages folder (a source checkout, or an
editable install), its folder is put first on `PYTHONPATH`. A test runs
the child in a folder holding a `jcast` package that cannot be imported.

**What the GUI chose that §4 did not settle.** The maintainer agreed to
each of these on 2026-11-15:
- The cutoff is `-r 1` until the model's is chosen, as on the command
  line. The model's cutoff is not chosen by itself once fitted.
- Another rMATS folder resets the cutoff to `-r`. The model's cutoff of
  one folder is not kept for another.
- The model is fitted when the Model page opens, not when the folder is
  chosen, so typing a path does not start fits.
- `-p` can only be set with `--uniprot`, and is 1 otherwise.
- Fields left empty show "required" once a run is started, not before.
- Stop asks first. Quitting during a run asks, then stops it as Stop
  does.
- Delete removes the files a run writes (`main.existing_outputs`), not
  the folder. It is offered only if the run wrote its log: a run that
  fails before starting leaves an earlier run's output in place.
- Settings are not remembered between sessions.

**Pipeline changes, for the pages; no output change:**
- `outcomes.fate_table`: the summary's table, as numbers. The log's
  summary is formatted from it.
- `model.draw_model` draws without a fit (the histogram and the `-r`
  cutoff) and returns the cutoff's line and label. `model.move_cutoff`
  moves them, with the label's text written in one place.

**Packaging and CI:**
- The `gui` extra installs `PySide6-Essentials>=6.7` (QtCore, QtGui and
  QtWidgets) instead of PySide6. 6.6.0 and 6.6.3.1 print NumPy's warning
  for a module compiled with NumPy 1.x on import; 6.7.0 does not.
- The GUI's tests pass with `PySide6-Essentials` 6.7.0 and pytest-qt
  4.4.0, the floors, and with 6.12.0 and 4.5.0.
- A new `test-gui` extra (`jcast[gui,test]` and `pytest-qt>=4.4`).
  pytest-qt is not in `test`: without a Qt binding it stops every test
  run with a usage error.
- In a fresh virtual environment, `pip install '.[test-gui]'` installs
  PySide6-Essentials 6.12.0 and pytest-qt 4.5.0, and the GUI's tests pass
  against the installed wheel.
- **CI:** a `gui` job on macos-latest, Python 3.13, runs
  `tests/test_gui.py` without a display (`QT_QPA_PLATFORM=offscreen`).
  The other jobs skip those tests: they do not install PySide6.
- The GUI is not tested on Linux. There, Qt's offscreen platform needs
  system libraries, which could not be checked without pushing.

**Tests:**
- 15 in [tests/test_gui.py](../tests/test_gui.py), 3 of them real runs
  (marked slow): a run whose output equals its command line's, a stop,
  and quitting during a run.
- 11 more in `test_command_line.py`, `test_model.py` and
  `test_outcomes.py`:
  - the GUI's command line parses back into the same arguments, for 6
    settings;
  - the child runs the GUI's `jcast`;
  - the tooltips;
  - `jcast gui` without PySide6;
  - `draw_model` without a fit;
  - `fate_table`.

**Documentation:**
- README: a section on the GUI, the installation, the `test-gui` extra,
  and a known issue (no junction browser yet).
- The website's installation note, and `jcast -h`, whose `gui` line no
  longer says "a preview".

### 12.2 Verification

`verify_13.py` compares the code before (`d9494a3`) with the code after:
`905f359`, a snapshot of the working tree whose code is `12e0576`'s. Each
process runs from a scratch directory, and the GUI offscreen
([verify_13.txt](gui_scope/verify_13/verify_13.txt)).

**The pipeline: no output change.**
- chr15 runs with `-c --uniprot` and with `-m --uniprot`: the output files
  (14 and 11) are byte-identical, and the logs identical but for
  timestamps, the output directory and the time taken.
- `model.pdf` of the 12 bundled folders with a fit, drawn as PNG: the
  same images.

**A GUI run is its command line's run.** Each GUI run, then the command
the Run page showed, run with another `-o`:

| Settings (Setup and Model pages) | GUI | Command line | Output |
|---|---|---|---|
| all splice types, `--uniprot` | 9.2 s | 9.2 s | 14 files byte-identical; logs identical |
| `-m` chosen on the Model page, `--uniprot -c` | 4.1 s | 4.2 s | 11 files byte-identical; logs identical |
| `-r 10` by a click on the plot, `-s RI A5SS`, `-q 0 0.05`, `--uniprot` | 2.0 s | 2.0 s | 6 files byte-identical; logs identical |
| `--uniprot -c -p 2 -v` | 6.2 s | 6.3 s | 14 files byte-identical; logs identical |

- The child process imported the snapshot's `jcast`: its `PYTHONPATH`
  starts with that checkout.
- Each run showed its stages, and every splice type's progress to the
  end.
- The Results page read each cutoff back from the log: `1 read (-r)`,
  `26 reads (-m)`, `10 reads (-r)`.

**Stop.** A run of every chr15 junction (`--uniprot -p 2`) was stopped
after 40 of the 90 MXE junctions:
- exit status 143;
- both worker processes ended;
- the log ends with the stop warning;
- Delete removed the 7 files written, and nothing else was in the folder.

**The Model page, on each of the 14 bundled folders:**
- The 12 with a fit show the cutoff that `model.preview` gave in §10.2.
- `encode_human_stomach` shows "the fit was degenerate", and the
  header-only folder "the fit failed: the rMATS files read hold no
  junctions". `-m` cannot be chosen for either.
- Each takes 0.02 to 1.3 s.

**The whole genome** (599,644 junctions):

| Model page | Time |
|---|---|
| The fit, in the background | 9.1 s |
| The plot | 0.09 s |
| A cutoff moved and drawn again | 0.05 s at most |
| The table | 0.17 s |

The Results page read the whole-genome output folder of 2026-11-13 (`-m`,
`--uniprot`, 101,734 records) in 2.4 s, in the background.

**Screenshots** of the four pages after a `-m` run:
[Setup](gui_scope/verify_13/setup.png),
[Model](gui_scope/verify_13/model.png),
[Run](gui_scope/verify_13/run.png),
[Results](gui_scope/verify_13/results.png).
- The paths are relative to the scratch directory.
- The Run page's console box shows gtfparse's
  `INFO:root:Extracted GTF attributes …`. The command line prints it too,
  whenever it parses a GTF that is not yet in the cache.

**Tests:** 226 pass (223 offline): the 197 before and the 26 new.
- The new tests test new code (`jcast.gui`'s modules, `fate_table`,
  `move_cutoff`, `draw_model` without a fit), so they fail on the code
  before.
- The GUI's tests also pass with the floors (PySide6-Essentials 6.7.0,
  pytest-qt 4.4.0), and after `pip install '.[test-gui]'` into a fresh
  virtual environment.

### 12.3 Not done

- **The junction browser** (§4): v0.5.x.
- **Remembering the settings** between sessions.
- **The GUI on Linux and Windows**, which CI does not test. On Windows,
  `QProcess.terminate` does not stop a console program; Stop kills the
  run after 15 s, so its log does not say that it stopped.
- **Screenshots in the README or on the website.**

## 13. The GUI after its first run on a Mac (PR 21)

### 13.1 Requested (2026-11-14)

The maintainer ran the GUI on a Mac and asked for three changes before
v0.5.0 is pushed:

1. The settings' help shown on hover, perhaps in a fixed place at the
   bottom of the window.
2. More than a filling bar on the Run page, such as counts or the gene
   being translated, if that works with `-p`.
3. The JCAST logo shipped with the GUI.

The bars did show counts in the screenshots of §12.2, but those were
drawn offscreen in Qt's Fusion style. **The macOS style draws no text in
a progress bar**, so on a Mac the bar was all there was.

### 13.2 What PR 21 changes

**The help line** ([main_window.py](../jcast/gui/main_window.py)).
- A panel under the pages shows the help of what the mouse is over: a
  setting, a button, the plot or a table. Its text is the widget's status
  tip, which Qt sends to the window when the mouse enters or leaves the
  widget. The window shows it in the panel, not in the status bar, which
  has room for one line only.
- For a setting, the help line is its spellings and its help in
  `jcast run -h` (`command.option_tip`), for example
  **-p, --processes** — number of processes to translate junctions in; ….
- A splice type that cannot be selected says why: its rMATS file is not
  in the folder.
- The tooltips of PR 13 are gone. They were on the fields' labels only,
  not on the fields.

**Progress** ([run_page.py](../jcast/gui/run_page.py),
`main._each`).
- Each splice type's count is now beside its bar, e.g.
  "203 of 720 (28%)", and a line under the bars names the last junction
  done: "Last junction done: SE 6211, gene BNIP2 (ENSG00000140299)".
- `jcast run --progress` adds the junction (`junction`, its rMATS ID;
  `gene`; `gene_id`) to each report after a junction is done.
- **With `-p`**, the run writes the workers' results in file order, so
  the last junction written is well defined: the line names it. The
  workers may already be further on.
- Reports come about every 1 % of a splice type's junctions, as before,
  and now also at least every 0.5 s (`PROGRESS_SECONDS`). In the
  maintainer's whole-genome run of 2026-11-13 (`-m -c --uniprot -p 8`;
  8 min 20 s translating 599,644 junctions), 1 % of the SE junctions was
  2,602 junctions, about 2 s apart; with `-p 1`, several times longer.
- The console's progress bar (tqdm) and the output do not change.

**The icon** ([jcast/gui/icon.png](../jcast/gui/icon.png), 512 × 512,
121 kB).
- It is the website's fish logo, its transparent margins trimmed and
  centred on a square
  ([make_icon.py](gui_scope/make_icon.py) makes it again).
- The GUI sets it as the application's icon (on macOS, the Dock's) and
  the window's.
- It ships in the package as package data (`[tool.setuptools.package-data]`).
- Of the two logos offered, the fish, since it is close to square once
  trimmed. The wordmark (`jcast.png`) is twice as wide as it is high.

**Tests.**
- New:
  - the help line (`test_the_help_line`);
  - the icon (`test_the_icon`);
  - the counts and the last junction after a run (in `test_a_run`);
  - the junction named in each progress report, with `-p 1` and `-p 2`
    (`test_progress`);
  - a report at least every `PROGRESS_SECONDS`
    (`test_progress_at_least_every_progress_seconds`).
- Changed: `test_option_help`.
- All of them fail on the code before (`fe65069`). `tests/test_gui.py`
  fails there as it is imported: the help line and the icon do not exist.

### 13.3 Verification

`verify_21.py` compares the code before (`fe65069`) with the code after:
`89bd730`, a snapshot of the working tree whose code is `4cc72c2`'s. Each
process runs from a scratch directory
([verify_21.txt](gui_scope/verify_21/verify_21.txt)).

**No output change.** chr15 runs (`-c --uniprot`), with `-p 1` and with
`-p 2`, each without and with `--progress`: the output files (14) are
byte-identical, and the logs identical but for timestamps, the output
directory and the time taken.

**The progress lines,** each of the four runs with `--progress`:
- 499 lines, 493 of them about junctions, before and after.
- After, 488 name the last junction done: all but each splice type's
  first report, at 0 junctions done. Each of the 488 names the junction
  at that place in its splice type's rows of `junctions.tsv`, with `-p 1`
  and with `-p 2`.
- The longest time between two lines: 0.28 and 0.36 s before (`-p 1`,
  `-p 2`), 0.31 and 0.42 s after. On chr15, 1 % of a splice type comes
  sooner than 0.5 s, so the new rule adds no line there. A test sets the
  interval to 0 s and gets a line for every junction.

**The packages** built from the code after hold `jcast/gui/icon.png`
(121,080 bytes): the wheel and the sdist.

**Screenshots,** taken offscreen, in Qt's Fusion style:
- [the Setup page](gui_scope/verify_21/setup_help.png), the help line
  showing `--uniprot`;
- [the Run page during a run](gui_scope/verify_21/run_progress.png):
  "210 of 720 (29%)" beside the SE bar, and "Last junction done: SE 6230,
  gene ADAM10 (ENSG00000137845)".

On a Mac the bars have no text of their own, so the counts beside them
are what shows; neither screenshot shows the macOS style.

**Tests:** 230 pass (227 offline).

## 14. The GUI keeps its settings (PR 26)

### 14.1 Decided (2026-11-18)

The maintainer moved "remembering the settings between sessions" from
v0.5.x into v0.5.0 (read-count audit §13). The handoff's brief: the
Setup page's paths and options, restored at start, with a way back to
the defaults; the tests on a settings file of their own, never the
user's.

### 14.2 What PR 26 changes

[setup_page.py](../jcast/gui/setup_page.py),
[main_window.py](../jcast/gui/main_window.py).
- **What is kept:** every field of the Setup page. These are the six
  paths (rMATS folder, GTF, genome, output folder, cache folder,
  `--uniprot`), the splice types checked, `-q`, `-c`, `-p` and `-v`.
  `SetupPage.values()` gives them, and `DEFAULTS` names them, as a new
  page has them.
- **When:** each time the page checks a change (300 ms after the typing
  stops), and when the window closes, so that a change not yet checked
  is kept too.
- **Where:** Qt's `QSettings`, under the group `setup`: on macOS
  `~/Library/Preferences/com.jcast.JCAST.plist`, on Linux
  `~/.config/JCAST/JCAST.conf`, on Windows the registry
  (`main_window.user_settings`). A value that cannot be read is replaced
  by its default.
- **At start,** the window sets the page from the kept settings and
  checks them at once. A path that has gone since shows its problem
  under its field, as a typed one would. Empty fields are not marked
  "required" until a run is started, as before.
- **Restore the defaults,** a button at the bottom of the page, asks
  first, then clears every path and sets every option to its default;
  the fields marked "required" are no longer marked.
- **Not kept:** the read count cutoff. It is chosen for each rMATS folder
  on the Model page, which starts at `-r 1` (PR 13's choice, kept on
  2026-11-15: another folder resets the cutoff). The Results page's
  folder is not kept either.

**Tests** (`tests/test_gui.py`, 5 new). Every GUI test now has a settings
file of its own (the fixture `settings_file`, replacing
`main_window.user_settings`). The new tests check:
- a first session has the defaults;
- the settings come back as each change is checked, and as the window
  closes; the cutoff does not;
- a kept folder that has gone shows its problem;
- Restore the defaults, answered yes and no.

### 14.3 Verification

`verify_26.py`
([verify_26.txt](gui_scope/verify_26/verify_26.txt)):
- **No run changes.** The PR changes no file outside `jcast/gui/`, the
  tests and the documents.
- **Qt's native format.** On a Mac the settings are a property list, not
  the INI file the tests use. In a `.plist` of the scratch directory, never
  the user's (whose file was not created):
  - a first session has the defaults;
  - a second has the settings the first filled in, with no problem, and
    `-r 1`;
  - after Restore the defaults, a third has the defaults again;
  - the file keeps each value with its type: a boolean, an integer, a
    float, a string.
- **Screenshot:** [the second session's Setup page](gui_scope/verify_26/setup_kept.png),
  with the help line over Restore the defaults.

**Tests:** 255 pass (252 offline).
