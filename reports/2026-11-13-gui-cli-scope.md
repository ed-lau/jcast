# GUI and command-line scope for v0.5.0 (PRs 18, 19, 13)

| | |
|---|---|
| **Date** | 2026-11-13 |
| **Status** | Decisions taken 2026-11-13 (§7). To be implemented in PRs 18, 19 and 13, in that order (§8). |
| **Question** | What should the v0.5.0 GUI do: start a run only, or also show results, and which? How is it launched? What does moving the command line to Typer change for users? |
| **Reproduce** | `python reports/gui_scope/scope_evidence.py --work <scratch> [--python <interpreter for the venvs>] [--rmats <whole-genome rMATS folder>]` (≈ 2 min; needs PyPI) |
| **Evidence** (code at `f37faf9`) | [gui_scope/scope_evidence/scope_evidence.txt](gui_scope/scope_evidence/scope_evidence.txt) |

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
