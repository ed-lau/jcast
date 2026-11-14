# Handoff — v0.5.0 under way (2026-11-14)

For the next working session. It covers where things stand, how the work
has been run, and what to do next.

## State

- **Released:** v0.4.0, on 2026-11-13 (PR #29, merge commit `e4d9117`,
  tag `v0.4.0`, GitHub release, PyPI 0.4.0 built from `258d8b3`). The
  website was updated for it in PR #30 (merged, `f37faf9`).
- **Branch:** `v0.5.0-gui`, from `main` at `f37faf9`. **Not pushed.** It
  collects the v0.5.0 PRs, as `v0.4.0-hardening` did for v0.4.0:
  - `fff9981`: the [GUI and command-line scope](reports/2026-11-13-gui-cli-scope.md),
    with the maintainer's decisions (§7) and its evidence;
  - `5479ee3`: PR 18, the command line on Typer; `7d36121` records its
    hash; `8975d02` the decision not to restore abbreviated long options;
  - `ec1068a`: PR 19, what the GUI needs from the pipeline; `b3c0852`
    records its hash;
  - `2f0e54c`: PR 20, `-m` without a fit falls back to `-r`, and the commit
    recording its hash;
  - `12e0576`: PR 13, the GUI, and the commit recording its hash;
  - `4cc72c2`: PR 21, the GUI after the maintainer's first run on a Mac
    (the help line, progress counts, the icon), and the commit recording
    its hash.
- **Tests:** 230, all passing; 227 offline (`python -m pytest -q -m 'not network'`,
  about 2.5 minutes). Run only the fast ones with `-m 'not slow'`. The
  GUI's 17 (`tests/test_gui.py`) need PySide6 and pytest-qt
  (`pip install -e '.[test-gui]'`; `.conda` has both) and are skipped
  without them.
- **CI:** [.github/workflows/tests.yml](.github/workflows/tests.yml): Python
  3.11 to 3.15 on Linux and macOS, the lowest dependency versions, the
  GUI's tests on macOS (`gui`, new in PR 13), and the installed wheel. It
  has not run on this branch yet (not pushed). The lowest-versions job
  will be the first to run Typer 0.26.0 on Linux, and the `gui` job is
  new.
- **Golden tests:** `tests/test_pipeline.py::test_chr15_golden` (GTF-only
  canonicals) and `tests/test_uniprot.py::test_chr15_uniprot`
  (`--uniprot`). Unchanged by PR 18.
- **Done:** roadmap PRs 0–12 and 14–17 (v0.4.0), and 18, 19, 20, 13 and 21
  ([REVIEW.md §9](REVIEW.md#v050-the-gui-and-the-command-line)): the
  whole of the v0.5.0 plan.

## The session of 2026-11-14

- **PR 13, the GUI** (report §12, `12e0576`). Four pages in
  [jcast/gui/](jcast/gui/), 1,220 lines:
  - Setup checks with `main.validate`;
  - Model fits `model.preview` in a background thread (`tasks.Background`);
  - Run runs `python -P -m jcast run … --progress` in a `QProcess`
    (`command.process`);
  - Results reads any output folder (`results_page.read_output`).

  [README](README.md#the-graphical-interface) describes them; screenshots
  are in [reports/gui_scope/verify_13/](reports/gui_scope/verify_13/).
- **Choices the maintainer has not seen yet** (report §12.1, "What the GUI
  chose that §4 did not settle"):
  - `-r 1` until the model's cutoff is chosen;
  - another folder resets the cutoff to `-r`;
  - the fit starts when the Model page opens;
  - Delete removes only a run's files, and only if the run wrote its log;
  - settings are not remembered.
- **`-P` and `PYTHONPATH`.** The child imports the GUI's own `jcast`
  wherever it runs: `python -m` imports from the working directory first.
- **Packaging.** The `gui` extra is `PySide6-Essentials>=6.7`; 6.6 prints
  NumPy 2's warning about a module compiled with NumPy 1.x. The new
  `test-gui` extra adds `pytest-qt>=4.4`, which cannot go in `test`: with
  no Qt binding it stops every test run.
- **PR 21** (report §13, `4cc72c2`). The maintainer ran the GUI on a Mac
  and asked for three changes:
  - a help line at the bottom of the window, from status tips: for a
    setting, its help in `jcast run -h`;
  - counts beside the progress bars, since the macOS style draws no text
    in a bar, and the last junction done, which `--progress` now names
    (with `-p`, the last one written). It also reports at least every
    0.5 s;
  - the fish logo as the icon (`jcast/gui/icon.png`, package data, made by
    `reports/gui_scope/make_icon.py`).
- **Verified** (`reports/gui_scope/verify_13.py`):
  - the pipeline changes leave the chr15 output and `model.pdf`
    identical;
  - four GUI runs give the same output as their command line;
  - a stopped run exits 143 and leaves no worker; Delete removes its
    files;
  - on the whole genome, the Model page fits in about 10 s, in the
    background.

## The session of 2026-11-13 (evening)

- **Scope of v0.5.0** ([report](reports/2026-11-13-gui-cli-scope.md), `fff9981`).
  Decided:
  - `jcast` is a group of commands, `run` and `gui`; a command line without
    a command name runs `run`. `-s SE MXE` keeps working.
  - `typer>=0.26,<1`; options keep their spellings.
  - The GUI runs `python -m jcast run …` in a child process (`QProcess`),
    so `-p` works and cancelling stops the process.
  - It shows the `-m` model before the run (the cutoff set there) and the
    run's summary after (for any output folder). A junction browser in
    v0.5.x; no standalone bundles in v0.5.0.
  - The model preview takes 26.7 s on the whole genome, 20.7 s of it in
    `get_junction_count_array` (a `Junction` per row).
- **PR 18, the command line on Typer** (`5479ee3`; report §9):
  - `validate(args)` returns every problem as (argument, message), for the
    GUI; `_as_before` maps the old forms (no command, `--gui`,
    `-s SE MXE`) onto the new command line.
  - `--mask` removed; `--version` added; the help grouped and its typos
    fixed; README and website updated, the site rendered again.
  - Verified (`reports/cli_checks/verify_18.py`): no output change; the
    log's `Arguments` line loses `mask` and `func`.
  - argparse took abbreviated long options (`--overw`); Click does not.
    The maintainer decided not to restore them (report §9.3).
- **PR 19, what the GUI needs from the pipeline** (report §10), no output
  change:
  - `model.preview(rmats_dir)`, `model.filter_counts`, `model.draw_model`
    (on a caller's `Axes`); `fates.filtered`, the filters shared by the
    run and the preview;
  - the read counts from the tables' columns (`junctions.read_count`,
    `RmatsResults.junction_counts`): 20.6 s → 1.7 s on the whole genome;
  - `outcomes.read_counts` and `outcomes.summary`: any output folder's
    summary from its `junctions.tsv`;
  - `jcast run --progress` (hidden): JSON lines on stdout;
  - SIGTERM: `Stopped`, the workers terminated, the log line, exit 143.
  - Verified with `reports/gui_scope/verify_19.py`. Its first run was
    wrong: started from the repository root, both sides imported the
    working tree (see Useful facts).
- **PR 20** (report §11; the maintainer's decision): when `-m` has no fit,
  for any reason, the run warns and falls back to `-r` (`model.fit`).
  rMATS files without junctions raised `ZeroDivisionError`. The cutoff's
  log line names `-r` then. Its verification found that PR 19 had also
  made runs on header-only rMATS files complete (they failed with
  `KeyError: 0`), unnoticed at the time.

## How the work is run (agreed with the maintainer)

1. **Audit first.**
   - Write a dated report in `reports/`, with a script that regenerates
     its evidence.
   - Measure on the chr15 test set, and on the user's whole genome when
     the effect is genome-wide. A prototype behind an audit-only switch
     (kept as a `.diff` in the evidence folder) can measure an option
     before the decision.
   - End the report with the decisions requested.
2. **The maintainer decides.** Record the decisions in the report (a
   decisions section), in the reports/README.md decisions log, and in the
   REVIEW.md roadmap. Commit the audit with them.
3. **Implement, then verify against the audit.**
   - Use a verification script, comparing the commit before with the
     working tree.
   - Freeze the code under test in a `git worktree` while it runs, so
     later edits don't leak in.
   - Add regression tests that fail on the old code where behaviour
     changes.
   - Update the golden values and their history when output changes.
4. **Document** in the report's implementation section, CHANGELOG, README
   (user-facing), REVIEW.md roadmap and reports/README.md.
5. **Commit:**
   - Commit only with the maintainer's approval. Decisions on a PR count
     as approval to implement and commit it.
   - End messages with
     `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
   - Follow each PR with a small commit recording its hash in the report
     and index.
   - Push only when asked.

## Next priorities

1. **Push `v0.5.0-gui` when the maintainer asks**, and check CI: the
   lowest-versions job (Typer 0.26.0 on Linux) and the `gui` job (macOS)
   run for the first time.
2. **Release v0.5.0:** `__version__` (still 0.4.0), the CHANGELOG's date,
   the website (`_quarto.yml`'s title, `change.qmd`), and, before or
   with it:
   - PyPI's trusted publishing: a workflow that uploads when a GitHub
     release is published;
   - delete the merged `v0.4.0-hardening` branch: the commits cited in
     `reports/` are in `main`'s history.
3. **GUI follow-ups:**
   - the junction browser (v0.5.x; report §4);
   - remembering the settings between sessions;
   - the GUI's tests on Linux: Qt's offscreen platform needs system
     libraries there, not tried;
   - screenshots in the README or on the website (they are in
     `reports/gui_scope/verify_13/`).
4. **REVIEW.md §4 LOW items still open:** the terminal-stop tolerance in
   `make_pep`, `chr` never added to chromosome names, the header-only
   fetal-heart test data, the date-stamped comments.
5. **Planned follow-ups** (CHANGELOG `[Unreleased]`):
   - an option to add T2/T4 to `jcast_uniprot_merged.fasta`;
   - a UniProt download helper that names the file by its release
     (`X-UniProt-Release` header), as a third command;
   - joining by genomic coordinates was evaluated and deferred
     ([UniProt evaluation §4](reports/2026-11-11-uniprot-fasta-evaluation.md#4-joining-by-coordinates-question-1)).
     It would also be the better join for readthrough genes' events
     ([`Ensembl_canonical` audit §4](reports/2026-11-12-ensembl-canonical-audit.md#4-readthrough-genes-and-rmats-events)).
   - Possibly: prefer the gene with an HGNC id when a gene name is shared
     ([mapping audit §8](reports/2026-11-12-uniprot-mapping-audit.md#8-addendum-the-unique-name-rule-on-the-whole-genome-pr-8)).
     Not decided.
   - Possibly: print the end-of-run summary on the console too. Not
     decided.

## Useful facts

- **Verification scripts run each process outside the repository.**
  `python -m jcast` and `python -c` put the current directory on
  `sys.path` before `PYTHONPATH`: run from the repository root, the code
  "before" is the working tree's. Set `cwd` to a scratch directory, and
  record `jcast.__file__` in the evidence.
- **The command line (PR 18):**
  - Typer prints usage errors in a box as wide as the console (80 columns
    when not a terminal). `tests/conftest.py` sets `COLUMNS=1000` for
    every test; set it for any script that reads the messages.
  - The help in the README and on the website is
    `COLUMNS=88 NO_COLOR=1 python -m jcast run -h`, trailing spaces
    stripped, as Typer 0.29.1 prints it. Render it again when an option
    changes.
  - Test another Typer version without touching `.conda`:
    `.conda/bin/python -m venv --system-site-packages <dir>`, then
    `<dir>/bin/python -m pip install typer==<version>`.
  - `uv` is not on the PATH of Claude's shell.
- **`quarto render`** copies `jcast_website/images` and `styles.css` into
  `docs/` with mode 755; `chmod 644 docs/images/*.png docs/styles.css`
  before committing, or git records a mode change.
- `run_jcast` creates the output directory and refuses one holding JCAST
  output unless `args.overwrite`. It closes its log handlers when it ends,
  and leaves the `jcast` logger not propagating: a test that reads its
  records with `caplog` must set `propagate` itself (pytest before 9.1
  does not capture them).
- **Tests:**
  - `tests/conftest.py` points the caches at a session temp directory (the
    `cache_dir` fixture).
  - Tests that call `run_jcast` must pass `cache_dir=str(cache_dir)` in
    `args`, or they write to the user's cache.
  - `tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf.cached`
    (git-ignored, from before PR 7) is no longer read and can be deleted.
    `MANIFEST.in` keeps it out of the sdist.
- **Packaging checks** (`reports/packaging/check_install.py`) need `uv` and
  an interpreter per Python version. `uv python install` gives them; set
  `UV_PYTHON_INSTALL_DIR` and `UV_CACHE_DIR` to a scratch directory.
  **Name the target of every `uv pip` command with `--python`:**
  `--system` follows `CONDA_PREFIX` into the conda base environment.
- **The maintainer's environments:**
  - `.conda` is arm64 since 2026-11-13.
  - `~/miniconda3` is the Intel build (osx-64). Environments it creates,
    including those made by VS Code's "Create Environment", are x86_64
    unless built with `CONDA_SUBDIR=osx-arm64` and pinned with
    `conda config --env --set subdir osx-arm64`.
- **The website:**
  - `cd jcast_website && quarto render` writes the site into `docs/`.
    GitHub Pages serves `main:/docs` at https://ed-lau.github.io/jcast/.
  - Quarto warns that `docs/` is outside the project and that it won't
    clean `docs/site_libs`; harmless.
  - `about.qmd` and `faqs.qmd` repeat the README's sections (output files,
    `junctions.tsv`, header format, `--uniprot`, `-p`, known issues, FAQ),
    and `about.qmd` shows `jcast --help` verbatim. Update them with the
    README.
  - Pandoc prints `\|` inside code in a table literally. Write such a cell
    as `<code>…&#124;…</code>`.
- **Building release files:**
  - Build from the commit, not the working tree: `git archive <commit>`,
    then `uv build`.
  - `twine check --strict` them. setuptools ≥ 77 writes Metadata-Version
    2.4: twine ≥ 6.1 reads it, older twine only with pkginfo ≥ 1.12.
  - PyPI accepts each version once.
- **Verification scripts:**
  - They take input paths as arguments and write only to the scratch
    directory.
  - For whole-genome runs, pass `--cache_dir`, or set `HOME`, to keep the
    rebuilt GTF cache out of the user's cache.
  - `junctions.tsv` makes run comparisons cheap: fates, canonicals and
    records per junction, without parsing the log.
  - `reports/runs_audit/verify_6.py` parses the log format from before
    PR 16; use `reports/cli_checks/verify_12.py` or
    `reports/canonical_audit/verify_17.py` as a model instead.
- **The GUI:**
  - Offscreen (`QT_QPA_PLATFORM=offscreen`), Qt draws the Fusion style;
    the macOS style draws nothing there, and on a Mac draws no text in a
    progress bar. Screenshots taken offscreen are not what a Mac shows.
  - `widget.grab().save(<png>)` takes a screenshot without a display.
  - A help line comes from `setStatusTip`; `MainWindow.event` shows it.
    Rich text: escape what goes in it (`command.option_tip`).
  - Test another PySide6 without touching `.conda`:
    `.conda/bin/python -m venv --system-site-packages <dir>`, then
    `<dir>/bin/python -m pip install PySide6-Essentials==<version>`.
- **Timings:** reading the GENCODE v50 GTF cache takes about 8 s, the
  genome about 35 s. In a whole-genome `-m -c --uniprot -p 16` run, the
  junction loop takes about 3 min, and the run about 4 min.
- **Line endings:** `uniprot_canonicals.tsv` is written with CRLF (the csv
  module's default); the other tables use LF.
