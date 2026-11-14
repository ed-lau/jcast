# Handoff — v0.5.0 under way (2026-11-13)

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
  - `5479ee3`: PR 18, the command line on Typer;
  - the commit recording that hash.
- **Tests:** 165, all passing; 162 offline (`python -m pytest -q -m 'not network'`,
  about 1.5 minutes). Run only the fast ones with `-m 'not slow'`.
- **CI:** [.github/workflows/tests.yml](.github/workflows/tests.yml): Python
  3.11 to 3.15 on Linux and macOS, the lowest dependency versions, and the
  installed wheel. It has not run on this branch yet (not pushed); the
  lowest-versions job will be the first to run Typer 0.26.0 on Linux.
- **Golden tests:** `tests/test_pipeline.py::test_chr15_golden` (GTF-only
  canonicals) and `tests/test_uniprot.py::test_chr15_uniprot`
  (`--uniprot`). Unchanged by PR 18.
- **Done:** roadmap PRs 0–12 and 14–17 (v0.4.0), and 18
  ([REVIEW.md §9](REVIEW.md#v050-the-gui-and-the-command-line)).

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
  - **Open, for the maintainer** (report §9.3): argparse took abbreviated
    long options (`--overw`); Click does not. Not restored.

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

1. **The maintainer's answer on abbreviated long options** (report §9.3).
2. **PR 19, what the GUI needs from the pipeline** (report §8), no output
   change:
   - the model plot drawn on a supplied `Axes`; `model.pdf` written with
     `matplotlib.figure.Figure`, without pyplot;
   - a model-preview function on the run's code (with `-m`, every rMATS
     file in the folder);
   - the read counts from the tables' columns, equal to
     `get_junction_count_array`'s (20.7 s on the whole genome now);
   - a hidden progress option (JSON lines) for the GUI to read;
   - a clean end on SIGTERM: the pool terminated, the cancel logged;
   - the run summary from `junctions.tsv`, by the code that writes it.
3. **PR 13, the GUI** (report §4): setup, model, run (`QProcess`),
   results; pytest-qt without a display; one CI job; the `gui` extra on
   `PySide6-Essentials`; README and website.
4. **Small, before or with v0.5.0:**
   - PyPI's trusted publishing: a workflow that uploads when a GitHub
     release is published.
   - Delete the merged `v0.4.0-hardening` branch: the commits cited in
     `reports/` are in `main`'s history.
5. **REVIEW.md §4 LOW items still open:** the terminal-stop tolerance in
   `make_pep`, `chr` never added to chromosome names, the header-only
   fetal-heart test data, the date-stamped comments.
6. **Planned follow-ups** (CHANGELOG `[Unreleased]`):
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
- **Timings:** reading the GENCODE v50 GTF cache takes about 8 s, the
  genome about 35 s. In a whole-genome `-m -c --uniprot -p 16` run, the
  junction loop takes about 3 min, and the run about 4 min.
- **Line endings:** `uniprot_canonicals.tsv` is written with CRLF (the csv
  module's default); the other tables use LF.
