# Handoff — v0.5.0: two items before the release (2026-11-14)

For the next working session. It covers where things stand, how the work
has been run, and what to do next. **Start with
[Next priorities](#next-priorities)**: the maintainer added two items to
v0.5.0 on 2026-11-14, each to begin with an audit.

## State

- **Released:** v0.4.0, on 2026-11-13 (PR #29, merge commit `e4d9117`,
  tag `v0.4.0`, GitHub release, PyPI 0.4.0 built from `258d8b3`). The
  website was updated for it in PR #30 (merged, `f37faf9`).
- **Branch:** `v0.5.0-gui`, from `main` at `f37faf9`. **Pushed** to
  `origin` on 2026-11-14, when the maintainer asked, after their check of
  PR 21 on a Mac. No pull request yet. It collects the v0.5.0 PRs, as
  `v0.4.0-hardening` did for v0.4.0:
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
    its hash;
  - `0986558`: a fix to the tests for CI (report §9.4);
  - then the CHANGELOG and roadmaps updated for the GUI and the two new
    items, with this handoff. Not pushed.
- **Tests:** 230, all passing; 227 offline
  (`python -m pytest -q -m 'not network'`, about 2.5 minutes). Run only
  the fast ones with `-m 'not slow'`. The GUI's 17 (`tests/test_gui.py`)
  need PySide6 and pytest-qt (`pip install -e '.[test-gui]'`; `.conda`
  has both) and are skipped without them.
- **CI:** [.github/workflows/tests.yml](.github/workflows/tests.yml): Python
  3.11 to 3.15 on Linux and macOS, the lowest dependency versions, the
  GUI's tests on macOS (`gui`, new in PR 13), and the installed wheel.
  - **Passes** on `0986558`: run 43663216086, all 13 jobs.
  - The first run on the branch (43662548893, at `b0a806e`) failed 7
    command-line tests in every job but `gui`: Typer styles its text when
    `GITHUB_ACTIONS` is set. Fixed in `tests/conftest.py` (report §9.4).
  - The lowest-versions job prints about 118 deprecation warnings from
    the old matplotlib and urllib3 it installs. They are harmless.
- **Golden tests:** `tests/test_pipeline.py::test_chr15_golden` (GTF-only
  canonicals) and `tests/test_uniprot.py::test_chr15_uniprot`
  (`--uniprot`). Unchanged by PR 18.
- **Done:** roadmap PRs 0–12 and 14–17 (v0.4.0), and 18, 19, 20, 13 and 21
  ([REVIEW.md §9](REVIEW.md#v050-the-gui-and-the-command-line)): the
  v0.5.0 plan of 2026-11-13.
- **Open before the release:** rows 22 and 23, added on 2026-11-14:
  - a view of each junction and its peptides;
  - the `-m` read-count model, assessed.

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
- **Choices the scope did not settle** (report §12.1, "What the GUI
  chose that §4 did not settle"), agreed by the maintainer on 2026-11-15:
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
- **Pushed** after the maintainer ran PR 21 on their Mac. On CI's first
  run, 7 command-line tests failed in every job but `gui`: Typer forces a
  terminal, with styled text, when `GITHUB_ACTIONS` is set. Fixed in the
  tests (report §9.4).
- **v0.5.0 widened** at the maintainer's request, since v0.4.0 is recent.
  Two items come before the release, rows 22 and 23 of the
  [REVIEW.md roadmap](REVIEW.md#v050-the-gui-and-the-command-line):
  - the results view, a view of each junction and its peptides, to be
    designed after a feasibility analysis (it was the v0.5.x junction
    browser);
  - the science of the `-m` read-count model:
    - does it fit RNA counts?
    - a cutoff per splice type?
    - does the cutoff lose bona fide protein-coding isoforms?

    The maintainer noted that the last is hard to test without mass
    spectrometry data.

  The CHANGELOG (one GUI entry for PRs 13 and 21; the plan), REVIEW.md
  (status, §5.6, §9, §10), the scope report, the reports index and the
  README's known issue were updated for it. The briefs are in
  [Next priorities](#next-priorities).
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

Rows 22 and 23 are independent and can run in either order. Each follows
[the agreed workflow](#how-the-work-is-run-agreed-with-the-maintainer):
an audit report with an evidence script, then the maintainer's decisions,
then a PR. The release comes after both.

**The maintainer's answers (2026-11-15):**
- **Mass spectrometry: deferred**, as a future research item. No RNA-seq
  here has spectra from the same samples, not even ENCODE's, and the
  result would depend on the spectra's depth. Instead, row 23 measures at
  each cutoff the known protein-coding isoforms (Swiss-Prot, GENCODE)
  that a JCAST database reproduces.
- **The whole genome's samples** are iPSC-derived cardiomyocytes, 3 DMSO
  (`CM_DMSO1–3`) and 3 treated with the NMD inhibitor NMDI14
  (`CM_NMD1–3`). In `rmats/`, `b1.txt` is DMSO, `b2.txt` NMDI14. Their
  analysis is in the folder's `R/` (`01_DESeq.Rmd`; NMD targets in
  `01_out/`).
- **PR 13's choices** (scope report §12.1): agreed, every one.

### 1. Row 22: a view of each junction and its peptides (feasibility first)

**The goal.** In the GUI's results, pick a junction by gene, splice type,
fate or read count. See its two slices, translated, beside the canonical
protein:
- where the alternative region lies, and the tier;
- the joins to the canonical;
- a T4's premature stop, a T2's frameshifted stretch;
- the tryptic peptides only one isoform has: what a search could
  identify.

**What exists:**
- **`junctions.tsv`:** one row per junction: fate, read count, frames,
  coding region, canonical, files, `new_records`. The Results page reads
  the whole-genome one (599,644 rows) in 2.4 s.
- **The FASTA headers** name the junction: fields 5 and 6, the splice type
  with the slice, and the rMATS id. They also give the exons'
  coordinates, the strand and frame, the read count and the tier.
- **The canonical and known isoforms:** `psq_canonical.fasta` (with
  `-c`); with `--uniprot`, `uniprot_canonicals.tsv` and
  `psq_isoforms.tsv`, the records identical to a UniProt isoform.
- **The exons** are in the rMATS row; the log names the rMATS folder, the
  GTF and the genome.
- Biopython is a dependency already (pairwise alignment, if needed).

**The obstacles** (scope report §4):
- A record whose sequence is already in its file is not written again.
  So a junction's sequence can sit under an earlier junction's header
  (`new_records` is then smaller than its records).
- Nothing records the residues where the alternative region starts and
  ends.

**To measure in the feasibility report:**
- **(a) From the output files alone.** Index the records by splice type
  and rMATS id. Show what was written, and say when a sequence came from
  an earlier junction. No output change, but incomplete.
- **(b) A new output,** a file or columns in `junctions.tsv`: each slice's
  record (or a sequence hash) and the alternative region's residues. This
  is an output change, so it needs an audit and a decision. Measure its
  size on the whole genome and the run time it adds.
- **(c) Translate again in the GUI,** on demand, with the run's own code
  (`Junction`, `sequences`). Loading the GTF cache takes about 8 s and the
  genome about 35 s. This is complete and exact, but needs the inputs
  still in place. Measure the memory.
- **The peptides:** an in-silico tryptic digest (K or R, not before P;
  lengths and missed cleavages to choose). Each peptide's uniqueness
  against the canonical, or against the search database
  (`jcast_uniprot_merged.fasta`). Measure the time per junction and for a
  whole database.
- **The display:**
  - a junction table over `junctions.tsv` (a `QTableView` with a model
    that filters 600,000 rows);
  - an exon diagram (QPainter, or matplotlib as on the Model page);
  - a sequence panel in a monospace font, with the alternative region
    and the unique peptides marked.

  Prototype two or three layouts offscreen and show the screenshots
  (Fusion style; see Useful facts).
- **Out of scope unless asked:** showing the peptides that a search
  identified (for example, FragPipe's `psm.tsv`). This would tie into
  row 23.

**Deliverable:** `reports/<date>-results-view-feasibility.md`, with:
- the measurements and screenshots;
- a recommended design;
- the decisions requested.

Then a PR with a verification script, as for PRs 13 and 21.

### 2. Row 23: the `-m` read-count model, assessed

**The model now** ([jcast/model.py](jcast/model.py), REVIEW.md §5.6):
- **The count** (`junctions.read_count`): min(summed IJC, summed SJC),
  over both samples and every replicate.
- **The model:** a Gamma ("trace") + Normal ("expressed") mixture on
  ln(count + 1), fitted by EM. The cutoff is where the weighted densities
  cross.
- **One fit for the whole folder:** every junction of every rMATS file in
  it, even those `-s` leaves out (PR 12), so every splice type gets the
  same cutoff.
- **Already known to be weak:**
  - the low counts form a lattice (ln 2, ln 3, …);
  - `encode_human_stomach` has no fit;
  - the lung cutoff moves from 11 to 8 when only the EM tolerance
    changes.
- **Pinned in tests:** `tests/test_model.py::test_general_mixture_model`
  pins each bundled dataset's cutoff, and the GUI's tests pin chr15's 26.

**Data:**
- **Bundled:** 13 rMATS folders in `tests/data/model/` (ENCODE human and
  mouse tissues, iPSC day 0 and 14, old mouse heart and kidney), and
  chr15 (`tests/data/rmats`).
- **The maintainer's whole genome,** at
  `/Volumes/Omics-Shared-2/Edward/20240514_Boomi_NMDI_RNASeq/`:
  - `rmats/`: JC and JCEC files;
  - `salmon/`: `CM_DMSO1–3` and `CM_NMD1–3`, each with `quant.sf` (TPM
    per GENCODE transcript) and `quant.genes.sf`;
  - `genome/`: GENCODE v50 GTF, GRCh38, and the Swiss-Prot isoform FASTA
    of 2026-11-11.
- **Two runs on those inputs, both with `--uniprot`:**

  | | `jcast/` | `jcast_guitest/` (from the GUI) |
  |---|---|---|
  | cutoff | `-m`, 40 reads | `-r 7` |
  | skipped as low count | 468,411 of 599,644 | 290,694 |
  | T1 junctions | 24,965 | 50,972 |

**The questions, and ways to answer them:**
- **Does it fit RNA counts?**
  - Compare on the count scale, where the data are. Take the model's
    probability of each count c, which is its density integrated over
    [ln(c + ½), ln(c + 1½)], against the observed frequencies: in
    rootograms, and by log-likelihood or BIC.
  - The models of counts to compare: a two-component negative-binomial
    mixture; a negative binomial (trace) with a Poisson-lognormal
    (expressed).
  - Measure the cutoff's stability by bootstrap, and with other EM starts
    and tolerances.
  - Check what the count is:
    - it is summed over every sample and replicate, so the cutoff grows
      with depth and with the number of samples;
    - it is the smaller of the two forms' counts;
    - IJC and SJC count different numbers of junctions (two inclusion
      junctions for SE). rMATS's `IncFormLen` and `SkipFormLen` allow for
      that; JCAST does not.

    Would a count per sample, or normalized for depth, compare better
    across datasets?
- **A cutoff per splice type?**
  - Fit each type alone, on the bundled folders and the whole genome.
    Compare cutoffs and fit quality.
  - The small types (MXE, A5SS) may have too few junctions for a fit. A
    middle way is a shared trace component with weights per type.
  - Known so far: on chr15, SE alone gives 6 reads, against 26 for every
    type (command-line checks §3).
  - If adopted, it changes the Model page (a plot per type?), the log's
    cutoff line and the pinned tests.
- **Does the cutoff lose bona fide protein-coding isoforms?** Proxies
  that need no mass spectrometry:
  - **Known isoforms.** Preliminary (2026-11-14, from the two runs'
    `psq_isoforms.tsv`): `-r 7` reproduces 5,393 Swiss-Prot isoform
    accessions and `-m 40` 3,999. 1,394 (26%) are reproduced only at
    `-r 7`; none only at `-m 40`.

    Next, run the whole genome at `-r 1` (every junction; about
    15 minutes with `-p 16`). Then, for each cutoff, measure:
    - the isoforms reproduced: recall of the known ones;
    - the new T1 records: the database's growth, which costs search FDR;
    - the T1, T2 and T4 fractions.
  - **Annotation:** junctions whose two forms both match annotated
    protein-coding transcripts (TSL ≤ 2, `Ensembl_canonical`, MANE if the
    GTF tags it), by read count.
  - **Transcript expression:** Salmon's TPM for the transcripts that
    carry each form, in the same samples.
  - **Frame preservation:** a random event keeps the frame about a third
    of the time (T1 against T2). Alternative exons in use are enriched in
    lengths divisible by 3; cite the literature in the report. The count
    above which the frame-keeping fraction levels off is a cutoff that
    comes from the data.
  - **NMD:** the samples are DMSO and the NMD inhibitor NMDI14. Genuine
    premature-stop junctions (T4) should gain reads under the inhibitor.
    Low-count T4 junctions that respond are transcripts, not noise.
  - **With mass spectrometry,** the decisive test. **Deferred** on
    2026-11-15, as a future research item:
    - build databases at several cutoffs and search the same spectra
      with each;
    - count the junction-spanning peptides identified.

    This needs spectra matched to the RNA-seq, and a search engine
    outside JCAST.

**Deliverable:** `reports/<date>-read-count-model-audit.md`, with:
- an evidence script: the bundled folders offline, the whole genome by
  paths;
- the decisions requested, among them:
  - keep the model and document it as a heuristic, or change it;
  - cutoffs per splice type;
  - a recommended `-r`.

Output changes only after a decision.

### 3. Release v0.5.0

Once rows 22 and 23 are done, CI passes, and the maintainer asks. Open a
pull request to `main`, as #29 was for v0.4.0. Then:
- `__version__` (still 0.4.0), the CHANGELOG's date;
- the website: `_quarto.yml`'s title, `change.qmd`, and the GUI's pages;
- before or with the release:
  - PyPI's trusted publishing: a workflow that uploads when a GitHub
    release is published;
  - delete the merged `v0.4.0-hardening` branch: the commits cited in
    `reports/` are in `main`'s history.

### 4. Later

- **GUI follow-ups:**
  - remembering the settings between sessions;
  - the GUI's tests on Linux: Qt's offscreen platform needs system
    libraries there, not tried;
  - screenshots in the README or on the website (they are in
    `reports/gui_scope/verify_13/` and `verify_21/`).
- **REVIEW.md §4 LOW items still open:** the terminal-stop tolerance in
  `make_pep`, `chr` never added to chromosome names, the header-only
  fetal-heart test data, the date-stamped comments.
- **Planned follow-ups** (CHANGELOG `[Unreleased]`):
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
  - Typer styles its text, as on a terminal, when `GITHUB_ACTIONS`,
    `FORCE_COLOR` or `PY_COLORS` is set: read once, at import, into
    `typer.rich_utils.FORCE_TERMINAL`. `tests/conftest.py` sets it to
    `False` for every test.
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
