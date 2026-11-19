# Handoff — v0.5.0: what remains before the release (2026-11-19)

For the next working session. It covers where things stand, how the work
has been run, and what to do next. **Start with
[Next priorities](#next-priorities)**: the maintainer added two items to
v0.5.0 on 2026-11-14, each to begin with an audit, and a third (the GUI
remembers its settings) on 2026-11-18. On 2026-11-19 the first, a view of
each junction and its peptides (row 22), moved to v0.6.0, and the release
waits for the search against spectra that settles row 23's pair of
cutoffs.

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
    items, with this handoff (`3bfb12a`, `da06619`), and row 23's audit
    §§3–7 (`f368783`) and §§8–16 (`5a1c425`);
  - `36f0fb6`: PR 24, JCAST reads rMATS folders made with `--statoff` or
    from one sample (report §17), and the commit recording its hash. Not
    pushed.
- **Tests:** 250, all passing; 247 offline
  (`python -m pytest -q -m 'not network'`, about 2 minutes). Since PR 25,
  a test not marked `network` fails if it connects to the internet. Run only
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
  v0.5.0 plan of 2026-11-13. And PR 24, from row 23's audit.
- **Open before the release:** row 23, added on 2026-11-14, the GUI's
  settings, added on 2026-11-18, and two small PRs chosen on 2026-11-19
  (row 22, a view of each junction and its peptides, moved to v0.6.0):
  - the `-m` read-count model, assessed. Its audit
    ([report](reports/2026-11-15-read-count-model-audit.md)) is measured
    on three data sets (§§3–10, 14, 15). The maintainer answered §13 on
    2026-11-18; the pair of cutoffs waits for a search against spectra
    (§16.1); JCAST reads `--statoff` folders since PR 24 (§17), and the
    maintainer's reruns on `--statoff` folders give its numbers (§17.4);
  - the GUI remembers its settings between sessions;
  - the Model page's two charts that need no run (report §16.2);
  - the three tests that call the UniProt API without `--uniprot`.

## The session of 2026-11-18

- **The maintainer's answers to §13** (report §13, "The maintainer's
  answers"):
  - the `-m` model is kept; open is whether its cutoff should apply to
    novel events only;
  - the pair of cutoffs: leaning yes, to be decided after a search with
    quantms. The maintainer will run JCAST on human heart and prepare
    the spectra of two deep heart proteomes (probably the Mann and
    Lundby labs'). They asked for help with two databases: `-m` on every
    event, and `-m` on novel events only;
  - a cutoff per splice type: dropped for now;
  - the Model page: the database-size charts (§§3.6, 4.5, 5.6), the
    alternative region by read count (§4.4), and the records by kind at
    each cutoff (§4.3) if quick;
  - NMDI14: recorded as inconclusive; a new science report is to follow
    (under [Later](#5-later));
  - the GUI is to remember its settings, in v0.5.0
    ([item 3](#3-the-gui-remembers-its-settings)).
- **They also asked** how the count treats replicates and conditions,
  and how many junctions are lost when one is missing from a condition
  (rMATS gives no p-value), with fetal and adult hearts in mind. A
  fourth data set, ENCODE human adult and fetal heart, is being
  downloaded to `ENCODE_HumanHeart_RNASeq/`, for the search.
- **The third data set (report §14):** `2019_iPSC_Diff_RNA`, 12 samples
  of iPSC to cardiomyocyte differentiation, rMATS 4.4.0, 1.4 million
  events. No conclusion changes:
  - `-m` (39) keeps 37% of the records and 70% of the Swiss-Prot
    isoforms;
  - novel events keep the frame at chance;
  - novel-only wins at every size;
  - the GTF call agrees 98.2%.

  The `-r 1` records reproduce the maintainer's `-m` run exactly.
- **Events rMATS drops (report §15):**
  - JCAST's count sums both conditions. rMATS's post step, with its
    statistics on, drops every event with no reads in one group
    (`filter_countfile`); none reaches JCAST with an NA p-value.
  - Measured against a post step with `--statoff`, JCAST loses 1.9%
    (NMDI14), 2.6% (mouse) and 4.8% (iPSC) of events.
  - In the iPSC series those include the cardiac genes that days 0–2 do
    not express. They would add 148 Swiss-Prot isoforms at `-r 1` and
    12 at `-m`.
  - JCAST cannot read a `--statoff` folder (every FDR `NA` is
    `skipped_fdr`) or a one-group folder (the empty fields make every
    count 0). Decisions requested in §15.7.
- **For the search (report §16.1):**
  [novel_only.diff](reports/model_audit/novel_only.diff) is an
  audit-only switch, `JCAST_AUDIT_ANNOTATED_CUTOFF=<n>`, that gives the
  events rMATS calls annotated their own cutoff. On the iPSC series, a
  run with it at `-m` (39) equals the `-r 1` records cut by the novel-only
  rule, file by file (481,191 junctions; 7,261 Swiss-Prot isoforms
  against 5,556 at `-m` on every event).
- **Tools:** the maintainer's `rmats440` environment is not on this
  machine. rMATS 4.4.0 was installed from bioconda (osx-arm64) into the
  session's scratch folder, for the `--statoff` post steps.
- **New:**
  - [dropped_events.py](reports/model_audit/dropped_events.py);
  - `cutoff_recovery.py --check_novel_only`;
  - the evidence folder `reports/model_audit/ipsc_diff_rmats440/`, and
    `dropped*` files in the three 4.4.0 folders.
- **Later the same day, the maintainer:**
  - **the Model page:** only the charts that need no run (report §16.2:
    junctions kept, annotated and novel; the alternative region by read
    count). Those from a run's records wait for later decisions; they
    would go on the results page or a Model page shown after a run;
  - **§15.7, 1 and 2 adopted:** JCAST is to read `--statoff` and
    one-group folders, and the documentation is to recommend
    `--statoff`. The maintainer will make the three folders again with
    `--statoff` and rerun JCAST after the fix;
  - **asked how far the zero counts move `-m`** (report §15.8). Fitted as
    JCAST fits today, a `--statoff` folder's counts of 0 leave no cutoff
    (1, 1 and 5 reads). Fitted on counts ≥ 1, as the fix will, `-m` falls
    3–5 reads (42 → 37, 25 → 22, 39 → 34), and the database at `-m`
    grows 4.5–4.8%, mostly from events already in the folder.
    `dropped_events.py` now measures at that cutoff too;
- **Committed** as `5a1c425` with the maintainer's approval: everything since
  `f368783` (report §§8–16, the scripts, the evidence, these notes), and
  `.gitignore`'s `runs/` in place of `out/` and `mzml/`.
- **PR 24 (`36f0fb6`, report §17),** at the maintainer's request: JCAST reads rMATS
  folders made with `--statoff` or from one sample.
  - An `NA` FDR is kept by `-q 0 1` only; an empty count field counts 0
    for its sample only; events with no inclusion or no skipping read in
    any sample are left out as the folder is read, and logged.
  - Verified with `reports/model_audit/verify_24.py` against `5a1c425`:
    folders made with statistics (chr15, the 65 bundled files, the three
    whole-genome folders) read and run as before; the `--statoff`
    folders give the statistics run's events plus the one-group events,
    with `-m` 37, 22 and 34; a mouse `--statoff` run at `-r 1` is exactly
    the union of §15.4's two runs.
  - For the check, rMATS's post step was run on one group of the mouse
    samples (`--b1` only): it leaves the second sample's count fields
    empty.
  - 20 new tests (`tests/test_statoff.py`), 14 failing on the code
    before.
- **The maintainer renamed the rMATS folders** while this ran:
  `rmats/` is now `rmats_4_4_0_stat/` (and `jcast/` →
  `jcast_rmats_4_4_0_stat/`, `jcast_4_4_0_stat/` for the iPSC series), to
  make the `--statoff` folders next. The report's commands name the old
  paths.
- **Three tests call the UniProt API without the `network` mark,** as
  they run without `--uniprot`:
  `test_command_line.py::test_the_model_cutoff_with_fewer_splice_types`,
  `test_model.py::test_filter_counts_are_the_runs` and
  `test_model.py::test_the_cutoff_names_its_option` (its `-m` case). On
  2026-11-18 the EBI Proteins API returned HTTP 500, and they sat in
  their retries: 7 and 2 minutes for the last two, which passed; the
  first was stopped and passed on its own with canonicals from the GTF.
  A small fix to propose: give them `--uniprot` with the bundled chr15
  Swiss-Prot, as the other run tests do. Done in PR 25 (2026-11-19), with
  a check in `tests/conftest.py` that fails any test not marked `network`
  that connects to the internet.

## The session of 2026-11-16

- **The maintainer's input:**
  - the model went from a Gaussian mixture to Gamma plus Normal around
    2021, because the counts stop abruptly at 0. Traced: `8c6f015`
    (2021-09-03), `3aff172` (2021-09-16), `ba11fc5` (2021-09-23); report
    §2;
  - the human cells re-run through rMATS 4.4.0 (now `rmats/`; 4.1.0's
    output moved to `rmats_4_1_0/`, and the 2026-11-13 `-m` run to
    `jcast_rmats_4_1_0/`);
  - a second data set, mouse ventricles at five ages (GSE49906, Giudice
    et al. 2014), at `/Volumes/Omics-Shared-2/Edward/GSE49906_MouseHeart_RNASeq/`;
  - asked for canonical genes by cutoff against PaxDb's heart proteins.
- **Row 23's audit, §§4–7** (same report, evidence in
  `reports/model_audit/human_rmats440/` and `mouse_rmats440/`; §3's moved
  to `human_rmats410/`):
  - **rMATS 4.1.0 mislabelled novel events.** 4.4.0 calls 61,246 of
    4.1.0's 151,331 novel-junction events annotated (4.1.1's "novel event
    classification" fix). Per-cutoff results barely change; the novel and
    annotated split of §3 is superseded (§4.1).
  - **Annotation separates better than read count,** more clearly on
    4.4.0. Novel events keep the frame at chance at 1–2 reads in both
    species, and at every count below 4,096 in the human cells (§§4.4,
    5.4).
  - **A cutoff on novel events only** now wins at equal size on every
    measure in both species, but cannot go below about twice
    Swiss-Prot's size (§§4.5, 5.6). Two cutoffs (low for annotated, higher
    for novel) are suggested, not measured (§7).
  - **`-m` lands at the same point:** 40, 42 and 25 reads each keep 35% of
    the database records (§7). A fixed `-r` costs more in mouse (§5.3).
  - **PaxDb (§6):** 0.8% (mouse) and 1.5% (human) of the heart proteins
    have their best junction below 7 reads; 2.7% and 6.7% below `-m`.
    Whether PaxDb finds a gene's protein is flat with read count below
    about 256 reads. At `-m`, a quarter of them get no alternative
    sequence.
  - **Script changes:** `cutoff_recovery.py` gained `--paxdb`,
    `tradeoff.tsv`, `novel.png` and a required `--out`; it reproduces
    §3's tables byte for byte. New: `rmats_versions.py`. The two new
    `-r 1` runs used `da06619` with the audit patch, `-p 10` each.
  - **Committed** as `f368783`, when the maintainer asked. The three
    `-r 1` logs are force-added: `*.log` is ignored.
- **The maintainer's answers, later the same day:**
  - commit, and go on with §8;
  - a cutoff per splice type looks less promising than one for novel
    events and one for annotated events; the novel-only cutoff has strong
    potential;
  - agreed: the GUI's Model page should show what a cutoff keeps, with
    charts like the audit's, so users can choose (other data sets may
    behave differently);
  - a third data set is being prepared: the 2019 paper's deep iPSC
    differentiation series, through rMATS 4.4.0.
- **Row 23's audit, §§8–13** (same report; new scripts `model_fit.py`,
  `two_cutoffs.py`, `nmd_response.py`):
  - **The fit (§8).** The counts have one mode (2 reads) and a long
    tail; `-m` misfits them where its cutoff falls; fitted to the counts
    the same family does not split. The cutoff keeps the 22–26% of
    junctions with the most reads at any depth. JCAST's EM stops before
    converging (converged: 45, 43, 26 against 42, 40, 25).
  - **Two cutoffs (§9).** A low cutoff for annotated events and a higher
    one for novel events beats one cutoff in both species at every size;
    novel-only is best at large sizes. From the GTF, JCAST could call
    novel events itself: 98.0% and 98.8% agreement with rMATS 4.4.0.
  - **NMDI14 (§10): inconclusive.** Neither JCAST's T4 stop forms nor
    GENCODE's NMD forms gain under NMDI14, and in the folder's DESeq2
    results NMD transcripts are no more often up than coding ones.
  - **Per-type cutoffs (§11)** recorded, not pursued; **the Model page
    (§12)**: what each chart would need; **decisions (§13)**.
  - **Not committed:** §§8–13, their evidence and scripts, and these
    notes.

## The session of 2026-11-15

- **The maintainer's answers** (`da06619`; listed under
  [Next priorities](#next-priorities)):
  - PR 13's choices are kept;
  - mass spectrometry is deferred, as a future research item;
  - the whole genome's samples are DMSO and NMDI14.
- **Row 23's audit, §3: which known isoforms a cutoff loses**
  ([report](reports/2026-11-15-read-count-model-audit.md),
  [evidence](reports/model_audit/)).
  - **One run at `-r 1` gives every cutoff.** It is made with an
    audit-only switch, [audit_records.diff](reports/model_audit/audit_records.diff):
    `JCAST_AUDIT_RECORDS=<file>` logs every record of every junction,
    written or not. Cut at 7 and at 40, the records give exactly the two
    earlier whole-genome runs (§3.2). The run took 8 min with `-p 16`.
  - **Recall falls steadily with the cutoff.** `-m` (40) keeps 63% of
    the Swiss-Prot isoforms `-r 1` reproduces, 47% of GENCODE's
    alternative proteins and 81% of those at ≥ 1 TPM (Salmon). It keeps
    28% of the records found in neither. No read count is a boundary.
  - **Annotated against novel events** (rMATS's
    `fromGTF.novelJunction` and `novelSpliceSite` lists):
    - annotated events with 1 read look like those with 100;
    - novel events with 1–2 reads keep the frame at chance (33.3%).
  - **A cutoff on novel events only** keeps 175–376 more Swiss-Prot
    isoforms for the same database size. It ties on expressed GENCODE
    isoforms up to about 20 reads and loses some above that (§3.6).
  - **The paper's model** (a Gaussian mixture, cutoff at a posterior of
    0.95) is not the code's (Gamma plus Normal, densities crossing; §2).
  - **Not committed yet:** the report, its evidence and
    `reports/model_audit/cutoff_recovery.py`. Their index, REVIEW.md and
    handoff entries are not committed either. Commit them with the
    decisions, as for other audits, or before if the maintainer asks.

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

**The order, decided 2026-11-19:**
1. Three small PRs: the tests that call the UniProt API get `--uniprot`;
   the GUI remembers its settings (item 3); the Model page's two charts
   that need no run (report §16.2).
2. The search against spectra (row 23, report §16.1): the maintainer's
   human heart folder (ENCODE, `ENCODE_HumanHeart_RNASeq`, downloading
   on 2026-11-19) is expected within a day. Two databases from it: `-m`
   as JCAST runs it, and `-m` on novel events only
   (`novel_only.diff`, not in JCAST), for the maintainer's quantms
   search.
3. The maintainer decides from the search whether the pair of cutoffs
   (`-m` on novel events, every annotated event kept) goes into v0.5.0.
4. Then the release (item 4).

Row 22 moved to v0.6.0: it needs a feasibility audit and some thought.
v0.5.0 is to make its results solid, and to answer the questions the
read-count audit raised.

Each item follows
[the agreed workflow](#how-the-work-is-run-agreed-with-the-maintainer):
an audit report with an evidence script, then the maintainer's decisions,
then a PR.

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

**The maintainer's answers (2026-11-16):** row 23 favours a cutoff for
novel events (with one for annotated events) over a cutoff per splice
type; the GUI's Model page is to show what a cutoff keeps; a third data
set (the 2019 paper's iPSC differentiation series, rMATS 4.4.0) is to
come.

**The maintainer's answers (2026-11-18):** listed in
[the session of 2026-11-18](#the-session-of-2026-11-18) and the report's
§13. In short:
- `-m` is kept;
- the pair of cutoffs waits for the search against spectra;
- no per-type cutoffs;
- the Model page charts are chosen;
- the GUI remembers its settings in v0.5.0.

### 1. Row 22: a view of each junction and its peptides (moved to v0.6.0)

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

**Progress (2026-11-18).** The audit report is
[reports/2026-11-15-read-count-model-audit.md](reports/2026-11-15-read-count-model-audit.md).
Measured on three data sets (§§3–10, 14):
- the known isoforms a cutoff loses: on rMATS 4.1.0 (§3), 4.4.0 (§4),
  mouse ventricles (§5) and iPSC differentiation (§14);
- the PaxDb heart proteins (§6), the fit (§8), and two cutoffs (§9);
- NMDI14 (§10, inconclusive);
- the events rMATS's statistics drop (§15).

The maintainer answered §13 on 2026-11-18. Next:
- **§15.7's fix: done in PR 24 (2026-11-18, report §17).** JCAST reads
  `--statoff` and one-sample folders, and the README recommends
  `--statoff`. The maintainer remade the three folders with `--statoff`
  and reran JCAST on 2026-11-18, checked on 2026-11-19 (§17.4,
  `check_reruns.py`): the cutoffs (37, 22, 34), database records and
  Swiss-Prot isoforms are §15.8's, and every junction and sequence of
  the earlier runs is kept;
- **the search against spectra (§16.1):** when the maintainer's heart
  folder is ready, two runs from it: `-m`, and `-m` with
  `novel_only.diff` and `JCAST_AUDIT_ANNOTATED_CUTOFF=1`. The folder is
  to be made with `--statoff`, and both runs made after §15.7's fix;
- **the ENCODE fetal and adult heart folder:** §15's measurement, run
  both ways;
- **then the PRs:**
  - the §15.7 fixes (done, PR 24);
  - the Model page charts that need no run (§16.2); those from a run's
    records wait for later decisions;
  - the pair of cutoffs, if the search supports it.

The `-r 1` runs' records (166, 174, 93 and 375 MB) were in a scratch folder.
To regenerate them, apply `reports/model_audit/audit_records.diff` in a
worktree and run as the report's Reproduce row says. GENCODE translations
and PaxDb files are downloaded as its Inputs row says.

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
  - `rmats/`: rMATS 4.4.0 with `--statoff` (2026-11-18), and `jcast/`
    the `-m` run on it; the folders made with statistics (2026-11-16)
    are `rmats_4_4_0_stat/` and `jcast_rmats_4_4_0_stat/`, and 4.1.0's
    `rmats_4_1_0/`. The mouse folder is laid out the same way; the iPSC
    series' earlier run is `jcast_4_4_0_stat/`. In the repository,
    `runs/<name>/` (ignored by git) links to each folder's `rmats/`,
    `jcast/`, `genome/` and `multiqc/`;
  - `salmon/`: `CM_DMSO1–3` and `CM_NMD1–3`, each with `quant.sf` (TPM
    per GENCODE transcript) and `quant.genes.sf`;
  - `genome/`: GENCODE v50 GTF, GRCh38, and the Swiss-Prot isoform FASTA
    of 2026-11-11.
- **Two runs on the rMATS 4.1.0 inputs, both with `--uniprot`** (the
  first now in `jcast_rmats_4_1_0/`, the second removed; `jcast/` holds
  the maintainer's `-m` run on 4.4.0, cutoff 42):

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

### 3. The GUI remembers its settings

Decided 2026-11-18 (it was under Later): the GUI keeps its settings
between sessions. A small PR: `QSettings` for the Setup page's paths and
options, restored at start, with a way back to the defaults. Its test
should use a temporary settings file, never the user's.

### 4. Release v0.5.0

Once row 23 (with the decision on the pair of cutoffs) and the small
PRs are done, CI passes, and the maintainer asks. Open a
pull request to `main`, as #29 was for v0.4.0. Then:
- `__version__` (still 0.4.0), the CHANGELOG's date;
- the website: `_quarto.yml`'s title, `change.qmd`, and the GUI's pages;
  PR 24 changed `about.qmd` and `faqs.qmd` (not rendered), and
  `install.qmd`'s rMATS section still describes rMATS turbo 0.1;
- before or with the release:
  - PyPI's trusted publishing: a workflow that uploads when a GitHub
    release is published;
  - delete the merged `v0.4.0-hardening` branch: the commits cited in
    `reports/` are in `main`'s history.

### 5. Later

- **A science report on NMDI14** (the maintainer's to-do, 2026-11-18): what
  the treatment does to splicing, from rMATS (`rmats/`), and to gene
  expression, from the folder's Salmon output, which the maintainer will
  link in. The read-count audit's §10 found no NMD block; this would ask
  what NMDI14 does instead. Not part of v0.5.0.
- **GUI follow-ups:**
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
