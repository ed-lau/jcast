# Handoff — after the v0.4.0 release (2026-11-13)

For the next working session. It covers where things stand, how the work
has been run, and what to do next.

## State

- **Released:** v0.4.0, on 2026-11-13.
  - PR #29 merged `v0.4.0-hardening` into `main` with a merge commit,
    `e4d9117`, which the tag `v0.4.0` points to.
  - GitHub release "JCAST v.0.4.0".
  - PyPI 0.4.0: the wheel and the sdist built from `258d8b3` (sha256
    `97edb11d…` and `bd5a4b4c…`), whose tree is the tag's. They were
    uploaded three minutes before the merge, with no effect on the
    contents.
- **Branch:** `docs-after-0.4.0`, from `main`, with the documentation
  changes made after the release (below). Pushed, with a PR into `main`
  (`gh pr list`). The website updates when it is merged. Start v0.5.0 work
  on a new branch from `main` after that merge.
- **Tests:** all 153 pass in about 2.5 minutes (`python -m pytest -q`; run
  only the fast ones with `-m 'not slow'`, offline with `-m 'not network'`).
- **CI:** [.github/workflows/tests.yml](.github/workflows/tests.yml) runs on
  every push and pull request: Python 3.11 to 3.15 on Linux and macOS, the
  lowest dependency versions, and the installed wheel. Green since its
  first run (43455685837).
- **Golden tests:**
  - `tests/test_pipeline.py::test_chr15_golden` (GTF-only canonicals). Its
    `read=1` keeps every junction since PR 6.
  - `tests/test_uniprot.py::test_chr15_uniprot` (`--uniprot`, on the
    Swiss-Prot subset in `tests/data/uniprot/`).
- **Done:** roadmap PRs 0–12 and 14–17 ([REVIEW.md §9](REVIEW.md#9-v040-hardening-roadmap-revised)).
  Every output-changing PR has an audit, an evidence folder and a
  verification script, indexed in [reports/README.md](reports/README.md)
  with its decisions and commit hashes.

## The session of 2026-11-13

- **PR 10, packaging** (`4a039a9`; [packaging check](reports/2026-11-13-packaging-check.md)):
  - The metadata is in `pyproject.toml`; `setup.py` and `requirements.txt`
    are gone. The version is set once, in `jcast/__init__.py`.
  - Python ≥ 3.11. The tests pass on 3.11 to 3.15 with the newest
    dependencies, and on 3.11 with the lowest the bounds allow. The
    numpy, scipy, matplotlib and requests floors were raised to versions
    that install and work.
  - The wheel installs only `jcast` (60 kB). The sdist holds the tests and
    their data (60 MB).
- **PR 11, CI** (`dee41b6`): as above, 12 jobs.
- **PR 12, command-line checks** (`d84f49d`; [report](reports/2026-11-13-command-line-checks.md)):
  - Only the `-s` files must exist (MED-8). With `-m`, the others present
    are read too, so the cutoff is unchanged.
  - A missing file is a usage error. `--g_or_ln`, `get_rest_nuc` and the
    scikit-learn Gaussian mixture are removed.
  - No output change for any run that worked before.
- **The release:**
  - The release commit `258d8b3` dated the CHANGELOG and made the README's
    links absolute, since the README is also the PyPI page.
  - PR #29, then the tag, the GitHub release and the PyPI upload, which
    the maintainer did.
- **After the release** (branch `docs-after-0.4.0`, PR into `main`):
  - The README says the GUI is a preview that cannot run JCAST. Two FAQ
    answers were out of date: the rMATS files needed (PR 12) and the GTF
    cache next to the GTF (PR 7).
  - **The website** was updated for v0.4.0. Its pages had described v0.3.5
    (Python 3.7+, `tox`, `--g_or_ln`, Tier 3). Its Quarto source,
    `jcast_website/`, is now tracked (maintainer's decision). The
    Documentation and FAQ pages take their sections from the README.
  - REVIEW.md's status records the release.
- **The maintainer's environment:** `.conda` was rebuilt as native arm64
  (Python 3.11.18, numpy 2.4, pandas 3.1, `jcast` installed editable,
  `subdir: osx-arm64` pinned in `.conda/.condarc`). `~/miniconda3` itself
  is still the Intel build, so its other environments run under Rosetta.
- **A mistake, undone:** a `uv pip install --system` went into
  `~/miniconda3` through `CONDA_PREFIX` (see Useful facts).

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

1. **PR 13, GUI (v0.5.0)** ([REVIEW.md §10](REVIEW.md#10-gui-plan--notes)):
   - The GUI released in 0.4.0 is a skeleton: three placeholder pages,
     157 lines in `jcast/gui/`. It cannot run JCAST.
   - Planned: a wizard (inputs, filters), the model preview, and a run
     panel that runs `run_jcast` in a `QThread`. Expose `--uniprot`,
     `--overwrite`, `--cache_dir` and `--verbose`.
   - `-p` must not fork from the GUI's process, because Qt runs threads.
     Run a subprocess or use another start method. On Python ≥ 3.12,
     forking a multi-threaded process warns ([packaging check §4.2](reports/2026-11-13-packaging-check.md#42-notes)).
   - Also needed: a cancellation hook in `run_jcast`, and the model plot
     drawn on a caller-supplied `Axes`.
   - The next release is v0.5.0 (the maintainer, 2026-11-13).
2. **Before or with v0.5.0**, all small:
   - PyPI's trusted publishing: a workflow that uploads when a GitHub
     release is published, so the upload follows the tag and needs no
     token.
   - The `-h` text has two typos from before v0.4.0: "junctionsand" and
     "transcriptslices". The website shows the `-h` text verbatim, so
     render the site again after fixing them.
   - The merged `v0.4.0-hardening` branch can be deleted: the commits
     cited in `reports/` are in `main`'s history. The 0.4.0 GitHub release
     (fixed by the maintainer) links to the tag, not the branch.
3. **REVIEW.md §4 LOW items still open:** the terminal-stop tolerance in
   `make_pep`, `chr` never added to chromosome names, the header-only
   fetal-heart test data, the date-stamped comments.
4. **Planned follow-ups** (CHANGELOG `[Unreleased]`):
   - an option to add T2/T4 to `jcast_uniprot_merged.fasta`;
   - a UniProt download helper that names the file by its release
     (`X-UniProt-Release` header);
   - removing `--mask` in the next release;
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
