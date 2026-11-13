# Handoff — v0.4.0 hardening (2026-11-13, after PRs 10–12)

For the next working session on the `v0.4.0-hardening` branch. It covers
where things stand, how the work has been run, and what to do next.

## State

- **Branch:** `v0.4.0-hardening`, last commit at handoff in `git log`.
  Check `git status -sb` for what `origin` has; push only when the
  maintainer asks.
- **Tests:** all 153 pass in about 2.5 minutes (`python -m pytest -q`; run
  only the fast ones with `-m 'not slow'`, offline with `-m 'not network'`).
- **CI:** [.github/workflows/tests.yml](.github/workflows/tests.yml) (PR 11)
  runs on every push and pull request. Its first run on GitHub
  (43455685837, at `efed80b`) passed all 12 jobs.
- **Golden tests:**
  - `tests/test_pipeline.py::test_chr15_golden` (GTF-only canonicals). Its
    `read=1` keeps every junction since PR 6.
  - `tests/test_uniprot.py::test_chr15_uniprot` (`--uniprot`, on the
    Swiss-Prot subset in `tests/data/uniprot/`).
- **Done:** roadmap PRs 0–12 and 14–17 ([REVIEW.md §9](REVIEW.md#9-v040-hardening-roadmap-revised)).
  Every output-changing PR has an audit, an evidence folder and a
  verification script, indexed in [reports/README.md](reports/README.md)
  with its decisions and commit hashes.
- **Recent (2026-11-13):**
  - **PR 10** ([packaging check](reports/2026-11-13-packaging-check.md)):
    - The metadata is in `pyproject.toml`; `setup.py` and
      `requirements.txt` are gone. The version is set once, in
      `jcast/__init__.py`: 0.4.0.
    - Python ≥ 3.11. The tests pass on 3.11 to 3.15 with the newest
      dependencies, and on 3.11 with the lowest the bounds allow.
    - Only `jcast` is installed (the wheel is 62 kB); the sdist holds the
      tests and their data (60 MB).
  - **PR 11:** CI, as above: 10 test jobs, a lowest-versions job and a
    package job.
  - **PR 12** ([command-line checks](reports/2026-11-13-command-line-checks.md)):
    - Only the `-s` files must exist (MED-8). With `-m`, the others
      present are read too, so the cutoff is unchanged.
    - A missing file is a usage error; `--g_or_ln`, `get_rest_nuc` and the
      scikit-learn Gaussian mixture are removed.
    - No output change for any run that worked before.

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

1. **The v0.4.0 release** (the maintainer's plan, 2026-11-13):
   - Done: the release commit (CHANGELOG `[0.4.0] — 2026-11-13`, its
     Planned items moved to a new `[Unreleased]`; README links made
     absolute, since the README is also the PyPI page); the pull request
     from `v0.4.0-hardening` into `main` (`gh pr view`).
   - After the merge: tag `v0.4.0` on `main`, and the GitHub release, with
     the wheel and sdist attached.
   - PyPI: the maintainer uploads, with their own token. Build the files
     from the tagged commit, not the working tree (`git archive`, then
     `uv build`), and run `twine check --strict`. setuptools ≥ 77 writes
     Metadata-Version 2.4: twine ≥ 6.1 reads it, older twine only with
     pkginfo ≥ 1.12. A version can be uploaded to PyPI only once.
   - The website (`docs/`) still says Python 3.7+ and `tox`. It is built
     from `jcast_website/`, which is not in the repository.
2. **PR 13, GUI:**
   - Expose `--uniprot`, `--overwrite`, `--cache_dir` and `--verbose`.
   - `-p` must not fork from the GUI's Qt thread (REVIEW.md §10). Run a
     subprocess or use another start method. On Python ≥ 3.12, forking a
     multi-threaded process warns ([packaging check §4.2](reports/2026-11-13-packaging-check.md#42-notes)).
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
- **The maintainer's `.conda` environment is x86_64**, running under
  Rosetta on an arm64 Mac (Python 3.11, pandas 2.3, numpy 1.26). A native
  arm64 environment would likely be faster.
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
