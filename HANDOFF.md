# Handoff — v0.4.0 hardening (2026-11-12)

For the next working session on the `v0.4.0-hardening` branch. It covers
where things stand, how the work has been run, and what to do next.

## State

- **Branch:** `v0.4.0-hardening`, last commit at handoff in `git log`.
- **Tests:** all 82 pass in about 2 minutes (`python -m pytest -q`; run
  only the fast ones with `-m 'not slow'`).
- **Golden tests:**
  - `tests/test_pipeline.py::test_chr15_golden` (GTF-only canonicals);
  - `tests/test_uniprot.py::test_chr15_uniprot` (`--uniprot`, on the
    Swiss-Prot subset in `tests/data/uniprot/`).
- **Done:** roadmap PRs 0–5, 9, 14 and 15 ([REVIEW.md §9](REVIEW.md#9-v040-hardening-roadmap-revised)).
  Every output-changing PR has an audit, an evidence folder and a
  verification script, indexed in [reports/README.md](reports/README.md)
  with its decisions and commit hashes.
- **Recent:**
  - `--uniprot FASTA` (PR 14): canonicals from a supplied UniProt FASTA,
    GTF fallback, built-in Swiss-Prot merge.
  - A per-gene GTF index and in-memory de-duplication of writes (PR 9).
  - `-p/--processes` (PR 15).
  - A whole-genome run (about 600,000 junctions) now takes about 15 min
    with `-p 16`, where it took days.
- **In progress on the user's side:** a full whole-genome run with
  `--uniprot -p 16` (NMDI RNA-seq data on
  `/Volumes/Omics-Shared-2/Edward/20240514_Boomi_NMDI_RNASeq`). Ask how it
  went; anything unexpected is the first thing to look at.

## How the work is run (agreed with the maintainer)

1. **Audit first.**
   - Write a dated report in `reports/`, with a script that regenerates
     its evidence.
   - Measure on the chr15 test set. For performance, use a whole-genome
     sample.
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

1. **Follow up on the user's full run** (above).
2. **PR 7: robustness of real runs.** Most pressing now that full runs are
   routine.
   - Create the output directory (MED-4). Today it must already exist.
   - Fresh output per run (MED-12). Output is appended, and sequences
     already in a file are skipped, so a rerun into the same directory
     mixes runs.
   - Catch an exception per junction, record it as a fate and continue
     (decided in the masked-bases audit §7). With `-p`, an exception in
     a worker currently aborts the run, as it does serially.
   - Logger handlers accumulate across calls (MED-13); the tests clear
     them by hand.
   - The UniProt API cache is written to the working directory (MED-5);
     the GTF `.cached` file sits next to the GTF (MED-6).
3. **PR 6: threshold semantics** (MED-11: read-count cutoffs off by one)
   and `sum_sjc` naming (MED-1 / §5.2). This changes output, so audit
   first.
4. **PR 8:**
   - `-c` writes the canonical on every path (MED-14; 2 chr15 genes).
   - `write_canonical`'s no-op API retry (MED-9).
   - API network failures (MED-10; avoidable with `--uniprot`).
   This changes output, so audit first.
5. **Smaller items:**
   - PR 12, argparse hygiene: drop `--g_or_ln`; `-s` reads only the
     selected files (MED-8); `parser.error`, not `assert`.
   - The rest of PR 9: filter junctions before the canonical lookup.
     Minor now.
   - PR 10 (packaging, Python ≥ 3.10) and PR 11 (CI).
6. **PR 13, GUI:**
   - Expose `--uniprot`.
   - `-p` must not fork from the GUI's Qt thread (REVIEW.md §10). Run a
     subprocess or use another start method.
7. **Planned follow-ups** (CHANGELOG "Planned"):
   - an option to add T2/T4 to `jcast_uniprot_merged.fasta`;
   - a UniProt download helper that names the file by its release
     (`X-UniProt-Release` header);
   - removing `--mask` in the next release;
   - joining by genomic coordinates was evaluated and deferred
     ([UniProt evaluation §4](reports/2026-11-11-uniprot-fasta-evaluation.md#4-joining-by-coordinates-question-1)).

## Useful facts

- `run_jcast` does not create the output directory. Tests and scripts
  create it first.
- In tests, `tests/conftest.py` moves the working directory to a temp
  dir, because of the UniProt cache. A run of `run_jcast` in-process should
  clear `logging.getLogger('jcast').handlers` before and after.
- The GENCODE v50 `.cached` GTF on the user's volume takes about 55 s to
  read, and the genome about 35 s. Verification scripts take input paths
  as arguments and write only to the scratch directory.
