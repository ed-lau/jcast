# Handoff — v0.4.0 hardening (2026-11-12, after PRs 6–8)

For the next working session on the `v0.4.0-hardening` branch. It covers
where things stand, how the work has been run, and what to do next.

## State

- **Branch:** `v0.4.0-hardening`, last commit at handoff in `git log`.
- **Tests:** all 114 pass in about 2.5 minutes (`python -m pytest -q`; run
  only the fast ones with `-m 'not slow'`).
- **Golden tests:**
  - `tests/test_pipeline.py::test_chr15_golden` (GTF-only canonicals). Its
    `read=1` keeps every junction since PR 6.
  - `tests/test_uniprot.py::test_chr15_uniprot` (`--uniprot`, on the
    Swiss-Prot subset in `tests/data/uniprot/`).
- **Done:** roadmap PRs 0–9, 14 and 15 ([REVIEW.md §9](REVIEW.md#9-v040-hardening-roadmap-revised)).
  Every output-changing PR has an audit, an evidence folder and a
  verification script, indexed in [reports/README.md](reports/README.md)
  with its decisions and commit hashes.
- **Recent** ([PR 6–8 audit](reports/2026-11-12-thresholds-runs-canonical-audit.md)
  §8–10):
  - **PR 7:**
    - `-o` is created. A directory holding JCAST output is refused unless
      `--overwrite`.
    - A failing junction gets `FAILURE (exception)` and the run continues.
    - Caches go in the user cache directory (`platformdirs`) or
      `--cache_dir`; the GTF cache is a pickle, about 5 min to build for
      GENCODE.
    - Logger handlers no longer stack.
  - **PR 6:**
    - Keep a junction if its count ≥ the cutoff. The default `-r 1` keeps
      every junction (the maintainer's choice).
    - `-m` cuts at the model's boundary, in reads.
    - `sum_sjc` is now `min_jc`.
    - Filters run first: a whole-genome `-m` run's junction loop takes
      240 s, down from 748 s.
  - **PR 8:**
    - `-c` writes every gene's canonical.
    - `write_canonical`'s API fallback is removed.
    - A failed or refused UniProt API request fails the junction.
    - `--uniprot` uses gene names only when unique in the GTF. That
      changed 77 genes, 72 of them non-coding
      ([mapping audit §8](reports/2026-11-12-uniprot-mapping-audit.md#8-addendum-the-unique-name-rule-on-the-whole-genome-pr-8)).
- **For the user's next whole-genome run:**
  - Rerunning into an existing output directory needs `--overwrite`.
  - The first run rebuilds the GTF cache once, in `~/Library/Caches/jcast`.
  - The old `gencode.v50…gtf.cached` next to the GTF can be deleted.
  - With `-m`, expect cutoff 40 and 131,233 junctions translated.

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

1. **PR 16: the run record** (decided;
   [log audit §6](reports/2026-11-12-log-audit.md#6-decisions-2026-11-12)):
   - a per-junction `junctions.tsv`;
   - a small log with an end-of-run summary (fates per splice type,
     records per file, failures, time);
   - sequences only with `--verbose`;
   - no `~module~function` suffix.
2. **PR 17: the GTF canonical from GENCODE's `Ensembl_canonical` tag**
   when no transcript passes the TSL ≤ 2 / start-and-stop filter. It
   covers 173 of the 384 coding genes left without a canonical. It changes
   output: audit first
   ([mapping audit §5–7](reports/2026-11-12-uniprot-mapping-audit.md)).
3. **Smaller items:**
   - PR 12, argparse hygiene: drop `--g_or_ln`; `-s` reads only the
     selected files (MED-8); `parser.error`, not `assert`.
   - PR 10 (packaging, Python ≥ 3.10) and PR 11 (CI).
4. **PR 13, GUI:**
   - Expose `--uniprot`, `--overwrite` and `--cache_dir`.
   - `-p` must not fork from the GUI's Qt thread (REVIEW.md §10). Run a
     subprocess or use another start method.
5. **Planned follow-ups** (CHANGELOG "Planned"):
   - an option to add T2/T4 to `jcast_uniprot_merged.fasta`;
   - a UniProt download helper that names the file by its release
     (`X-UniProt-Release` header);
   - removing `--mask` in the next release;
   - joining by genomic coordinates was evaluated and deferred
     ([UniProt evaluation §4](reports/2026-11-11-uniprot-fasta-evaluation.md#4-joining-by-coordinates-question-1)).
   - Possibly: prefer the gene with an HGNC id when a gene name is shared
     ([mapping audit §8](reports/2026-11-12-uniprot-mapping-audit.md#8-addendum-the-unique-name-rule-on-the-whole-genome-pr-8)).
     Not decided.

## Useful facts

- `run_jcast` creates the output directory and refuses one holding JCAST
  output unless `args.overwrite`. It closes its log handlers when it ends.
- **Tests:**
  - `tests/conftest.py` points the caches at a session temp directory (the
    `cache_dir` fixture).
  - Tests that call `run_jcast` must pass `cache_dir=str(cache_dir)` in
    `args`, or they write to the user's cache.
- **Verification scripts:**
  - They take input paths as arguments and write only to the scratch
    directory.
  - `reports/runs_audit/verify_6.py --pair 'label|repo|args|repo|args'`
    compares two CLI runs in any two checkouts.
  - For whole-genome runs, pass `--cache_dir`, or set `HOME`, to keep the
    rebuilt GTF cache out of the user's cache.
- **Timings:** reading the GENCODE v50 GTF cache takes about 8 s, the
  genome about 35 s. In a whole-genome `-m -c --uniprot -p 16` run, the
  junction loop takes 4–5 min.
