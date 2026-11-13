# Handoff — v0.4.0 hardening (2026-11-12, after PRs 16–17)

For the next working session on the `v0.4.0-hardening` branch. It covers
where things stand, how the work has been run, and what to do next.

## State

- **Branch:** `v0.4.0-hardening`, last commit at handoff in `git log`.
  `origin` has it up to `8b845aa` (before PRs 6–8); later commits are
  local until the maintainer asks for a push.
- **Tests:** all 148 pass in about 2.5 minutes (`python -m pytest -q`; run
  only the fast ones with `-m 'not slow'`).
- **Golden tests:**
  - `tests/test_pipeline.py::test_chr15_golden` (GTF-only canonicals). Its
    `read=1` keeps every junction since PR 6.
  - `tests/test_uniprot.py::test_chr15_uniprot` (`--uniprot`, on the
    Swiss-Prot subset in `tests/data/uniprot/`).
- **Done:** roadmap PRs 0–9 and 14–17 ([REVIEW.md §9](REVIEW.md#9-v040-hardening-roadmap-revised)).
  Every output-changing PR has an audit, an evidence folder and a
  verification script, indexed in [reports/README.md](reports/README.md)
  with its decisions and commit hashes.
- **Recent:**
  - **PR 16** ([log audit §7](reports/2026-11-12-log-audit.md#7-implementation-and-verification-pr-16-a85f8fc)):
    - `junctions.tsv`: one row per junction, in rMATS order, the same
      with any `-p` (fate, frames, coding region, canonical, files
      written).
    - The log is small: settings, inputs, warnings, errors and an
      end-of-run summary (5 kB for a whole genome, was 346 MB).
      `-v/--verbose` writes the per-junction detail (300 MB).
    - No output file changed.
  - **PR 17** ([`Ensembl_canonical` audit §10](reports/2026-11-12-ensembl-canonical-audit.md#10-implementation-and-verification-pr-17-8ddf8c3)):
    - When no transcript passes the GTF rule, the GTF canonical is the
      `Ensembl_canonical` transcript, if it is protein_coding, complete
      and not a readthrough transcript.
    - 17 genes of the whole-genome run gain a canonical. No fate
      changes, and the merged database is unchanged.
    - Readthrough genes get no GTF canonical, by design (the
      maintainer's choice: they are not JCAST's purpose, and their slices
      would join to fusion proteins). They are 370 of the 371 coding
      genes still without one.
- **For the user's next whole-genome run:** the log is now small; the
  per-junction outcome is in `junctions.tsv`. Pass `-v` only to debug.

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

1. **Smaller items:**
   - PR 12, argparse hygiene: drop `--g_or_ln`; `-s` reads only the
     selected files (MED-8); `parser.error`, not `assert`.
   - PR 10 (packaging, Python ≥ 3.10) and PR 11 (CI).
   - REVIEW.md §4 LOW items still open: the terminal-stop tolerance in
     `make_pep`, `chr` never added to chromosome names, the glob-built
     `__all__`, the unused `get_rest_nuc`, the header-only fetal-heart
     test data.
2. **PR 13, GUI:**
   - Expose `--uniprot`, `--overwrite`, `--cache_dir` and `--verbose`.
   - `-p` must not fork from the GUI's Qt thread (REVIEW.md §10). Run a
     subprocess or use another start method.
3. **Planned follow-ups** (CHANGELOG "Planned"):
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
  output unless `args.overwrite`. It closes its log handlers when it ends.
- **Tests:**
  - `tests/conftest.py` points the caches at a session temp directory (the
    `cache_dir` fixture).
  - Tests that call `run_jcast` must pass `cache_dir=str(cache_dir)` in
    `args`, or they write to the user's cache.
  - `tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf.cached`
    (git-ignored, from before PR 7) is no longer read and can be deleted.
- **Verification scripts:**
  - They take input paths as arguments and write only to the scratch
    directory.
  - `reports/runs_audit/verify_6.py --pair 'label|repo|args|repo|args'`
    compares two CLI runs in any two checkouts.
  - For whole-genome runs, pass `--cache_dir`, or set `HOME`, to keep the
    rebuilt GTF cache out of the user's cache.
  - `junctions.tsv` makes run comparisons cheap: fates, canonicals and
    records per junction, without parsing the log.
- **Timings:** reading the GENCODE v50 GTF cache takes about 8 s, the
  genome about 35 s. In a whole-genome `-m -c --uniprot -p 16` run, the
  junction loop takes about 3 min, and the run about 4 min.
- **Line endings:** `uniprot_canonicals.tsv` is written with CRLF (the csv
  module's default); the other tables use LF.
