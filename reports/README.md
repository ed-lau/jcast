# JCAST audit reports

Scientific-integrity audits behind the v0.4.0 hardening series. Each
output-changing PR has an audit, written before or alongside the code
change, with:

- the evidence files and the script that regenerates them;
- regression tests that fail before the change;
- the commit at which the evidence reproduces exactly.

The roadmap is [REVIEW.md §9](../REVIEW.md#9-v040-hardening-roadmap-revised).
User-facing changes are in [CHANGELOG.md](../CHANGELOG.md).

All audits use the bundled chr15 test set (1,105 rMATS junctions, Ensembl
GRCh38.89) in GTF-only mode with UniProt stubbed out, so they need no
network.

## Audits

| Report | PR | Status | Evidence reproduces at |
|---|---|---|---|
| [tx0/tx1 translation bounds](2026-11-03-tx1-translation-bounds-audit.md) | 4 | Implemented in `da72a0a` | `a555446` (`python reports/tx1_audit/audit_tx1.py`) |
| [Trimming to the coding region (UTR exons)](2026-11-10-utr-exon-trimming-audit.md) | 4a | Implemented in `89e1a82` | `89e1a82` (`python reports/utr_audit/audit_utr.py`) |
| [RI retained-intron slices](2026-11-10-ri-slice-audit.md) | 4a2 | Implemented in `bc7d894` | `bc7d894` (`python reports/ri_audit/audit_ri.py`) |
| [Reading frame and three-frame fallback](2026-11-10-phase-audit.md) | 4b | Candidate `9e27815`; final `2679edb` | four variants at `9e27815`; final check with `--final` |
| [Frameshift test from coding lengths](2026-11-10-frameshift-audit.md) | 4b2 | Implemented (see the reports' commit log) | the PR 4b2 commit (`python reports/frameshift_audit/audit_frameshift.py`) |
| [Alternative C-termini (first evaluation)](2026-11-10-alt-cterm-tiering-audit.md) | 4c | Superseded by the next row. Measured on PR 4 alone. | `a555446` (`python reports/cterm_audit/audit_cterm.py`) |
| [Alternative C- and N-termini on the current pipeline](2026-11-10-alt-termini-audit.md) | 4c (and 4b2) | **Decisions approved 2026-11-10** (§7). Also found a PR 4b regression in the frameshift test (§6), proposed as PR 4b2. | `919f6c7` (`python reports/termini_audit/audit_termini.py`) |

Each script compares a baseline commit with the working tree, or patches
the candidate in at runtime, so later code changes its output. To reproduce
a report's evidence files exactly, check out the commit in the last column
first.

## Decisions log

| Date | Decision | Evidence |
|---|---|---|
| 2026-11-10 | Drop exons lying entirely in the UTR; skip events entirely outside the coding region (new fate) | [4a audit](2026-11-10-utr-exon-trimming-audit.md) |
| 2026-11-10 | Defer one coding region per gene vs per transcript | [tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions) item 2 |
| 2026-11-10 | Second-start-site retry: deferred, then **removed in PR 4a**, because with UTR-exon removal it deleted exons coding in TSL 1–2 transcripts (20 junctions) | [4a audit §4.3](2026-11-10-utr-exon-trimming-audit.md#43-why-the-retry-had-to-be-removed) |
| 2026-11-10 | Read each slice's frame at its 5′ end from the first retained coding exon; translate in one frame only | [4b audit](2026-11-10-phase-audit.md) |
| 2026-11-10 | **Remove the three-frame fallback and Tier 3** (switch `three_frame_fallback`, audit only, never released) | [4b audit §5.3, §7](2026-11-10-phase-audit.md#53-effect-of-the-three-frame-fallback-fallback--one_frame) |
| 2026-11-10 | **Keep a frame for slices lacking the start codon** (switch `frame_for_start_less_slice`, audit only, never released); add an annotated-start rule in PR 4c | [4b audit §5.4, §7](2026-11-10-phase-audit.md#54-slices-that-lack-the-start-codon) |
| 2026-11-10 | Alternative C-terminal isoforms go to a separate `psq_T1_alt_cterm.fasta` | [C-terminal audit §7](2026-11-10-alt-cterm-tiering-audit.md) |
| 2026-11-10 | PR 4c: annotated sites at TSL ≤ 2; one file per terminus (`psq_T1_alt_cterm`, `psq_T1_alt_nterm`, `psq_T1_alt_both`, each with `_orphan`); one fate `SUCCESS TIER 1 ALT TERMINUS`; frameshift test waived only for a slice ending at an internal annotated stop | [termini audit §7](2026-11-10-alt-termini-audit.md#7-decisions-needed) |
| 2026-11-10 | PR 4b2: fix the frameshift test (coding lengths) before PR 4c | [termini audit §6](2026-11-10-alt-termini-audit.md#6-regression-found-the-frameshift-test-ignores-per-slice-frames-proposed-pr-4b2) |

## Current tiers (after PR 4b)

Before tiering, each junction is trimmed to the gene's annotated coding
region, and exons entirely in the UTR are removed. Junctions whose event
lies outside that region are skipped. Each slice is then translated once,
in its own annotated frame; a slice with no annotated frame is not
translated.

| Tier | Meaning |
|---|---|
| **T1** | Both slices translate with no premature stop codon, and their coding lengths (from each slice's own frame) differ by a multiple of 3. |
| **T2** | As T1, but with a frameshift. |
| **T4** | One slice translates. The other is translated up to its first stop codon; that product is joined N-terminally if it is ≥ `ptc_threshold` × the translated slice's length, otherwise it is written as an orphan. |

T3 no longer exists. The full definitions are in
[REVIEW.md §3](../REVIEW.md#tier-classification-current-after-pr-4b).

## Golden-test history

`tests/test_pipeline.py::test_chr15_golden` pins every junction's fate and
every output file's contents. Its values at each stage:

| Stage | Commit | FAIL | T4 | Skip: low count | Skip: outside CDS | T1 | T2 | T3 | `psq_T1` records |
|---|---|---|---|---|---|---|---|---|---|
| v0.3.5a code (PR 1 golden) | `3ed17b0` | 360 | 236 | 297 | – | 185 | 25 | 2 | 152 |
| PR 4: stop-codon bounds | `da72a0a` | 330 | 261 | 297 | – | 189 | 26 | 2 | 164 |
| PR 4a: trim to the coding region | `89e1a82` | 240 | 259 | 297 | 90 | 190 | 27 | 2 | 169 |
| PR 4a2: RI slices | `bc7d894` | 238 | 258 | 297 | 92 | 189 | 28 | 3 | 169 |
| PR 4b candidate: per-slice frames, fallback kept | `9e27815` | 165 | 289 | 297 | 92 | 214 | 45 | 3 | 170 |
| PR 4b final: fallback and T3 removed | `2679edb` | 157 | 302 | 297 | 92 | 212 | 45 | – | 170 |
| PR 4b2: frameshift test from coding lengths | PR 4b2 commit | 157 | 302 | 297 | 92 | 227 | 30 | – | 171 |

`psq_canonical.fasta` has 268 records at every stage. At no stage is a
`psq_T1.fasta` sequence matching an annotated protein lost.

## New defects found during the audits

- **HIGH-11**, found in the PR 4a audit: RI slices contained both flanking
  exons twice. Fixed in PR 4a2.
- **A5SS/A3SS on the `−` strand were trimmed at the wrong bound**, found in
  the PR 4a audit. Fixed in PR 4a.
- **The frameshift test ignored per-slice frames after PR 4b,** found in
  the PR 4c audit (§6). It misclassifies T1/T2 for 17 chr15 junctions.
  Fixed in PR 4b2.
- **The T4 length threshold was mis-documented.** It decides whether the
  truncated product is joined to the canonical, not whether it is written.
  Corrected in the README and REVIEW.md with PR 4b.
