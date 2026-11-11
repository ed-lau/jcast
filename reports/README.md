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
| [Frameshift test from coding lengths](2026-11-10-frameshift-audit.md) | 4b2 | Implemented in `a0f3338` | `a0f3338` (`python reports/frameshift_audit/audit_frameshift.py`) |
| [Alternative C-termini (first evaluation)](2026-11-10-alt-cterm-tiering-audit.md) | 4c | Superseded by the next row. Measured on PR 4 alone. | `a555446` (`python reports/cterm_audit/audit_cterm.py`) |
| [Alternative C- and N-termini on the current pipeline](2026-11-10-alt-termini-audit.md) | 4c (and 4b2) | **Implemented in PR 4c** (§10), after PR 4b2, which fixes the regression it found (§6). | evaluation: `919f6c7` (`audit_termini.py`); expected after 4b2: `a0f3338` (`audit_termini.py --out …/expected_after_4b2`); verification: `e11cc86` (`verify_4c.py`) |
| [Soft-masked and ambiguous genome bases](2026-11-11-masked-bases-audit.md) | 3 | **Implemented in `fcdc0ee`** (§8); decisions approved 2026-11-11 (§7) | audit: `819cc53`; verification: `fcdc0ee` (`python reports/mask_audit/audit_mask.py --softmasked …dna_sm.chromosome.15.fa.gz [--verify]`); genome survey with `survey_genome.py` (GENCODE GRCh38 + v50 GTF, not bundled) |
| [Joining slices to the canonical; anchor length](2026-11-11-stitch-anchor-audit.md) | 5 | **Implemented in `2142711`** (§9); decisions approved 2026-11-11 (§8) | audit: `887a98a` (`python reports/stitch_audit/audit_stitch.py`); verification: `2142711` (`verify_5.py`); proteome survey with `survey_anchors.py` (GENCODE v50, not bundled) |
| [A supplied UniProt FASTA; joining by coordinates](2026-11-11-uniprot-fasta-evaluation.md) | 14 | **Implemented in `0f3b409`** (§8); decisions approved 2026-11-11 (§7) | evaluation: `8b0c024` (`python reports/uniprot_eval/evaluate_uniprot_fasta.py --uniprot <Swiss-Prot FASTA>`); verification: `0f3b409` (`verify_14.py --uniprot <Swiss-Prot FASTA>`); both need the local UniProt cache |

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
| 2026-11-11 | PR 3: uppercase genome sequence at retrieval; skip a slice with any base outside `ACGT` (as for N); keep `--mask` but deprecate it (no effect; remove in the next release); catch an exception per junction, record it and continue (in PR 7) | [masked-bases audit §7](2026-11-11-masked-bases-audit.md#7-decisions-requested) |
| 2026-11-11 | PR 5: join a peptide only at an end that is not already complete (stop codon, gene start codon); use an anchor only if it occurs once, start before end; anchor length 7; terminus comparison stays at 10 | [stitch audit §8](2026-11-11-stitch-anchor-audit.md#8-decisions-requested) |
| 2026-11-11 | PR 14: add `--uniprot FASTA` (canonical by sequence, then gene name; GTF canonical if neither) and a built-in Swiss-Prot merge of `psq_T1` + `psq_T1_alt_*` without `tr|` records; keep the API path without `--uniprot`; defer joining by coordinates | [UniProt evaluation §7](2026-11-11-uniprot-fasta-evaluation.md#7-decisions-2026-11-11) |

## Current tiers (after PR 4c)

Before tiering, each junction is trimmed to the gene's annotated coding
region, and exons entirely in the UTR are removed. Junctions whose event
lies outside that region are skipped. Each slice is then translated once,
in its own annotated frame; a slice with no annotated frame is not
translated.

| Tier | Meaning |
|---|---|
| **T1** | Both slices translate with no premature stop codon, and their coding lengths (from each slice's own frame) differ by a multiple of 3. |
| **T1 ALT TERMINUS** | At least one slice ends at an annotated alternative stop codon (C-anchored, joined at its N-terminus only) or starts at an annotated alternative start codon (N-anchored, joined at its C-terminus only). Files: `psq_T1_alt_cterm`, `psq_T1_alt_nterm`, `psq_T1_alt_both` (each with `_orphan`). |
| **T2** | As T1, but with a frameshift. |
| **T4** | One slice translates. The other is translated up to its first stop codon; that product is joined N-terminally if it is ≥ `ptc_threshold` × the translated slice's length, otherwise it is written as an orphan. |

T3 no longer exists. The full definitions are in
[REVIEW.md §3](../REVIEW.md#tier-classification-current-after-pr-4c).

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
| PR 4b2: frameshift test from coding lengths | `a0f3338` | 157 | 302 | 297 | 92 | 227 | 30 | – | 171 |
| PR 4c: alternative termini (+ T1 ALT 25) | `e11cc86` | 153 | 289 | 297 | 92 | 219 | 30 | – | 173 |
| PR 5: joining to the canonical (+ T1 ALT 24) | `2142711` | 153 | 289 | 297 | 92 | 220 | 30 | – | 179 |

`psq_canonical.fasta` has 268 records at every stage. Up to PR 4c, no
`psq_T1.fasta` sequence matching an annotated protein is lost at any stage.
PR 5 removes one: CD276-201. It came only from a join at the wrong copy of
a repeated anchor
([stitch audit §9](2026-11-11-stitch-anchor-audit.md#9-implementation-and-verification-pr-5)).

PR 14 leaves this test unchanged: it does not use `--uniprot`. The new
`tests/test_uniprot.py::test_chr15_uniprot` pins the same run with
`--uniprot` on a Swiss-Prot 2026_04 subset: `psq_T1` 186
records, `psq_canonical` 275, and `jcast_uniprot_merged.fasta` 837 (711
UniProt entries + 126 JCAST records)
([UniProt evaluation §8](2026-11-11-uniprot-fasta-evaluation.md#8-implementation-and-verification-pr-14)).

## New defects found during the audits

- **HIGH-11**, found in the PR 4a audit: RI slices contained both flanking
  exons twice. Fixed in PR 4a2.
- **A5SS/A3SS on the `−` strand were trimmed at the wrong bound**, found in
  the PR 4a audit. Fixed in PR 4a.
- **The frameshift test ignored per-slice frames after PR 4b,** found in
  the PR 4c audit (§6). It misclassifies T1/T2 for 17 chr15 junctions.
  Fixed in PR 4b2.
- **`-c` does not write the canonical for genes whose junctions all
  succeed** (REVIEW.md MED-14; 2 chr15 genes), found while verifying
  PR 4c. Scheduled for PR 8.
- **Hard-masked genomes silently drop junctions**, found in the PR 3 audit
  (§4.2). The README recommended them as a workaround for soft-masking;
  on chr15, 206 of 1,105 junctions (20 successes) are skipped as masked.
- **One exception aborts the whole run**, found in the PR 3 audit: the
  per-junction loop catches nothing. Approved for PR 7.
- **Joins at an already-complete end returned the canonical**, found in
  the PR 5 audit. A T4 peptide ending at its stop codon was joined
  C-terminally; a peptide starting at the gene's start codon was joined
  N-terminally. 10 of 12 wrong joins on chr15 gave exactly the canonical.
  Fixed in PR 5.
- **The T4 length threshold was mis-documented.** It decides whether the
  truncated product is joined to the canonical, not whether it is written.
  Corrected in the README and REVIEW.md with PR 4b.
