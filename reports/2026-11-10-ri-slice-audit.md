# Audit: retained-intron slices contain their flanking exons twice (PR 4a2)

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Audit of the PR 4a2 implementation, run before it is committed. Baseline is PR 4a (commit `89e1a82`); candidate is the PR 4a2 working tree. Implemented in `bc7d894`. **Later changes:** see the [reports index](README.md); Tier 3 was removed in PR 4b, so T3 figures here are historical. |
| **Defect** | REVIEW.md HIGH-11, found during the [PR 4a audit](2026-11-10-utr-exon-trimming-audit.md#6-new-finding-ri-slice-1-duplicates-its-flanking-exons-reviewmd-high-11) |
| **Code audited** | `Sequence.make_slice_localgenome` ([jcast/sequences.py](../jcast/sequences.py)); `Junction.event_outside_cds` ([jcast/junctions.py](../jcast/junctions.py)) |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions (88 RI), Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/ri_audit/audit_ri.py` (≈ 2.5 min, no network). See [§7](#7-reproducing-this-audit). |
| **Evidence** | [ri_audit/summary.txt](ri_audit/summary.txt) · [ri_audit/junctions.tsv](ri_audit/junctions.tsv) (one row per RI junction) · [ri_audit/fasta_diff.tsv](ri_audit/fasta_diff.tsv) (every changed output sequence, with exact annotation matches) |

Every number in this report is printed in `summary.txt` or can be read
directly from the two TSV files.

---

## 1. Summary

**The defect.**

- rMATS reports a retained intron's *retained exon* (`riExonStart_0base`..
  `riExonEnd`) as one span. That span runs from the upstream exon's start to
  the downstream exon's end. This holds for all 88 RI rows on chr15.
- JCAST maps this span to `alt1`.
- `make_slice_localgenome` treated RI like SE and MXE, building the
  intron-retaining slice as `anc + alt1 + down`.

So both flanking exons appeared **twice**. In UNC45A RI 916, slice 1 was
37 + 515 + 324 = 876 nt; the retained transcript is 515 nt. The defect has
been present since at least v0.2.6 (2020-07). The earlier code used the same
concatenation for every junction type.

**The fix.**

- For RI, slice 1 is the retained exon alone, and slice 2 (spliced) is
  `anc + down`, as before.
- `event_outside_cds` uses the same composition.

**Output on chr15.**

- No non-RI junction changes fate, coordinates or phase.
- 5 of 88 RI junctions change fate. Two are artefacts now correctly skipped.
- The content of RI output changes far more: 58 RI sequences are removed and
  54 added, mostly in `psq_T4.fasta`.
  - 29 of the 56 removed retained-intron records contain a repeated 8-aa
    peptide, the signature of the duplication.
  - 2 of the 51 added ones do, consistent with natural repeats.
- `psq_T1.fasta` and `psq_canonical.fasta` are unchanged.

**What this means for past results.**

- Every retained-intron (RI slice 1) sequence JCAST has produced was
  translated from a sequence with duplicated exons.
- Wherever the anchor exon's length is not a multiple of 3, the duplication
  shifted the reading frame after the anchor. The output was then mostly
  out-of-frame sequence up to a spurious stop.
- On this data the affected records were almost all in `psq_T4.fasta` and
  `psq_T4_orphan.fasta`. Spliced RI slices (slice 2) were correct.

---

## 2. Evidence for the defect

- **Coordinates.** In all 88 RI rows of `tests/data/rmats/RI.MATS.JC.txt`,
  `riExonStart_0base == upstreamES` and `riExonEnd == downstreamEE`.
- **Slice length.** `tests/test_junction.py::test_retained_intron_slices`
  builds UNC45A RI 916:
  - shipped code: slice 1 is 876 nt, failing the test;
  - fixed code: slice 1 is exactly the 515-nt genomic retained exon, and
    slice 2 is anchor + downstream (37 + 324 nt).
- **Size of the excess.** Across the 75 analysed (non-skipped) RI junctions,
  slice 1 totals 81,472 nt under PR 4a and 59,365 nt under PR 4a2. The
  difference is the duplicated flanks.

---

## 3. Method

[ri_audit/audit_ri.py](ri_audit/audit_ri.py) reuses the harness from the
[PR 4a audit](utr_audit/audit_utr.py). It runs `main._translate_one` on all
1,105 junctions:

- once in a git worktree of `89e1a82` (PR 4a);
- once in the working tree (PR 4a2).

Each run has its own process, with the same settings as
`tests/test_pipeline.py`. The script then:

- checks that every non-RI junction is identical in fate, trimmed
  coordinates and phase;
- tabulates RI fates and slice lengths;
- diffs every output FASTA;
- compares each changed sequence with all proteins translated from the GTF,
  recording the matching transcripts and their biotypes.

**Limit of annotation-based validation.** Ensembl annotates most
retained-intron transcripts with the `retained_intron` biotype, which has no
CDS. Annotation therefore cannot confirm most RI translations. The evidence
for correctness is structural: the slice equals the genomic retained exon
(§2), and the duplication signature disappears (§4.3).

---

## 4. Results

### 4.1 Fates

| RI fate | PR 4a | PR 4a2 |
|---|---|---|
| FAILURE | 19 | 17 |
| PARTIAL TIER 4 | 57 | 56 |
| SKIPPED (low coverage) | 10 | 10 |
| SKIPPED (outside coding region) | 1 | 3 |
| SUCCESS TIER 1 | 1 | 0 |
| SUCCESS TIER 2 | 0 | 1 |
| SUCCESS TIER 3 | 0 | 1 |

| Junction | Strand | PR 4a → PR 4a2 | Explanation |
|---|---|---|---|
| FAH RI 943 | + | T1 → SKIPPED (outside CDS) | The anchor and the intron lie upstream of the start codon, which falls at the downstream exon's start. After trimming, the retained exon equals the downstream exon. The PR 4a "T1" was that exon translated twice (orphan slices of 53 and 26 aa: 53 = 2 × 26 + 1). Now correctly identified as an event outside the coding region. |
| FAM219B RI 963 | − | FAIL → SKIPPED (outside CDS) | Same configuration on the `−` strand. |
| FRMD5 RI 4475 | − | FAIL → T4 | The retained slice now hits its genuine in-intron stop, and the spliced slice translates. A typical retained-intron outcome. |
| FURIN RI 924 | + | T4 → T2 | The "intron" is a single nucleotide (anchor ends at 90,879,892, downstream starts at 90,879,894), so the retained slice is the spliced slice plus 1 nt. That is a frameshift, so T2. It is probably an rMATS alignment artefact, and the read count is 2. |
| MAN2C1 RI 953 | − | T4 → T3 | The retained slice no longer translates in the GTF frame, but both slices translate in another frame. T3 output is not used downstream of JCAST (REVIEW.md §5.1 decisions), and PR 4b will revisit it. |

**Why T4 dominates.** 56 of 75 analysed RI junctions remain T4. That is the
expected biology: a retained intron usually introduces an in-frame stop
(only 18 of 75 intron lengths are a multiple of 3). The T4 record is now the true
truncated retained-intron product, not a translation of duplicated sequence.

### 4.2 Output files

| File | PR 4a | PR 4a2 |
|---|---|---|
| `psq_canonical.fasta` | 268 | 268 |
| `psq_T1.fasta` | 169 | 169 |
| `psq_T1_orphan.fasta` | 37 | 35 |
| `psq_T2.fasta` | 22 | 22 |
| `psq_T2_orphan.fasta` | 13 | 14 |
| `psq_T3_orphan.fasta` | 4 | 6 |
| `psq_T4.fasta` | 217 | 213 |
| `psq_T4_orphan.fasta` | 75 | 74 |

All changed sequences are RI records. In detail:

- **Removed (58):** 48 from `psq_T4`, 8 from `psq_T4_orphan`, 2 from
  `psq_T1_orphan` (FAH). None matches an annotated protein.
- **Added (54):** 44 in `psq_T4`, 7 in `psq_T4_orphan`, 2 in
  `psq_T3_orphan`, 1 in `psq_T2_orphan`. One (in `psq_T4`) exactly matches
  an annotated `nonsense_mediated_decay` protein. None matched any
  annotation before.

### 4.3 The duplication signature

A duplicated exon whose length is a multiple of 3 repeats its peptide. Counting
records that contain any 8-aa peptide twice:

| Records | With a repeated 8-mer |
|---|---|
| Retained-intron records removed by the fix | 29 of 56 |
| Retained-intron records added by the fix | 2 of 51 |

No repeat is expected in the other 27 removed records in two situations,
which this audit does not separate:

- the anchor's length is not a multiple of 3, so the duplicate is read in
  another frame;
- translation stops before reaching the duplicate.

---

## 5. Regression tests

| Test | Asserts |
|---|---|
| `test_retained_intron_slices` | UNC45A RI 916: slice 1 is the 515-nt genomic retained exon; slice 2 is anchor + downstream. Fails on PR 4a (876 nt). |
| `test_retained_intron_in_utr_is_skipped` | FAH RI 943 is `SKIPPED_OUTSIDE_CDS`. On PR 4a it was T1, from the doubled exon. |
| `test_chr15_golden` | New fate counts: FAIL 238, T4 258, SKIP low coverage 297, SKIP outside CDS 92, T1 189, T2 28, T3 3. New digests for the 5 changed files. |

---

## 6. Follow-ups

1. **PR 4c script.** The C-terminal audit
   ([audit_cterm.py](cterm_audit/audit_cterm.py)) mirrors the old RI slice
   composition in `segments()`. Update it when that audit is re-run.
2. **PR 4b.** Taking the phase from the first retained coding exon will also
   govern RI. The fix here keeps both RI slices starting at the same base, so
   the phase is shared.
3. **1-nt "introns"** (FURIN RI 924) could be filtered as rMATS artefacts.
   This is low priority; noted only.

---

## 7. Reproducing this audit

From the repository root, with the package dependencies installed:

```bash
python reports/ri_audit/audit_ri.py               # baseline 89e1a82 vs the working tree
python reports/ri_audit/audit_ri.py --baseline X  # any other baseline commit
```

The script creates and removes its own git worktree for the baseline, stubs
out UniProt, and rewrites the three evidence files in `reports/ri_audit/`.
The inputs are the same as for the
[tx1 audit, §9](2026-11-03-tx1-translation-bounds-audit.md#9-reproducing-this-audit).
