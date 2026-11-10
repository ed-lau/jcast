# Audit: reading frame from the first retained coding exon, and the three-frame fallback (PR 4b)

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Audit of the PR 4b candidate (committed with its two audit switches so the four-variant comparison is reproducible; see §8). Both recommendations in §7 were approved on 2026-11-10 and implemented in the final PR 4b commit. |
| **Baseline** | PR 4a2, commit `bc7d894` |
| **Code audited** | `Junction.get_translated_phase` (rewritten), `Junction.slice_exons`, `Sequence.translate` / `translate_forced` (per-slice frames), the T3 path in `main._translate_one` |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions, Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/phase_audit/audit_phase.py` (≈ 5 min, no network) against the candidate working tree. See [§8](#8-reproducing-this-audit). |
| **Evidence** | Four-variant comparison (candidate commit `9e27815`): [phase_audit/summary.txt](phase_audit/summary.txt) · [phase_audit/junctions.tsv](phase_audit/junctions.tsv) (one row per junction, all variants) · [phase_audit/fasta_diff.tsv](phase_audit/fasta_diff.tsv) (every output sequence differing from the baseline, per variant, with exact annotation matches). Final verification (§9): [phase_audit/final/summary.txt](phase_audit/final/summary.txt) · [phase_audit/final/junctions.tsv](phase_audit/final/junctions.tsv) |

Every number in this report is printed in `summary.txt` or can be read
directly from the two TSV files.

---

## 1. Summary

**Problem.** `get_translated_phase` took the reading frame from one fixed
exon, chosen by junction type and strand, for example the anchor on `+`. It
required that exon to match a CDS row *exactly at both ends*. This fails
whenever trimming has cut or removed that exon (PR 4/4a), or rMATS's exon
boundaries differ from the CDS row. It also assumes both slices start in the
same exon.

With no phase, JCAST fell back to guessing all three frames. That fallback
keeps the *last* frame tried (MED-2), so it mostly ends in FAILURE.

**Candidate.** Each slice gets its own frame, read at its 5′ end (the start
of its first retained exon) from the protein-coding CDS rows covering that
position:

- first those at the bounds' TSL threshold, then any TSL;
- the majority frame wins if transcripts disagree.

The frame is computed from the CDS row's GTF `frame` and the distance to the
row's 5′ end, so no exact boundary match is needed.

**Results on chr15** (candidate without the three-frame fallback):

- **The new lookup agrees wherever the old one worked.** Both found a frame
  for 532 junctions, and the frames are identical in all 532.
- **104 more junctions get a frame:**
  - 75 with the same frame for both slices;
  - 17 whose slices start in different exons and need *different* frames,
    which the old single phase could not represent;
  - 12 where only one slice gets a frame.
- **FAILURE falls from 238 to 157:**
  - FAILURE → T1: 25 junctions;
  - FAILURE → T2: 17;
  - FAILURE → T4: 41.
- **No sequence matching an annotated protein is lost.**
  - Annotated protein-coding matches in `psq_T1.fasta` go from 67 to 68
    (GTF2A2-208), and in `psq_T2.fasta` from 2 to 3 (LINS1-207, built in
    per-slice frames 2 and 0).
  - `psq_canonical.fasta` is unchanged.

**The three-frame fallback is nearly unused once frames are read correctly.**
Removing it changes 15 junctions:

| Change | Junctions | Detail |
|---|---|---|
| T3 → T4 | 3 | |
| FAILURE → T4 | 10 | Junctions where only one slice has a frame, which the fallback made fail |
| T1 → FAILURE | 2 | GOLGA2P10, a pseudogene with no CDS, translated in a guessed frame |

No annotated sequence is lost. **Recommended: remove the fallback and the T3
tier** (decision 1, §7).

**New finding: slices without the start codon.** In 42 junctions, the gene's
start codon lies in an alternative exon, so only one slice contains it.

- The other slice is read in-frame from mid-CDS.
- 32 of these junctions become T1 or T2 under the candidate.
- Usually the start-less slice either re-stitches into exactly the canonical
  protein (ACTC1, SPG21) or becomes an orphan fragment.
- In a few it is a real isoform whose own transcript starts at a later start
  codon (GTF2A2-208, LINS1-207).
- Withholding a frame from such slices was tested and loses those real
  isoforms and one canonical protein. **Recommended: keep reading a frame for
  them, and handle alternative N-termini with an annotated-start rule
  alongside PR 4c** (decision 2, §7).

---

## 2. Background

The old lookup (`get_translated_phase` before this PR) worked as follows:

- **Exon choice.** It picked `anc` (+) or `down` (−) for MXE/SE/RI, `alt1`
  for A5SS, and `anc` for A3SS.
- **Exact match.** It queried protein-coding CDS rows with
  `start == exon start` and `end == exon end`.
- **Single phase.** It set one `phase` for both slices, using the first
  matching row's frame. TSL was not considered.

`Sequence.translate(use_phase=True)` used that phase if it existed, and
otherwise ran the three-frame loop. When translation in the phase failed for
either slice, `_translate_one` tried the three-frame loop again, producing T3
("GTF phase mismatch").

Under PR 4a2, 184 of the 716 non-skipped junctions (26 %) had no phase.

---

## 3. Candidate

1. **`Junction.slice_exons()`** returns each slice's retained exons. It is
   composed as in `make_slice_localgenome` and is shared with
   `event_outside_cds`.
2. **`Junction.get_translated_phase()`**, for each slice:
   - takes its 5′ end: the lowest start on `+`, the highest end on `−`;
   - collects every protein-coding CDS row covering that position, first at
     `transcript_support_level <= params.tsl_threshold`, then at any TSL;
   - computes the frame at that position:
     `(start + frame − pos) mod 3` on `+`, `(frame − (end − pos)) mod 3` on
     `−`;
   - takes the majority frame (ties go to the lower frame) and sets
     `phase_ambiguous` when rows disagree.

   It sets `slice_phases = (p1, p2)` and keeps `phase` only when the two
   are equal.
3. **`Sequence.translate(use_phase=True)`** translates each slice in its own
   frame. `translate_forced` (T4) uses the slice's own frame. FASTA headers
   report each slice's own frame.
4. **Temporary switches**, used only for this audit and to be removed in the
   final PR:
   - `params.three_frame_fallback` — when False, a slice without a frame is
     not translated and the T3 path is skipped;
   - `params.frame_for_start_less_slice` — when False, a slice lacking the
     start codon that the other slice contains gets no frame.

---

## 4. Method

[phase_audit/audit_phase.py](phase_audit/audit_phase.py) runs
`main._translate_one` on all 1,105 junctions in four configurations, each in
its own process:

| Variant | Code | Three-frame fallback | Frame for start-less slices |
|---|---|---|---|
| `baseline` | git worktree of `bc7d894` (PR 4a2) | yes (old behaviour) | – |
| `fallback` | candidate | yes | yes |
| `one_frame` | candidate | **no** | yes |
| `no_startless` | candidate | no | **no** |

The settings are those of `tests/test_pipeline.py`. Per junction, the script
records:

- the fate, the phases, and whether the frame was ambiguous;
- whether three-frame translation was invoked;
- whether each slice contains the gene's start codon (the most upstream one,
  as used for `tx0`/`tx1`).

Output sequences are compared with the baseline and matched exactly against
all proteins translated from the GTF, with transcript names and biotypes.

---

## 5. Results

### 5.1 Phase availability (716 junctions not skipped)

| Candidate frames | Baseline had a phase | Junctions |
|---|---|---|
| both slices, same frame | yes | 532 |
| both slices, same frame | no | 75 |
| both slices, different frames | no | 17 |
| one slice only | no | 12 |
| neither slice | no | 80 |

- **Agreement:** in all 532 junctions where both versions found a frame, the
  frames are identical.
- **Ambiguity:** for 2 junctions, transcripts disagree on the frame and the
  majority was used.
- **The 80 with no frame:**
  - 41 are in genes with no protein-coding CDS at all (lncRNAs and
    pseudogenes, e.g. LINC01197, HERC2P9, GOLGA2P10);
  - 39 are in coding genes whose CDS rows do not cover the slice's start.
- **The 12 with one frame:** the other slice starts in sequence that no CDS
  covers.

### 5.2 Fates

| Fate | Baseline (PR 4a2) | `fallback` | `one_frame` | `no_startless` |
|---|---|---|---|---|
| FAILURE | 238 | 165 | **157** | 160 |
| PARTIAL TIER 4 | 258 | 289 | **302** | 332 |
| SKIPPED (low coverage) | 297 | 297 | 297 | 297 |
| SKIPPED (outside coding region) | 92 | 92 | 92 | 92 |
| SUCCESS TIER 1 | 189 | 214 | **212** | 192 |
| SUCCESS TIER 2 | 28 | 45 | **45** | 32 |
| SUCCESS TIER 3 | 3 | 3 | **0** | 0 |

Baseline → `one_frame` changes 88 junctions:

| Transition | Junctions |
|---|---|
| FAILURE → T1 | 25 |
| FAILURE → T2 | 17 |
| FAILURE → T4 | 41 |
| T3 → T4 | 3 |
| T1 → FAILURE | 2 |

Nothing else moves down.

Where the gene's start codon lies, for the 83 FAILURE → T1/T2/T4 junctions:

| Start codon position | FAIL → T1 | FAIL → T2 | FAIL → T4 |
|---|---|---|---|
| upstream of the junction (in neither slice) | 5 | 5 | 29 |
| in both slices | – | – | 5 |
| **in one slice only** | **20** | **12** | 7 |

### 5.3 Effect of the three-frame fallback (`fallback` → `one_frame`)

Only 15 junctions differ:

| Junctions | `fallback` → `one_frame` | Explanation |
|---|---|---|
| SNRPN MXE 3732/3734/3736/3738, SNRPN SE 43853, ANXA2 SE 12435/12439/12442/12444/12446 | FAILURE → T4 | Only one slice has a frame. The fallback guessed frames for both slices together, keeping the last frame tried (MED-2), and failed. Without it, the framed slice is translated. |
| USP3 SE 6181, ANXA2 SE 12445, MAN2C1 RI 953 | T3 → T4 | One slice hits a stop codon in its annotated frame. The fallback found another frame in which both slices read through. That is a guess, and T3 output is not used downstream. |
| GOLGA2P10 SE 5796/5797 | T1 → FAILURE | A pseudogene with no CDS, translated in a guessed frame. |

**Removing the fallback loses no sequence matching an annotated protein.**
Compared with the baseline, `one_frame` removes 10 records, all without an
annotated match:

- 6 `psq_T3_orphan` records;
- 3 `psq_T1_orphan` records;
- 1 `psq_T2_orphan` record.

### 5.4 Slices that lack the start codon

In 42 non-skipped junctions, the gene's start codon (the most upstream
annotated, TSL ≤ 2) lies in an exon only one slice contains. Typically it is
in the skipped exon, and the exon upstream of it is 5′ UTR and removed by
PR 4a. The other slice then begins mid-CDS, with no start codon.

**What the candidate does.** It reads that slice in the annotated frame of
its first exon. Then:

- **SE pattern (ACTC1 SE 6582, SPG21 SE 47956):** the start-less slice
  (anchor only) is a substring of the canonical protein. Stitching restores
  the canonical N-terminus, recreating exactly the canonical. The junction
  counts as T1 and writes only the canonical, which is already present.
- **Unstitchable start-less slices** become orphan N-terminal fragments.
- **A real isoform** can result when that slice contains its own
  transcript's start codon, downstream of the gene-wide first one:
  - GTF2A2 SE 11951: T1 record = GTF2A2-208;
  - LINS1 MXE 449: T2 record = LINS1-207, and LINS1's canonical.

**Does the start-less slice contain another annotated start codon?**

| | Junctions |
|---|---|
| yes, of a TSL ≤ 2 protein-coding transcript | 10 |
| only of TSL > 2 transcripts | 9 |
| no | 23 |

**The `no_startless` variant** gives such slices no frame. Compared with
`one_frame`, it turns 33 T1/T2 junctions into T4 and 3 T4 into FAILURE, and
removes:

- GTF2A2-208 from `psq_T1`;
- LINS1-207 from `psq_T2`;
- LINS1's canonical protein from `psq_canonical`.

So the blunt rule removes real isoforms.

**Effect on output.** On chr15, the T1/T2 fates of these junctions add
almost no records:

- `psq_T1.fasta`: 169 → 170;
- `psq_T2.fasta`: 22 → 23.

Their T1/T2 counts are therefore inflated, but the FASTA content is not
contaminated.

### 5.5 Output files

| File | Baseline | `fallback` | `one_frame` | `no_startless` |
|---|---|---|---|---|
| `psq_canonical.fasta` | 268 | 268 | 268 | 267 |
| `psq_T1.fasta` | 169 | 170 | 170 | 169 |
| `psq_T1_orphan.fasta` | 35 | 39 | 36 | 32 |
| `psq_T2.fasta` | 22 | 23 | 23 | 22 |
| `psq_T2_orphan.fasta` | 14 | 16 | 16 | 14 |
| `psq_T3_orphan.fasta` | 6 | 6 | 0 | 0 |
| `psq_T4.fasta` | 213 | 232 | 235 | 238 |
| `psq_T4_orphan.fasta` | 74 | 82 | 82 | 87 |

**Gained in `one_frame`, matching an annotated protein:**

| File | Protein-coding | NMD |
|---|---|---|
| `psq_T1` | 1 | – |
| `psq_T1_orphan` | 1 | – |
| `psq_T2` | 1 | – |
| `psq_T4` | 5 | 4 |
| `psq_T4_orphan` | – | 3 |

Lost in `one_frame`: 10 records, none matching an annotated protein.

### 5.6 Cases carried over from earlier audits

| Junction | PR 4a2 | Candidate (`one_frame`) | Note |
|---|---|---|---|
| ANXA2 SE 12446 | FAILURE | T4 | Only the skipping slice has a frame (2) |
| GTF2A2 SE 11945 | FAILURE | T2 | Per-slice frames 0 and 2 |
| ZFYVE19 SE 6502 | FAILURE | FAILURE | Both slices have frame 0 at the first start codon, yet both hit stop codons. Its anchor splicing belongs to a transcript starting at a later start codon. That is the gene-wide envelope issue ([tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions) item 2), not a phase problem. |
| MAN2C1 RI 953 | T3 | T4 | §5.3 |

---

## 6. What the final PR contains

Both recommendations in §7 were approved; §9 records the outcome.

1. Per-slice frames as in §3.
2. **Remove** the three-frame loop in `Sequence.translate`, the T3 path in
   `_translate_one`, `fates.TIER3`, and the `psq_T3*.fasta` outputs.
   - MED-2 is closed by removal.
   - Update the README's tier descriptions.
   - Remove both temporary switches.
3. **Regression tests,** landed first:
   - GTF2A2 SE 11945: per-slice frames (0, 2), fate T2;
   - LINS1 MXE 449: frames (2, 0), and the T2 record equals the LINS1-207
     protein;
   - ANXA2 SE 12446: only the skipping slice framed, fate T4;
   - GOLGA2P10 SE 5796: no CDS, so FAILURE (no guessed frame);
   - a frame-agreement test on an existing T1 junction.
4. **Golden update** (`one_frame` column of §5.2/§5.5):
   - fates: FAIL 157, T4 302, SKIP 297 + 92, T1 212, T2 45; T3 no longer
     exists;
   - new digests.
5. **Re-run this audit as baseline vs final.** It must reproduce the
   `one_frame` column exactly.

---

## 7. Decisions needed

**1. Three-frame fallback and T3: remove (recommended), or keep.**

- **Evidence for removal (§5.3).** Once frames are read correctly, the
  fallback changes only 15 junctions. Each change either is a guessed frame
  (T3; pseudogene T1s) or hides a translatable slice (10 FAILURE → T4). No
  annotated sequence depends on it, and T3 output is not used downstream.
- **Keeping it** would also require fixing MED-2.

**2. Slices lacking the start codon that the other slice contains: keep
reading a frame (recommended), or withhold it.**

- **Withholding** removes real isoforms (§5.4).
- **Keeping** inflates T1/T2 counts for about 32 junctions whose start-less
  slice is a fragment or recreates the canonical, but adds no wrong records
  on chr15.
- **The principled fix is an annotated-start rule,** the mirror of PR 4c's
  annotated-stop rule. A start-less slice that contains another annotated
  start codon (10 of 42) would be translated from that start and joined
  C-terminally only. Otherwise it would be marked as having no start.
  Recommended: fold this into PR 4c as "alternative termini" (both ends),
  with its own output file mirroring `psq_T1_alt_cterm.fasta`.

---

## 8. Reproducing this audit

```bash
# four-variant comparison: needs the candidate's audit switches, so check out 9e27815 first
git checkout 9e27815 && python reports/phase_audit/audit_phase.py

# final verification: baseline bc7d894 vs the PR 4b code, checked against the one_frame variant
python reports/phase_audit/audit_phase.py --final
```

The candidate's two temporary switches exist only at commit `9e27815`,
where their defaults reproduce the old behaviour and the golden test pins
the `fallback` variant. The final PR 4b commit removes the switches.
`--final` re-runs baseline vs final, and verifies that every junction's
fate and every changed output sequence equal the audited `one_frame`
variant.

The script creates and removes its own git worktree for the baseline, stubs
out UniProt, and rewrites the three evidence files in `reports/phase_audit/`.
The inputs are the same as for the
[tx1 audit, §9](2026-11-03-tx1-translation-bounds-audit.md#9-reproducing-this-audit).

---

## 9. Final implementation and verification

**Decisions (approved 2026-11-10).**

1. **Remove the three-frame fallback and Tier 3.** Implemented.
2. **Keep reading a frame for slices that lack the start codon.** Implemented
   as the default; there is no switch. The annotated-start rule is scheduled
   in PR 4c.

**Commits.**

| Commit | Content |
|---|---|
| `9e27815` | Candidate: per-slice frames plus the two audit switches, defaulting to the old behaviour. Golden test = `fallback` variant. |
| `2679edb` (final PR 4b) | Removes `Sequence.translate`'s three-frame loop and its `use_phase` argument, the T3 path in `_translate_one`, `fates.TIER3`, both switches, and the T3 entry in the GUI tier counters. Golden test = `one_frame` variant. |

**Regression tests,** each confirmed to fail on the code before the change:

| Test | Fails on | Asserts |
|---|---|---|
| `test_phase_unchanged_where_it_was_found` | `bc7d894` | MFGE8 SE 5739: frames (2, 2), equal to the old phase |
| `test_per_slice_phase` | `bc7d894` | GTF2A2 SE 11945: frames (0, 2); fate T2 |
| `test_per_slice_phase_reconstructs_annotated_isoform` | `bc7d894` | LINS1 MXE 449: frames (2, 0); fate T2; the T2 file contains the LINS1-207 protein, translated from the GTF |
| `test_slice_without_frame_is_not_guessed` | `9e27815` | ANXA2 SE 12446: frames (None, 2); fate T4 |
| `test_no_coding_annotation_is_not_translated` | `9e27815` | GOLGA2P10 SE 5796 (pseudogene): no frames; FAILURE |
| `test_annotated_frame_is_not_overridden` | `9e27815` | USP3 SE 6181 (formerly T3): T4; no `psq_T3*` file |

**Verification**
([final/summary.txt](phase_audit/final/summary.txt)).
`audit_phase.py --final` compares `bc7d894` with the final code:

- **Fates:** for 0 of 1,105 junctions does the final fate differ from the
  audited `one_frame` fate.
- **Output:** 49 sequences differ from the baseline, the identical set to
  the audited `one_frame` variant.
- **Fate counts:**

  | Fate | Junctions |
  |---|---|
  | FAILURE | 157 |
  | PARTIAL TIER 4 | 302 |
  | SKIPPED (low coverage) | 297 |
  | SKIPPED (outside coding region) | 92 |
  | SUCCESS TIER 1 | 212 |
  | SUCCESS TIER 2 | 45 |

- **Files:** `psq_T1.fasta` has 170 records and `psq_canonical.fasta` 268;
  no `psq_T3*` file is produced.

The full test suite (34 tests, including the network tests) passes.

**Current tier definitions** are in
[REVIEW.md §3](../REVIEW.md#tier-classification-current-after-pr-4c), the
[README](../README.md#output-files) and the
[CHANGELOG](../CHANGELOG.md#tiers-in-v040-summary-of-prs-44b2).
