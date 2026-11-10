# Audit: alternative C-termini are tiered as premature stops

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Pre-implementation evaluation. No production code changed. Proposed as PR 4c. It lands **after** PR 4 ([tx0/tx1 audit](2026-11-03-tx1-translation-bounds-audit.md)) and after PRs 4a (drop pure-UTR exons) and 4b (phase from the first retained coding exon). Those change the junction set analysed here, so the audit must be re-run before 4c is implemented and the figures below updated. |
| **Code audited** | JCAST v0.3.5a, commit `3ed17b0` (branch `v0.4.0-hardening`) with the PR 4 corrected bounds patched in at runtime. Tiering: [jcast/main.py:145-330](../jcast/main.py#L145-L330); translation: [jcast/helpers.py:157-243](../jcast/helpers.py#L157-L243); stitching: [jcast/sequences.py:465-525](../jcast/sequences.py#L465-L525) |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions, Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/cterm_audit/audit_cterm.py` (≈ 3 min, no network). See [§9](#9-reproducing-this-audit). |
| **Evidence** | [cterm_audit/summary.txt](cterm_audit/summary.txt) · [cterm_audit/slices.tsv](cterm_audit/slices.tsv) (one row per slice) · [cterm_audit/junctions.tsv](cterm_audit/junctions.tsv) (one row per junction) |

Every number in this report is printed in `summary.txt` or can be read directly
from the two TSV files. "Slice" means one of the two isoform sequences JCAST
builds per junction. Counts are of slices, before the FASTA writer removes
duplicate sequences.

---

## 1. Summary

**Question.** T1 output is joined ("stitched") to the canonical protein. Are
isoforms with an alternative C-terminus therefore excluded from T1? Would
marking them and joining them only at the N-terminus fix this?

**Answer.**

- **Stitching is not what excludes them.** Tier is decided before stitching,
  and an unstitched T1 slice still goes to `psq_T1_orphan.fasta`.
- **The PTC classification excludes them.** When a splice event makes
  translation end at a different stop codon, JCAST classifies that stop as a
  **premature termination codon (PTC)**. The junction then lands in T4 (one
  slice affected) or FAIL (both affected).
- **N-terminal-only joining already exists.** It is the T2/T4 fallback.
- **What is missing is a test that separates a real C-terminus from a PTC.**

**Rule evaluated.** A stop codon is a natural C-terminus if it sits at an
annotated `stop_codon` of a `protein_coding` transcript of the same gene.

- The annotated-stop rule discriminates sharply between real C-termini and
  PTCs (§5.1). Unannotated stops, including those in the slice's last exon,
  never yield an annotated protein (0 of 212).
- Applying the rule moves **46 junctions** out of T4 or FAIL. Of their 92
  slices, **67 %** reconstruct an annotated protein exactly, against **59 %**
  for today's T1.
- 49 slices go to the new file `psq_T1_alt_cterm.fasta`, 32 of them exact
  annotated proteins.
- **Nothing produced by today's T1 changes.** All 364 slices keep their file
  and sequence.
- No other junction changes fate.

**Your decision.** Alternative C-termini are written to a separate file,
`psq_T1_alt_cterm.fasta`, with unstitched ones in
`psq_T1_alt_cterm_orphan.fasta`. A slice is written there only if its
C-terminus differs from the canonical's. The specification is in §7.

---

## 2. Background: how tiering and stitching work

`_translate_one` ([main.py:145](../jcast/main.py#L145)) translates both slices
in the GTF frame with `make_pep(terminate=True)`. That call returns an empty
string if a stop codon occurs **anywhere except the last codon**. Tiers follow
from that result:

| Both slices non-empty? | Frameshift¹ | Tier | Stitching |
|---|---|---|---|
| yes | no | T1 | both ends must match the canonical; otherwise orphan file |
| yes | yes | T2 | both ends; fallback N-terminus only |
| no, one slice | – | T4 | PTC slice force-translated to its first stop, then N-terminus only, if ≥ 33 % of the other slice's length |
| no, neither | – | (T3 three-frame retry, else) FAIL | – |

¹ The slices' nucleotide lengths differ by a non-multiple of 3 ([sequences.py:66](../jcast/sequences.py#L66)).

Stitching ([sequences.py:465](../jcast/sequences.py#L465)) finds the slice's
first and last 10 aa in the canonical:

- **Both ends match:** the output is `canonical[:start] + slice + canonical[end:]`.
- **N-terminus only (`slice_has_ptc=True`):** the output is `canonical[:start] + slice`.

Consider a slice whose translation legitimately ends at a different stop than
the canonical's:

- **Stop before the last codon of the slice:** it is classed as a PTC, so the
  junction drops to T4 or FAIL.
- **Stop is the last codon:** it passes, but its last 10 aa are not in the
  canonical, so it becomes a T1 orphan.

In both cases the sequence is a correctly translated protein. Only the
labelling is wrong.

---

## 3. Options considered

| # | Option | Assessment |
|---|---|---|
| 1 | Join every PTC slice N-terminally and label it C-terminal | Already done for T4. It cannot tell real C-termini from PTCs: classes C + D below give 0 annotated proteins out of 212. |
| 2 | Treat a stop in the slice's last exon as natural (NMD-style rule) | Class C: 0 of 110 annotated. The slice's last exon is often not the transcript's last exon, so the rule has no basis here. |
| 3 | **Treat a stop at an annotated protein-coding stop codon as natural** | **Evaluated in this report.** It discriminates sharply (§5.1). |
| 4 | Assign each junction to specific transcripts and trim to that transcript's CDS | The proper long-term fix (REVIEW.md v1.0 splice-graph item; [tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions) item 2). Too large for v0.4. |

---

## 4. Method

Script: [cterm_audit/audit_cterm.py](cterm_audit/audit_cterm.py). It imports
`corrected_get_translated_region` from the tx1 audit.

1. **Actual fate.** Every junction is run through the shipped
   `main._translate_one` with the PR 4 bounds and the same settings as
   `tests/test_pipeline.py` (`-r 1`, `--mask NnXx`, `-q 0 1`, GTF canonical,
   UniProt stubbed). These fates match the tx1 audit's corrected column and
   the PR 4 golden values (FAIL 330, T4 261, SKIP 297, T1 189, T2 26, T3 2).
2. **Analysed set.** A junction is analysed if it passed the filters and has
   a GTF phase: **540** junctions. Those without a phase (259 FAIL, 7 T1,
   2 T2) take the three-frame route and are out of scope.
3. **Per slice.** The script rebuilds the slice nucleotides as
   `make_slice_localgenome` does and finds the first in-frame stop codon. It
   then maps that stop's 5′ base back to a genomic coordinate using the
   slice's exon segments.
4. **Stop class.** Each slice is classed by its first stop codon:

   | Class | First stop codon |
   |---|---|
   | open | none in the slice |
   | terminal | in the last codon (what `make_pep` already accepts) |
   | A | annotated `stop_codon` of a `protein_coding` transcript of the gene |
   | B | annotated only on non-protein-coding transcripts (all `nonsense_mediated_decay` here) |
   | C | unannotated, in the slice's 3′-most exon |
   | D | unannotated, upstream of a splice junction within the slice |

   Annotated stops come from all transcripts, with no TSL filter (see §8,
   item 4).
5. **Validation.**
   - **Truth set:** every CDS-bearing transcript in the GTF is translated
     (2,979 distinct proteins, any biotype).
   - **Match:** a stitched product counts only if it equals an annotated
     protein over its full length. The matching transcript names are recorded
     in the TSV.
   - **Labels:** "new" means stitched but matching no annotated protein;
     "unstitched" means the 10-aa anchor was not found in the canonical.
6. **Proposed routing (§7) is simulated.** A junction qualifies if both
   slices are open, terminal or A, and at least one is A. Junctions that are
   T1 today also qualify.
   - **C-anchored slices** (A, or terminal at an annotated protein-coding
     stop) are joined N-terminally only. They go to `psq_T1` if their last
     10 aa are the canonical's, otherwise to `psq_T1_alt_cterm`.
   - **Other slices** are stitched as today.
7. **Consistency check.** Some C-anchored slices already stitch at both ends
   today. For those, the N-only product must equal today's product, so the
   rule changes no existing output.

---

## 5. Results

### 5.1 The annotated-stop rule separates real C-termini from PTCs

Stitching each slice that has a stop codon N-terminally only gives:

| Stop class | Slices | Annotated protein-coding | Annotated, other biotype | New | Unstitched |
|---|---|---|---|---|---|
| A: annotated protein-coding stop | 67 | **37** | 3 | 22 | 5 |
| B: annotated only on NMD transcripts | 126 | 5 | **77** | 34 | 10 |
| C: unannotated, in the slice's last exon | 110 | 0 | 0 | 97 | 13 |
| D: unannotated, before a junction | 102 | 0 | 0 | 81 | 21 |

The slices with no stop (634) or a last-codon stop (41) are routed as today.

- **Class A (annotated protein-coding stop):** 37 of 62 stitched products
  (60 %) reproduce an annotated protein-coding protein exactly. These
  sequences are real C-terminal isoforms.
- **Class B (stop annotated only on NMD transcripts):** mostly reproduces the
  annotated NMD transcript's product. It should stay in T4.
- **Classes C and D (unannotated stops):** never reproduce an annotated
  protein, so the PTC label is appropriate for them.
- **All A slices are in T4 or FAIL today:** 35 are in FAIL junctions and 32
  in T4.
- **Where the stop sits:** 50 of the 67 A stops are in the slice's last exon
  and 17 are upstream of a junction, for example a stop inside an included
  exon that makes the downstream exon 3′ UTR. The 17 yield 8 annotated
  protein-coding products, so position alone would not have identified them.
- **12 A slices end at the canonical's own C-terminus.** The canonical's stop
  lies upstream of the gene-wide bound `tx1`, which is
  [tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions)
  item 2 ("one bound for all transcripts"). The rule rescues these too, and
  they go to `psq_T1`.

### 5.2 Fate changes

| Fate today (PR 4 bounds) | Proposed | Junctions |
|---|---|---|
| FAILURE | SUCCESS TIER 1 ALT CTERM | 11 |
| FAILURE | SUCCESS TIER 1 | 3 |
| PARTIAL TIER 4 | SUCCESS TIER 1 ALT CTERM | 32 |
| *all other analysed junctions* | *unchanged* | 494 |

The three FAILURE → T1 junctions are TMED3 SE 5845, FAM219B A3SS 1304 and
ADPGK A3SS 1551. In each, both C-anchored slices end at the canonical
C-terminus (the §5.1 "one bound for all transcripts" case).

No T1, T2 or T3 junction changes fate. No junction moves down.

### 5.3 Where the 92 slices of the 46 rescued junctions go

| Proposed file | Annotated protein-coding | Annotated NMD | New | Unstitched | Total |
|---|---|---|---|---|---|
| `psq_T1` | 30 | 0 | 5 | – | 35 |
| `psq_T1_alt_cterm` | **32** | 2 | 15 | – | 49 |
| `psq_T1_alt_cterm_orphan` | – | – | – | 4 | 4 |
| `psq_T1_orphan` | – | – | – | 4 | 4 |
| **Total** | **62 (67 %)** | 2 | 20 | 8 | 92 |

**Baseline.** The 364 slices of today's 182 T1 junctions:

| | Slices | Share |
|---|---|---|
| Annotated protein-coding | 215 | 59 % |
| New | 109 | |
| Unstitched | 40 | |

The rescued sequences are therefore at least as well supported by annotation
as existing T1 output.

**Unstitched slices.** All 8 are in PARP6 SE 13583/13584/13585 and HAUS2 A3SS
7146. In each, the partner slice, which has no stop codon, also fails to
stitch. The problem is with the canonical anchor in these genes, not with this
rule.

### 5.4 Existing T1 output is unchanged

- **Files and sequences.** All 364 slices of today's T1 junctions keep their
  file (324 in `psq_T1`, 40 in `psq_T1_orphan`) and their sequence under the
  proposed routing.
- **Consistency check (§4 step 7).** 35 C-anchored terminal slices already
  stitch at both ends today. For all 35, the N-only product is identical.
- **No residues after a stop.** No stitched slice today has canonical residues
  appended after its stop codon. This was checked because the both-ends rule
  could in principle do that, if a slice's last 10 aa occurred inside the
  canonical.

### 5.5 Rescued junctions

Each cell gives the slice's stop class, its translated length, the proposed
file and the validation. The validation is the matching transcript, "new" or
"unstitched"; "(NMD)" marks a match to an NMD transcript.

| Gene | Event | Strand | SJC | Today | Slice 1 | Slice 2 |
|---|---|---|---|---|---|---|
| BTBD1 | MXE 962 | − | 336 | T4 | A, 97 aa → alt_cterm (BTBD1-202) | terminal, 145 aa → T1 (new) |
| HAUS2 | MXE 3611 | + | 30 | T4 | A, 71 aa → alt_cterm (HAUS2-204, NMD) | open, 98 aa → T1 (HAUS2-208) |
| TARSL2 | SE 5590 | − | 2 | T4 | A, 51 aa → alt_cterm (TARSL2-202) | terminal, 87 aa → T1 (TARSL2-201) |
| PEX11A | SE 5716 | − | 19 | T4 | terminal, 247 aa → T1 (PEX11A-201) | A, 32 aa → alt_cterm (PEX11A-203) |
| SEC11A | SE 5773 | − | 29 | T4 | A, 60 aa → alt_cterm (SEC11A-209) | open, 59 aa → T1 (SEC11A-201) |
| TMED3 | SE 5845 | + | 2 | FAIL | A, 217 aa → T1 (TMED3-201) | A, 134 aa → T1 (new) |
| PML | SE 5987 | + | 2 | T4 | open, 127 aa → T1 (PML-215) | A, 28 aa → alt_cterm (PML-204) |
| PML | SE 5988 | + | 2 | FAIL | A, 117 aa → alt_cterm (new) | A, 28 aa → alt_cterm (PML-204) |
| SMAD6 | SE 6084 | + | 4 | T4 | A, 46 aa → alt_cterm (SMAD6-202, NMD) | open, 131 aa → T1 (SMAD6-201) |
| SMAD6 | SE 6085 | + | 5 | T4 | A, 45 aa → alt_cterm (new) | open, 131 aa → T1 (SMAD6-201) |
| TPM1 | SE 6201 | + | 136 | FAIL | A, 50 aa → alt_cterm (new) | A, 53 aa → alt_cterm (new) |
| TPM1 | SE 6203 | + | 136 | FAIL | A, 50 aa → alt_cterm (new) | A, 53 aa → alt_cterm (new) |
| FAM227B | SE 6333 | − | 12 | FAIL | A, 43 aa → alt_cterm (FAM227B-210) | A, 55 aa → alt_cterm (new) |
| FAM227B | SE 6334 | − | 2 | FAIL | A, 65 aa → alt_cterm (FAM227B-210) | A, 53 aa → alt_cterm (FAM227B-204) |
| MPI | SE 6375 | + | 4 | T4 | terminal, 199 aa → T1 (MPI-202) | A, 63 aa → alt_cterm (MPI-213) |
| EMC4 | SE 6592 | + | 17 | T4 | A, 116 aa → T1 (EMC4-202) | terminal, 82 aa → alt_cterm (EMC4-201) |
| EMC4 | SE 6593 | + | 17 | T4 | A, 97 aa → alt_cterm (new) | terminal, 81 aa → alt_cterm (EMC4-201) |
| UBE3A | SE 6649 | − | 3 | FAIL | A, 132 aa → alt_cterm (new) | A, 80 aa → alt_cterm (new) |
| ULK3 | SE 8573 | − | 3 | T4 | A, 41 aa → alt_cterm (ULK3-218) | open, 81 aa → T1 (ULK3-201) |
| IDH3A | SE 13424 | + | 29 | T4 | A, 28 aa → alt_cterm (IDH3A-215) | open, 49 aa → T1 (IDH3A-201) |
| PARP6 | SE 13583 | − | 11 | T4 | open, 74 aa → T1 orphan | A, 45 aa → alt_cterm orphan |
| PARP6 | SE 13584 | − | 6 | T4 | open, 77 aa → T1 orphan | A, 20 aa → alt_cterm orphan |
| PARP6 | SE 13585 | − | 6 | T4 | open, 78 aa → T1 orphan | A, 20 aa → alt_cterm orphan |
| MEIS2 | SE 18709 | − | 32 | T4 | A, 48 aa → alt_cterm (MEIS2-206) | terminal, 124 aa → T1 (MEIS2-222) |
| SERF2 | SE 26771 | + | 2 | FAIL | A, 80 aa → alt_cterm (SERF2-209) | A, 53 aa → alt_cterm (new) |
| FAM96A | SE 28030 | − | 45 | FAIL | A, 88 aa → alt_cterm (FAM96A-206) | A, 60 aa → alt_cterm (FAM96A-202) |
| PSMA4 | SE 37652 | + | 8 | T4 | open, 94 aa → T1 (PSMA4-201) | A, 48 aa → alt_cterm (PSMA4-205) |
| BMF | SE 37828 | − | 8 | T4 | terminal, 184 aa → T1 (BMF-201) | A, 129 aa → alt_cterm (BMF-206) |
| C15orf57 | SE 48189 | − | 2 | T4 | A, 103 aa → T1 (C15orf57-201) | terminal, 52 aa → alt_cterm (C15orf57-205) |
| RPS17 | RI 931 | − | 33 | T4 | A, 71 aa → T1 (new) | terminal, 48 aa → alt_cterm (RPS17-201) |
| PML | RI 978 | + | 8 | FAIL | A, 22 aa → alt_cterm (new) | A, 276 aa → alt_cterm (PML-202) |
| UBE3A | RI 1036 | − | 14 | T4 | A, 96 aa → alt_cterm (new) | terminal, 87 aa → T1 (UBE3A-207) |
| RPS27L | RI 1400 | − | 13 | T4 | A, 95 aa → alt_cterm (new) | terminal, 45 aa → T1 (RPS27L-205) |
| PML | A5SS 769 | + | 5 | T4 | A, 94 aa → alt_cterm (PML-208) | open, 104 aa → T1 (PML-201) |
| SHF | A5SS 810 | − | 7 | T4 | open, 96 aa → T1 (SHF-201) | A, 84 aa → alt_cterm (new) |
| HDDC3 | A3SS 1238 | − | 2 | FAIL | A, 84 aa → alt_cterm (HDDC3-201) | A, 123 aa → T1 (HDDC3-202) |
| FAM219B | A3SS 1304 | − | 30 | FAIL | A, 71 aa → T1 (FAM219B-201) | A, 70 aa → T1 (new) |
| SNAP23 | A3SS 1392 | + | 13 | T4 | A, 25 aa → alt_cterm (SNAP23-211) | open, 30 aa → T1 (SNAP23-201) |
| SNAP23 | A3SS 1393 | + | 13 | T4 | A, 34 aa → alt_cterm (SNAP23-207) | open, 30 aa → T1 (SNAP23-201) |
| ADPGK | A3SS 1551 | − | 6 | FAIL | A, 217 aa → T1 (new) | A, 216 aa → T1 (ADPGK-201) |
| IDH3A | A3SS 2695 | + | 10 | T4 | A, 28 aa → alt_cterm (IDH3A-215) | open, 49 aa → T1 (IDH3A-201) |
| TTC23 | A3SS 4584 | − | 15 | T4 | A, 62 aa → alt_cterm (TTC23-201) | open, 77 aa → T1 (TTC23-203) |
| SERF2 | A3SS 4710 | + | 22 | T4 | terminal, 130 aa → T1 (SERF2-202) | A, 56 aa → alt_cterm (SERF2-201) |
| DNAJC17 | A3SS 5490 | − | 4 | T4 | A, 32 aa → alt_cterm (DNAJC17-205) | open, 60 aa → T1 (DNAJC17-201) |
| NMB | A3SS 6113 | − | 3 | FAIL | A, 68 aa → alt_cterm (NMB-201) | A, 101 aa → T1 (NMB-202) |
| HAUS2 | A3SS 7146 | + | 2 | T4 | A, 38 aa → alt_cterm orphan | open, 71 aa → T1 orphan |

Patterns worth noting:

- **FAILURE junctions with two real isoforms.** FAM96A, FAM227B 6334, HDDC3
  and NMB are FAILURE today, yet both slices reconstruct annotated proteins.
  Each event switches between two annotated C-termini.
- **Terminal slices that are orphans today.** EMC4, C15orf57 and RPS17 have a
  slice whose stop is the last codon at the gene-wide `tx1`, yet that stop is
  not the canonical's. These slices are T1 orphans today because their last
  10 aa are not in the canonical. This is why the rule covers terminal
  annotated stops as well as internal ones.

---

## 6. Impact assessment

**Past results.**

- With the shipped `tx1` bug, these slices ran into the 3′ UTR and hit the
  real stop codon mid-slice. So they were also in T4 or FAIL, never in T1.
- Junctions with one PTC slice did write that slice to `psq_T4.fasta`,
  joined N-terminally, if it passed the 33 % length gate. Past `psq_T4.fasta`
  files therefore contain a number of correct alternative C-terminal
  isoforms.
- A search database built from T1 + canonical alone lacked them.

**After PR 4 + 4c (simulated on PR 4 alone; §7.5).**

- A T1 + canonical database gains the `psq_T1` records listed in §5.3.
- Users who opt in also gain `psq_T1_alt_cterm.fasta` (49 slices on chr15).
- On chr15 this file is 65 % exact annotated proteins. The rest are
  combinations of a known C-terminus with a splice event the annotation does
  not contain.

**Generalisation.**

- **chr15 only.** These figures come from one chromosome of one dataset.
- **Annotation dependence.** The rule's yield scales with annotation
  completeness. With sparse annotation (non-model organisms) fewer stops are
  recognised, and the remainder stay in T4 exactly as today. The rule never
  moves a junction down.

---

## 7. Implementation specification (PR 4c)

### 7.1 Rule

A slice is **C-anchored** if its first in-frame stop codon (GTF phase) starts
at the 5′ base of an annotated `stop_codon` of a `protein_coding` transcript of
the same gene. That base is `start` on the `+` strand and `end` on the `−`
strand. The stop may be internal or the last codon.

A C-anchored slice is translated up to that stop and is **not** a PTC.

### 7.2 Code changes

1. **`Junction._get_translated_region`** (after PR 4): also store
   `self.stop_sites`, the set of 5′ stop-codon bases of the gene's
   protein-coding transcripts. The `stop_codon` rows are already queried
   there.
2. **`Sequence.make_slice_localgenome`**: record each slice's exon segments,
   `(es, ee)` in genomic order, exactly as concatenated.
3. **`Sequence.translate(use_phase=True)`**: find the first in-frame stop and
   map it to the genome through those segments. If it is in
   `j.stop_sites`, set `sliceN_aa = make_pep(..., terminate=False)` and
   `sliceN_c_anchored = True`. The three-frame path (`use_phase=False`) is
   unchanged.
4. **`Sequence.stitch_to_canonical_aa`**: add a mode for C-anchored slices
   that joins the N-terminus only, without trying both ends first, so that
   canonical residues can never follow a stop. Classify the product as:
   - equal to the canonical → write it as the canonical;
   - last 10 aa equal to the canonical's last 10 aa → `T1`;
   - otherwise → `T1_alt_cterm`;
   - no N-terminal anchor → `T1_alt_cterm_orphan`.
5. **`_translate_one` gates.**
   - *Today's T1 gate* (both slices non-empty, no frameshift): unchanged,
     except that C-anchored slices use the mode from step 4. On chr15 this
     changes no output (§5.4).
   - *New gate, before T2:* both slices are non-empty after step 3, and at
     least one has an internal C-anchored stop. Route through step 4 without
     applying the frameshift test. Slices ending at different stops differ in
     length by arbitrary amounts, so that test means nothing for them.
   - *T2, T3, T4, FAIL:* unchanged.
6. **`write_slices`**: take the suffix per slice. The suffix already sets both
   the file name and the header's last field, so headers become
   `…|r{count}|T1_alt_cterm`. The number of `|` fields is unchanged, so
   existing parsers keep working.
7. **`fates.py`**: add
   `TIER1_ALT_CTERM = 'SUCCESS TIER 1 ALT CTERM. Retrieved phase: {0}. Used phase: {1}. Alternative C-terminus. \n\n'`.
   Use it when at least one slice is written to an `alt_cterm` file;
   otherwise use `TIER1`.

### 7.3 Regression tests (land before the code change; they must fail first)

| Junction | Today (PR 4) | Expected |
|---|---|---|
| FAM96A SE 28030 (−) | FAILURE | ALT CTERM; slice 1 = FAM96A-206 protein and slice 2 = FAM96A-202 protein, both in `psq_T1_alt_cterm` |
| NMB A3SS 6113 (−) | FAILURE | ALT CTERM; slice 1 = NMB-201 → `psq_T1_alt_cterm`; slice 2 = NMB-202 → `psq_T1` |
| PEX11A SE 5716 (−) | T4 | ALT CTERM; slice 2 = PEX11A-203 → `psq_T1_alt_cterm`; slice 1 = PEX11A-201 → `psq_T1` |
| EMC4 SE 6592 (+) | T4 | ALT CTERM; terminal slice 2 = EMC4-201 → `psq_T1_alt_cterm` (covers C-anchored terminal stops) |
| SEMA4B SE 5703 (+), slice 1 stop annotated only on NMD transcripts | T4 | T4, unchanged (negative control) |

The expected proteins should be taken from the GTF at test time rather than
hard-coded: translate the named transcript's CDS.

### 7.4 Golden test update

PR 4's golden values change as follows:

| Fate | After PR 4 | After PR 4c (simulated on PR 4 alone) |
|---|---|---|
| FAILURE | 330 | 316 |
| PARTIAL TIER 4 | 261 | 229 |
| SUCCESS TIER 1 | 189 | 192 |
| SUCCESS TIER 1 ALT CTERM | – | 43 |
| SUCCESS TIER 2 | 26 | 26 |
| SUCCESS TIER 3 | 2 | 2 |
| SKIPPED | 297 | 297 |

The implementation must reproduce these counts or explain every difference.
New FASTA digests, including `psq_T1_alt_cterm.fasta`, are recorded when the
PR is made. The test should also check that every `psq_T1.fasta` sequence
from PR 4 is still present.

### 7.5 Proposed CHANGELOG entry

> **Added.** Alternative C-terminal isoforms.
>
> - A stop codon at an annotated protein-coding stop of the gene is now
>   treated as a natural C-terminus rather than a premature stop codon. Such
>   slices are joined to the canonical at the N-terminus only.
> - Isoforms whose C-terminus differs from the canonical's are written to
>   `psq_T1_alt_cterm.fasta` (unstitched: `psq_T1_alt_cterm_orphan.fasta`).
>   The new fate is `SUCCESS TIER 1 ALT CTERM`.
> - On the chr15 test data, 46 junctions move up from T4 or FAIL. Existing T1
>   output is unchanged.
> - See `reports/2026-11-10-alt-cterm-tiering-audit.md`.

---

## 8. Limitations and open points

1. **Circular validation.** The rule and the truth set both come from the
   same GTF. The 67 % shows that the rule reconstructs known proteins
   faithfully. It is not independent evidence that the 20 new rescued
   sequences exist. Testing that would take proteomics evidence, for example
   searching `psq_T1_alt_cterm.fasta` against a dataset.
2. **UniProt mode is untested.** The audit uses the GTF canonical. With
   `--canonical` from UniProt, the N-terminal 10-aa anchor must match the
   UniProt canonical. Where Ensembl and UniProt disagree at the N-terminus,
   more slices become `alt_cterm_orphan`.
3. **NMD.** 17 class-A stops lie upstream of a junction within the slice. By
   the 50-nt rule some of these transcripts could be NMD targets, but Ensembl
   annotates them as protein-coding. The rule defers to the annotation.
   Stops annotated only on NMD transcripts (class B) stay PTCs.
4. **TSL filter.** Stop sites are taken from all protein-coding transcripts,
   whatever their TSL. The PR 4 bounds use `transcript_support_level ≤
   params.tsl_threshold`. Applying the same filter would be more
   conservative. Decide when implementing and re-run this audit; it is a
   one-line change in `annotated_stop_sites`.
5. **The mirror case: alternative N-termini.** Slices starting at an
   alternative annotated start codon would need joining at the C-terminus
   only. Not measured here.
6. **Junctions without a GTF phase** (268 on chr15) are not covered. They
   depend on the phase-lookup issue
   ([tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions)
   item 4).
7. **The 33 % length gate in T4** no longer applies to rescued junctions. That
   is intended, because a naturally short alternative C-terminus is not a
   fragment.

---

## 9. Reproducing this audit

From the repository root, with the package dependencies installed:

```bash
python reports/cterm_audit/audit_cterm.py
```

The script:

- copies the GTF to a temporary directory, so no `.cached` file is left in the
  repository;
- stubs out UniProt;
- patches the PR 4 bounds in from
  [tx1_audit/audit_tx1.py](tx1_audit/audit_tx1.py);
- rewrites the three evidence files in `reports/cterm_audit/`.

The inputs and their checksums are the same as for the
[tx1 audit, §9](2026-11-03-tx1-translation-bounds-audit.md#9-reproducing-this-audit).

**`slices.tsv` columns.**

| Column | Contents |
|---|---|
| `state` | stop class from §4 |
| `stop_offset_nt` | stop position in the slice, transcript orientation |
| `stop_genomic` | 5′ base of the stop, 1-based |
| `stop_biotypes` | biotypes of the transcripts annotating that stop |
| `stop_in_3prime_exon` | whether the stop is in the slice's 3′-most exon |
| `current_product_class` / `current_product_matches` | today's stitch (both ends if PTC-free, else N-only as in T4) |
| `current_file` | today's file, for T1 junctions |
| `n_stitch_class` / `n_stitch_matches` | N-only stitch of any slice with a stop |
| `c_anchored`, `ends_at_canonical_cterm`, `n_only_equals_current` | flags used by the proposed routing and the §5.4 check |
| `proposed_fate` / `proposed_file` | proposed routing |
| `aa_to_stop` | translated slice up to its first stop |

**`junctions.tsv`.** One row per junction: actual fate, whether it was
analysed, both slice states, the frameshift flag and the proposed fate.
