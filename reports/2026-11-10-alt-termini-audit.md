# Audit: alternative C- and N-termini on the current pipeline (PR 4c)

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Pre-implementation evaluation. No production code changed. Re-runs the [C-terminal audit](2026-11-10-alt-cterm-tiering-audit.md) on the current pipeline and extends it to N-termini. **All recommendations in §7 were approved on 2026-11-10:** PR 4b2 first; TSL ≤ 2 sites; separate files per terminus; one fate `SUCCESS TIER 1 ALT TERMINUS`; the refined frameshift rule. **It also found a regression from PR 4b in the frameshift test (§6), proposed as PR 4b2 to land before 4c.** |
| **Code audited** | Commit `919f6c7` (PR 4b final + docs): tiering in [jcast/main.py](../jcast/main.py), frames in `Junction.get_translated_phase`, stitching in `Sequence.stitch_to_canonical_aa`, frameshift test `Sequence.frameshift` |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions (716 not skipped), Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/termini_audit/audit_termini.py` at `919f6c7` (≈ 2 min, no network) |
| **Evidence** | [termini_audit/summary.txt](termini_audit/summary.txt) · [termini_audit/slices.tsv](termini_audit/slices.tsv) (one row per slice and TSL policy) · [termini_audit/junctions.tsv](termini_audit/junctions.tsv) (one row per junction and policy) |

Every number in this report is printed in `summary.txt` or can be read
directly from the two TSV files. The §6 cross-check is an exception: it
was computed from `slices.tsv` (columns `frame`, `coding_nt`), and the
query is given there.

---

## 1. Summary

**What 4c does.** It recognises two kinds of legitimate isoform end that
today's tiers misclassify.

- **Alternative C-terminus (C-anchored).** The slice's first in-frame stop
  codon is an annotated stop codon of a protein-coding transcript.
  - It is a natural end, not a premature stop.
  - The slice is joined to the canonical at its N-terminus only.
  - Today such junctions are T4 or FAILURE.
- **Alternative N-terminus (N-anchored).** The slice lacks the gene's start
  codon, which the other slice contains, but it contains another annotated
  protein-coding start codon.
  - It is translated from that start, which also sets the frame.
  - It is joined at its C-terminus only.
  - Today such slices are read from mid-CDS and usually recreate the
    canonical (phase audit §5.4).
- **Both.** A slice anchored at both ends is a self-contained protein, start
  to stop, with no joining. None occurs on chr15.

Each product goes to a file named for the terminus that differs from the
canonical:

- `psq_T1_alt_cterm.fasta` if the C-terminus differs;
- `psq_T1_alt_nterm.fasta` if the N-terminus differs;
- `psq_T1.fasta` if neither does (for example, an internal event joined
  back);
- `psq_canonical.fasta` if the product equals the canonical.

**Results on chr15** (recommended policy: annotated sites at TSL ≤ 2, the
same threshold as the coding bounds):

| | Value |
|---|---|
| Junctions routed | 36 |
| Junctions that change fate | 27: FAIL → T1-alt 3, FAIL → T1 1, T4 → T1-alt 13, T2 → T1-alt 5, T2 → T1 1, T1 → T1-alt 4 |
| Routed slices that exactly reproduce an annotated protein-coding protein | 62 of 72 (86 %) |
| Today's T1 slices that exactly reproduce one (baseline) | 265 of 424 (62.5 %) |
| `psq_T1_alt_cterm.fasta` | 19 slices, 16 annotated |
| `psq_T1_alt_nterm.fasta` | 7 slices, all annotated |
| Existing `psq_T1.fasta` records moved or lost | 0 |

**The annotated-site rule still separates well.** Joining each slice
N-terminally only:

| First stop codon of the slice | Exact annotated protein-coding |
|---|---|
| Annotated protein-coding stop (routed) | 19 of 26 (73 %) |
| Annotated only on NMD transcripts (control) | 11 of 93 (mostly NMD proteins: 55) |
| Unannotated (control) | 13 of 325; under any-TSL, 0 of 238 (§4.1) |

**The N-terminal rule.** All 8 N-anchored slices reproduce an annotated
protein-coding protein, for example ARPP19-202, DET1-201, DPP8-201 and
ZNF280D-207. Today, 7 of these 8 slices recreate the canonical.

**A wider TSL policy** (annotated sites from any TSL) roughly doubles the
yield, to 74 junctions routed. It lowers precision: class-A stops give 57 %
exact matches, and the additional class-A slices it admits match less often
than today's T1. See decision 1.

---

## 2. Changes since the first C-terminal audit

The [first audit](2026-11-10-alt-cterm-tiering-audit.md) ran on PR 4 alone.
Since then:

- PR 4a trims every exon to the coding region and skips events outside it.
- PR 4a2 fixes RI slice composition.
- PR 4b gives each slice its own frame and removes the three-frame fallback.

This audit uses the shipped code for all of that:

- `Junction.slice_exons` and `slice_phases`;
- the actual fates;
- the actual FASTA output.

It adds three things:

- the N-terminal rule;
- a comparison of two TSL policies for the annotated sites, which was an
  open point in the first audit (its §8 item 4);
- a refined frameshift rule (§3, step 5).

---

## 3. Method

[termini_audit/audit_termini.py](termini_audit/audit_termini.py) works as
follows.

1. **Actual baseline.** Every junction is run through `main._translate_one`
   with the settings of `tests/test_pipeline.py`. The script records the
   fate and every record written, attributed to junction and slice through
   the FASTA header.
2. **Per slice, for each of two TSL policies** — `tsl2`
   (`transcript_support_level <= params.tsl_threshold`) and `any`:
   - **N-side.** Does the slice contain the gene's start codon, the most
     upstream one as used for `tx0`/`tx1`? If not, and the other slice does,
     the slice is *start-less*. It is N-anchored if it contains another
     annotated protein-coding start codon that reads ATG in the genome; the
     first one in transcript order is used. 27 annotated starts that are not
     ATG at that position, such as split codons, are ignored.
   - **C-side.** The first stop codon is found from the translation start:
     the slice's frame, or the alternative ATG. It is classed as open (no
     stop), terminal (last codon), A (annotated protein-coding), B
     (annotated only on non-coding transcripts) or D (unannotated). The
     slice is C-anchored if A, or if terminal at an annotated
     protein-coding stop.
3. **Product.**

   | Anchoring | Join to the canonical |
   |---|---|
   | C-anchored only | N-terminus only |
   | N-anchored only | C-terminus only |
   | both | none (self-contained) |
   | neither | both ends, as today |

4. **File.** The file is chosen by which terminus differs from the
   canonical's first or last 10 aa; §1 lists the files.
5. **Routing.** A junction is routed if all of the following hold:
   - at least one slice is anchored;
   - both slices reach an acceptable end (open, terminal or A);
   - it is not a frameshift.

   The frameshift test is waived only when a slice ends at an internal
   annotated stop. The slices then end at different stop codons, so their
   lengths differ arbitrarily. Otherwise, the test compares the slices'
   lengths from their translation starts to the slice end. A routed
   junction's fate is `SUCCESS TIER 1 ALT TERMINUS` if any slice goes to an
   alt file, otherwise `SUCCESS TIER 1`. Unrouted junctions keep their fate.
6. **Validation.** Every product, today's and proposed, is compared with
   all proteins translated from the GTF (`audit_cterm.annotated_proteins`):
   *annotated protein_coding*, *annotated other* (NMD), *new*, or
   *unjoined*.

---

## 4. Results

### 4.1 Does the annotated-site rule still discriminate?

Each slice joined N-terminally only, by its first stop codon:

| Policy | Stop class | Slices | Exact annotated protein-coding | Annotated other (NMD) | New | Unjoined |
|---|---|---|---|---|---|---|
| tsl2 | **A**: annotated protein-coding (routed) | 26 | **19 (73 %)** | 0 | 7 | 0 |
| tsl2 | B: annotated only on non-coding (control) | 93 | 11 | 55 | 18 | 9 |
| tsl2 | D: unannotated (control) | 325 | 13 | 32 | 228 | 52 |
| any | **A** (routed) | 67 | **38 (57 %)** | 3 | 18 | 8 |
| any | B (control) | 139 | 5 | 84 | 38 | 12 |
| any | D (control) | 238 | **0** | 0 | 197 | 41 |

Under `tsl2`, the 13 annotated matches in class D are stops annotated only
on TSL > 2 transcripts. Under `any` they move to class A, and D is left
with no protein-coding matches at all.

The 41 class-A slices that `any` adds yield 19 exact matches (46 %). That is
below today's T1 baseline of 62.5 % (§4.4).

### 4.2 N-terminal anchoring

| Policy | Start-less slices | With an annotated alternative ATG (N-anchored) | Of those, exact annotated protein-coding | Today's product was the canonical |
|---|---|---|---|---|
| tsl2 | 42 | 8 | **8** | 7 |
| any | 42 | 17 | 13 (+2 new, 2 unjoined) | 10 |

Under `tsl2` every N-anchored product is a real annotated isoform:

- ARPP19-202, from ARPP19 MXE 579 and SE 7770/7772;
- DET1-201 (SE 5744);
- DPP8-201 (SE 6124);
- ZNF280D-207 (SE 6275/6276).

7 of these 8 slices currently recreate the canonical. For PKM SE 36967 the
alternative start is the canonical's own start, so the product is the
canonical; that junction is a frameshift and stays T2 (§6).

Start-less slices without an annotated alternative start (34 under `tsl2`)
are unchanged: they are still read from mid-CDS, as decided in PR 4b.

### 4.3 Fates (policy `tsl2`)

| Today | Proposed | Junctions |
|---|---|---|
| FAILURE | SUCCESS TIER 1 ALT TERMINUS | 3 |
| FAILURE | SUCCESS TIER 1 | 1 |
| PARTIAL TIER 4 | SUCCESS TIER 1 ALT TERMINUS | 13 |
| SUCCESS TIER 2 | SUCCESS TIER 1 ALT TERMINUS | 5 |
| SUCCESS TIER 2 | SUCCESS TIER 1 | 1 |
| SUCCESS TIER 1 | SUCCESS TIER 1 ALT TERMINUS | 4 |
| *everything else (689 non-skipped)* | *unchanged* | |

Nine more junctions are routed but keep T1 with the same products: MFGE8,
CASC4, COMMD4, CHP1, SNAPC5, BNIP2, BTBD1, FAM219B and ADPGK. In each, both
slices end at annotated stops that today's both-ends join already handles.

Under `any`, 65 junctions change fate: FAIL → T1-alt 11, FAIL → T1 1,
T4 → T1-alt 41, T2 → T1-alt 6, T2 → T1 1, T1 → T1-alt 5.

**Held back as frameshifts.** Under the refined test, five anchored
junctions are true frameshifts and keep their current fate: APH1B MXE 3368,
SEC11A MXE 434, SNURF MXE 499 and PKM SE 36967 (all T2), and GTF2A2 SE 11950
(T1, see §6). Under `any`, HAUS2 SE 42228 is also held back.

### 4.4 Where routed slices go, and validation (policy `tsl2`)

| Proposed file | Slices | Exact annotated protein-coding | New | Unjoined |
|---|---|---|---|---|
| `psq_canonical` | 29 | 29 | – | – |
| `psq_T1` | 15 | 10 | 5 | – |
| `psq_T1_alt_cterm` | 19 | **16** | 3 | – |
| `psq_T1_alt_nterm` | 7 | **7** | – | – |
| alt orphan | 2 | – | – | 2 |
| **Total** | 72 | **62 (86 %)** | 8 | 2 |

**Baseline.** Today's T1 slices, under today's rule, include the canonical
matches: 265 of 424 (62.5 %) are exact annotated protein-coding proteins.

**Moves from today's files:**

| Today | Proposed | Slices |
|---|---|---|
| `psq_T4` / `psq_T4_orphan` | `psq_T1_alt_cterm` | 13 |
| `psq_T4` | `psq_T1` | 2 |
| `psq_T1_orphan` / `psq_T2_orphan` | alt orphan | 2 |
| not written | various | 45 (29 of them the canonical, already present) |

No `psq_T1.fasta` record moves.

Under `any`, `psq_T1_alt_cterm` would hold 52 slices (35 annotated
protein-coding) and `psq_T1_alt_nterm` 14 (12).

### 4.5 Per-junction list (policy `tsl2`, 36 routed junctions)

Each cell gives the anchoring (C, N, N+C or –), the stop class, the slice
length, the proposed file and the matching annotated transcript.

| Gene | Event | Strand | SJC | Today → proposed | Slice 1 | Slice 2 |
|---|---|---|---|---|---|---|
| ARPP19 | MXE 579 | − | 5 | T1 → T1 ALT | N, open, 40 aa → alt_nterm (ARPP19-202) | –, open, 56 aa → T1 (ARPP19-201) |
| PEX11A | SE 5716 | − | 19 | T4 → T1 ALT | C, terminal, 247 aa → canonical (PEX11A-201) | C, A, 32 aa → alt_cterm (PEX11A-203) |
| MFGE8 | SE 5739 | − | 68 | T1 → T1 | C, terminal, 158 aa → canonical (MFGE8-201) | C, terminal, 106 aa → T1 (MFGE8-202) |
| DET1 | SE 5744 | − | 7 | T2 → T1 ALT | –, open, 372 aa → canonical (DET1-202) | N, open, 361 aa → alt_nterm (DET1-201) |
| TMED3 | SE 5845 | + | 2 | FAIL → T1 | C, A, 217 aa → canonical (TMED3-201) | C, A, 134 aa → T1 (new) |
| COMMD4 | SE 5921 | + | 2 | T1 → T1 | C, terminal, 98 aa → canonical (COMMD4-201) | C, terminal, 39 aa → T1 (COMMD4-202) |
| PML | SE 5987 | + | 2 | T4 → T1 ALT | –, open, 127 aa → T1 (PML-215) | C, A, 28 aa → alt_cterm (PML-204) |
| PML | SE 5988 | + | 2 | FAIL → T1 ALT | C, A, 117 aa → alt_cterm (new) | C, A, 28 aa → alt_cterm (PML-204) |
| DPP8 | SE 6124 | − | 2 | T2 → T1 ALT | –, open, 102 aa → canonical (DPP8-203) | N, open, 86 aa → alt_nterm (DPP8-201) |
| BNIP2 | SE 6210 | − | 15 | T1 → T1 | C, terminal, 47 aa → T1 (new) | C, terminal, 35 aa → canonical (BNIP2-201) |
| ZNF280D | SE 6275 | − | 5 | T2 → T1 ALT | –, open, 58 aa → canonical (ZNF280D-201) | N, open, 45 aa → alt_nterm (ZNF280D-207) |
| ZNF280D | SE 6276 | − | 5 | T2 → T1 ALT | –, open, 58 aa → canonical (ZNF280D-201) | N, open, 45 aa → alt_nterm (ZNF280D-207) |
| MPI | SE 6375 | + | 4 | T4 → T1 ALT | C, terminal, 199 aa → canonical (MPI-202) | C, A, 63 aa → alt_cterm (MPI-213) |
| CASC4 | SE 6426 | + | 5 | T1 → T1 | C, terminal, 135 aa → canonical (CASC4-201) | C, terminal, 79 aa → T1 (CASC4-202) |
| CHP1 | SE 6489 | + | 2 | T1 → T1 | C, terminal, 78 aa → canonical (CHP1-201) | C, terminal, 37 aa → T1 (CHP1-204) |
| EMC4 | SE 6592 | + | 17 | T4 → T1 ALT | C, A, 116 aa → canonical (EMC4-202) | C, terminal, 82 aa → alt_cterm (EMC4-201) |
| ARPP19 | SE 7770 | − | 18 | T1 → T1 ALT | –, open, 56 aa → T1 (ARPP19-201) | N, open, 40 aa → alt_nterm (ARPP19-202) |
| ARPP19 | SE 7772 | − | 18 | T1 → T1 ALT | –, open, 56 aa → T1 (ARPP19-201) | N, open, 40 aa → alt_nterm (ARPP19-202) |
| GTF2A2 | SE 11945 | − | 2 | T2 → T1 ALT | C, terminal, 50 aa → canonical (GTF2A2-202) | C, terminal, 7 aa → alt orphan (too short to join) |
| BTBD1 | SE 12765 | − | 4 | T1 → T1 | C, terminal, 130 aa → canonical (BTBD1-201) | C, terminal, 81 aa → T1 (new) |
| MEIS2 | SE 18709 | − | 32 | T4 → T1 ALT | C, A, 48 aa → alt_cterm (MEIS2-206) | C, terminal, 124 aa → canonical (MEIS2-222) |
| OAZ2 | SE 24775 | − | 2 | T2 → T1 | C, terminal, 60 aa → canonical (OAZ2-205) | C, terminal, 42 aa → canonical (OAZ2-205) |
| FAM96A | SE 28030 | − | 45 | FAIL → T1 ALT | C, A, 88 aa → alt_cterm (FAM96A-206) | C, A, 60 aa → alt_cterm (FAM96A-202) |
| SNAPC5 | SE 29748 | − | 3 | T1 → T1 | C, terminal, 98 aa → canonical (SNAPC5-202) | C, terminal, 68 aa → T1 (SNAPC5-201) |
| LYSMD4 | SE 43505 | − | 2 | T1 → T1 ALT | C, terminal, 296 aa → alt orphan (unjoined) | C, terminal, 202 aa → canonical (LYSMD4-202) |
| C15orf57 | SE 48189 | − | 2 | T4 → T1 ALT | C, A, 103 aa → canonical (C15orf57-201) | C, terminal, 52 aa → alt_cterm (C15orf57-205) |
| RPS17 | RI 931 | − | 33 | T4 → T1 ALT | C, A, 49 aa → canonical (RPS17-209) | C, terminal, 48 aa → alt_cterm (RPS17-201) |
| PML | RI 978 | + | 8 | FAIL → T1 ALT | C, A, 58 aa → alt_cterm (PML-207) | C, A, 276 aa → alt_cterm (PML-202) |
| RPS27L | RI 1400 | − | 13 | T4 → T1 ALT | C, A, 58 aa → alt_cterm (new) | C, terminal, 45 aa → canonical (RPS27L-205) |
| FRMD5 | RI 4475 | − | 6 | T4 → T1 ALT | C, terminal, 191 aa → canonical (FRMD5-202) | C, A, 129 aa → alt_cterm (new) |
| PML | A5SS 769 | + | 5 | T4 → T1 ALT | C, A, 94 aa → alt_cterm (PML-208) | –, open, 104 aa → canonical (PML-201) |
| HDDC3 | A3SS 1238 | − | 2 | T4 → T1 ALT | C, A, 84 aa → alt_cterm (HDDC3-201) | C, terminal, 123 aa → canonical (HDDC3-202) |
| FAM219B | A3SS 1304 | − | 30 | T1 → T1 | C, terminal, 71 aa → canonical (FAM219B-201) | C, terminal, 70 aa → T1 (new) |
| ADPGK | A3SS 1551 | − | 6 | T1 → T1 | C, terminal, 217 aa → T1 (new) | C, terminal, 216 aa → canonical (ADPGK-201) |
| SERF2 | A3SS 4710 | + | 22 | T4 → T1 ALT | C, terminal, 130 aa → T1 (SERF2-202) | C, A, 56 aa → alt_cterm (SERF2-201) |
| NMB | A3SS 6113 | − | 3 | T4 → T1 ALT | C, A, 68 aa → alt_cterm (NMB-201) | C, terminal, 101 aa → canonical (NMB-202) |

**Notes.**

- **OAZ2** is the antizyme gene, which uses programmed ribosomal
  frameshifting. Both slices give the canonical, and it moves from T2 to T1
  because the refined frameshift test (§6) finds the slices in frame at
  their shared stop. Its output is unchanged.
- **Not routed under `tsl2`.** SNAP23, IDH3A, TTC23, DNAJC17, PSMA4, BMF and
  others from the first audit are not routed. Their annotated alternative
  stops belong only to TSL > 2 transcripts. Under `any` they are routed.

---

## 5. Comparison with the first C-terminal audit

The first audit ran on PR 4 with no TSL filter on sites and simulated
routing. It found 46 junctions moving to T1 or T1-alt: 14 from FAIL and 32
from T4.

On the current pipeline, the comparable `any` policy routes 53 junctions
from FAIL or T4 (12 and 41). That fits PRs 4a–4b translating more junctions
correctly.

Under `tsl2`, 17 of these remain (4 from FAIL, 13 from T4).

---

## 6. Regression found: the frameshift test ignores per-slice frames (proposed PR 4b2)

**The problem.** `Sequence.frameshift` decides T1 vs T2 as
`(len(slice1_nt) - len(slice2_nt)) % 3 != 0`. That is valid only when both
slices are read from the same starting frame. Since PR 4b, each slice has
its own frame (`slice_phases`). Whether the event shifts the frame of the
shared downstream sequence depends on the *coding* lengths:
`(len1 − phase1) − (len2 − phase2)`.

**The two tests disagree for 17 junctions** — every T1/T2 junction whose
slices have different frames.

> **Erratum (2026-11-10).** This table first listed PKM SE 36967 as
> T1 → T2, giving 15 T2 → T1 and 2 T1 → T2. The query took `coding_nt`
> from `slices.tsv`. For an N-anchored slice, that column is measured from
> the alternative start codon, not from the slice's frame, so it was wrong
> for PKM's slice 2. Recomputed from the pipeline's actual slices in the
> [PR 4b2 audit](2026-11-10-frameshift-audit.md), PKM SE 36967 is T2 → T1.
> The total of 17 is unchanged. The corrected table follows.

| Today | Junctions | Coding-length test |
|---|---|---|
| T2 | ADAL MXE 2115, LINS1 MXE 449, GTF2A2 SE 11945, HYKK SE 17536, WDR61 SE 22528, OAZ2 SE 24775, MGA SE 31207, DET1 SE 5744, CD276 SE 5993, DPP8 SE 6124, ZNF280D SE 6275/6276, TMEM62 SE 6443, UBE3A SE 6665/6666, PKM SE 36967 (16) | in frame, so T1 |
| T1 | GTF2A2 SE 11950 (1) | frameshift, so T2 |

**Evidence for the coding-length test.** LINS1 MXE 449's T2 record is
exactly the annotated LINS1-207 protein (phase audit §5.4), which is
impossible for a genuine frameshift.

**Why the PR 4b audit missed it.** That audit validated frames and products,
not the T1/T2 split.

**Fix.** Compute `frameshift` from the coding lengths. It is independent of
4c and changes the fate of these 17 junctions only. Records may move between
`psq_T1*` and `psq_T2*`; sequences are unchanged. **Proposed as PR 4b2, with
its own short audit, before 4c.** This audit's routing already uses the
coding-length test, so the 4c figures above are consistent with it.

---

## 7. Decisions needed

**0. PR 4b2: fix the frameshift test first (recommended).** A regression
from PR 4b, §6.

**1. TSL policy for annotated start and stop sites: `tsl2` (recommended) or
`any`.**

- **`tsl2`:**
  - consistent with the coding bounds and frame lookup, which use the same
    `params.tsl_threshold`;
  - 86 % of routed slices are exact annotated proteins, against 62.5 % for
    today's T1;
  - 36 junctions routed.
- **`any`:**
  - 74 junctions routed;
  - the additional class-A stops match annotation at 46 %, below today's T1;
  - adds stops seen only in poorly supported transcripts.
- **Users who want more yield** can raise `params.tsl_threshold`. That also
  widens the coding bounds, which is consistent. It is not a CLI option yet.

**2. Output files: separate files per terminus (recommended).**

- `psq_T1_alt_cterm.fasta` (your earlier choice) and
  `psq_T1_alt_nterm.fasta`.
- Their unjoinable slices go to `psq_T1_alt_cterm_orphan.fasta` and
  `psq_T1_alt_nterm_orphan.fasta`.
- A slice whose two termini both differ goes to `psq_T1_alt_both.fasta`.
  There is none on chr15; the file is created only when needed.

**3. Fate name: one fate, `SUCCESS TIER 1 ALT TERMINUS` (recommended).** The
file records which terminus differs. Alternatively, there could be one fate
per terminus.

**4. Frameshift rule for anchored junctions: the refined rule (§3, step 5)
(recommended).**

- The frameshift test is waived only when a slice ends at an internal
  annotated stop.
- The first audit's rule waived it whenever a slice was anchored. That
  would have promoted five true frameshifts (APH1B, SEC11A, SNURF, PKM,
  GTF2A2 SE 11950).

---

## 8. Implementation outline (after decisions; PR 4c)

1. **`Junction`:**
   - store the gene's annotated protein-coding start and stop codon sites at
     the TSL policy chosen (decision 1);
   - expose each slice's genomic positions, from `slice_exons`.
2. **`Sequence.translate`:**
   - per slice, determine N-anchoring (alternative ATG) and C-anchoring
     (annotated stop);
   - an anchored slice ends at its annotated stop and is not a PTC;
   - set `sliceN_anchors`.
3. **`stitch_to_canonical_aa`:** add N-only (existing), C-only and none
   modes, chosen from the anchors.
4. **`_translate_one`:** add an anchored gate before T1/T2, with the
   frameshift rule from decision 4. `write_slices` takes a per-slice suffix
   (the file from §3, step 4).
5. **`fates`:** add `TIER1_ALT_TERMINUS` (decision 3).
6. **Regression tests,** landed first:
   - FAM96A SE 28030 (FAIL → T1 ALT; both alt_cterm = FAM96A-206/-202);
   - NMB A3SS 6113 (T4 → T1 ALT; NMB-201 alt_cterm);
   - EMC4 SE 6592 (terminal annotated stop → alt_cterm EMC4-201);
   - DET1 SE 5744 (N-anchored → alt_nterm DET1-201);
   - ARPP19 SE 7770 (T1 → T1 ALT; ARPP19-202 alt_nterm);
   - SEMA4B SE 5703 (NMD-only stop: stays T4);
   - APH1B MXE 3368 (frameshift: stays T2).
7. **Verification.** Implement, then run the real pipeline through this
   audit's harness. It must reproduce §4.3 and §4.4 for the chosen policy,
   after PR 4b2.
8. **Documentation.** Tier tables in the README, REVIEW.md §3 and the
   CHANGELOG; the reports index.

---

## 9. Reproducing this audit

```bash
git checkout 919f6c7
python reports/termini_audit/audit_termini.py
```

The script:

- runs the shipped pipeline once, stubbing out UniProt and copying the GTF
  to a temporary directory;
- simulates the proposed routing for both TSL policies;
- rewrites the three evidence files in `reports/termini_audit/`.

The inputs are those of the
[tx1 audit, §9](2026-11-03-tx1-translation-bounds-audit.md#9-reproducing-this-audit).
