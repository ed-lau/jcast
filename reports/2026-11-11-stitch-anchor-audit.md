# Audit: joining slices to the canonical, and the anchor length (PR 5; REVIEW.md MED-3)

| | |
|---|---|
| **Date** | 2026-11-11 |
| **Status** | **Implemented in PR 5** (§9). Decisions in §8 approved 2026-11-11: rule "unique + complete ends", anchor length 7, a separate terminus-comparison length kept at 10, no command-line option. Baseline is commit `887a98a`. |
| **Code** | `Sequence.stitch_to_canonical_aa` and `Sequence._join` ([jcast/sequences.py](../jcast/sequences.py)); `params.aa_stitch_length` |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | chr15: `python reports/stitch_audit/audit_stitch.py` at `887a98a` (≈ 4 min, no network; the UniProt part needs the local `cache/uniprot-cache.db`). Proteome: `python reports/stitch_audit/survey_anchors.py --gtf gencode.v50.annotation.gtf --transcripts gencode.v50.transcripts.fa.gz --out reports/stitch_audit/gencode_v50` (≈ 2 min) |
| **Evidence** | [stitch_audit/summary.txt](stitch_audit/summary.txt) · [stitch_audit/calls.tsv](stitch_audit/calls.tsv) · [stitch_audit/sweep.tsv](stitch_audit/sweep.tsv) · [stitch_audit/gencode_v50/](stitch_audit/gencode_v50/) |

## 1. Summary

- **On chr15, 12 joins are wrong.**
  - 10 of them produce exactly the canonical protein. The alternative
    isoform is then written as the canonical record and lost from the
    output.
  - The other 2 (CD276) produce a protein about 220 aa shorter than the real
    one.
  - The wrong joins are 12 of 998 distinct join attempts. 880 are correct,
    and 42 could have been joined but weren't.
- **MED-3's own proposal, a uniqueness check on each anchor, fixes only 2
  of the 12.** The main cause is different: the code joins a peptide at an
  end that is already a real protein end.
  - **8 cases:** a Tier 4 peptide ends at its premature stop codon. The
    code first tries joining at both ends; the last 10 residues are about 9
    shared residues plus 1 that matches by chance, so the peptide is
    "joined" back into the full canonical.
  - **2 cases (COMMD4):** a peptide starts at the gene's own start codon.
    Its first residues also match an internal Met24, so the canonical's
    first 24 residues are prepended.
  - **The remaining 2 (CD276):** an anchor occurs twice in the duplicated
    Ig domains, and the first occurrence is the wrong one.
- **Proposed rule:** an end that is already complete is never joined (a
  stop codon, or the gene's start codon). A repeated anchor is resolved by
  the other anchor, or the slice is not joined.
  - At the current length of 10, this gives 917 correct, 0 wrong and 17
    missed.
- **A 10-residue anchor is longer than necessary.** Most missed joins are
  slices whose shared part at one end is shorter than 10 residues. With the
  proposed rule:

  | Anchor length | Correct | Wrong | Missed |
  |---|---|---|---|
  | 10 | 917 | 0 | 17 |
  | 5–7 | 933 | 0 | 1 |

  Across the 19,168 MANE Select proteins (GENCODE v50):
  - a 7-residue anchor taken from an exon boundary of another protein
    occurs in a given protein in 0.007% of trials, against 0.005% at 10;
  - 0.71% of 7-residue anchors repeat within their own protein, against
    0.51% at 10.

  7 is the shortest length at which both rates have levelled off.

## 2. How joining works

The canonical protein comes from UniProt by default, or from the GTF. A
translated slice is placed in it by searching for the slice's first and
last `aa_stitch_length` (10) residues with `str.find`, which returns the
first occurrence.

| Caller | Mode |
|---|---|
| T1 (`stitch_to_canonical_aa`) | both ends |
| T2 | both ends, else the N-terminus only |
| T4 translated slice | both ends |
| T4 force-translated slice (ends at its premature stop) | both ends, else the N-terminus only (`slice_has_ptc`) |
| T1 ALT TERMINUS (`_join`, PR 4c) | `n` for a C-anchored slice, `c` for an N-anchored one, `both`, or none |

Nothing checks that an anchor occurs only once, or that the start anchor
comes before the end anchor (MED-3).

## 3. Method

[audit_stitch.py](stitch_audit/audit_stitch.py) runs every chr15 junction
through `main._translate_one` with the canonical from the GTF. It records
every join attempt: peptide, mode, canonical, product, and the translation
start. 1,095 attempts were recorded, of which 998 are distinct (same
peptide, canonical and mode).

**Ground truth.**
- **Coordinates:** with the GTF canonical, every canonical residue and every
  peptide residue has genomic coordinates. The audit rebuilds each
  canonical's coding sequence from its coordinates and re-translates it.
  All 998 reproduce the canonical exactly.
- **Correct product:** the canonical sequence before the residue whose
  codon starts where the peptide's first codon starts, then the peptide,
  then the canonical after the residue whose codon starts where its last
  codon starts.
- **Codons are placed by their first base, in frame.** This tolerates a
  splice site a few bases away (e.g. NAGNAG), which changes a codon's last
  base but not where it joins. Three − strand slices needed this.
- **Complete ends are not joined.** A peptide that ends at a stop codon
  ends there. One that starts at the gene's start codon (where trimming
  starts a slice that contains it) begins there.
- **No correct product:** the slice should not be joined at that end when
  the codon is not in the canonical.

Each attempt is classified as:
- **correct join:** the product equals the correct product;
- **wrong join:** the product differs, or a product was made where none
  should be;
- **missed join:** no product, though a correct one exists;
- **unjoined (correct):** no product, and none should be made.

The audit checks that its own replay of the current rule reproduces every
pipeline product.

## 4. Results with the current code (anchor length 10)

| Classification | Join attempts |
|---|---|
| correct join | 880 |
| missed join | 42 |
| unjoined (correct) | 64 |
| **wrong join** | **12** |

The 12 wrong joins:

| Junction | Slice | Cause | Product |
|---|---|---|---|
| SE 19972 UBE2Q2, RI 952 MAN2C1, RI 1036 UBE3A, A5SS 751 TM6SF1, A5SS 2612 HMG20A, A5SS 3645 WDR73, A5SS 4932 DNAJA4, A3SS 7415 JMJD7-PLA2G4B | T4, ends at its premature stop | Joined at both ends: the end anchor matches the canonical just after the shared region | exactly the canonical (8) |
| MXE 444 and SE 5917 COMMD4 | starts at the gene's start codon (exon 2 skipped) | Its first residues (MSSVKL…) also occur at canonical Met24, so canonical residues 1–24 are prepended | exactly the canonical (2) |
| A3SS 1316 CD276, both slices | — | The start anchor occurs twice (duplicated IgV–IgC domains of 4Ig-B7-H3); the first occurrence is taken | 316 and 308 aa instead of the correct join |

The 42 missed joins:
- **25 have a complete end.** The code needs both anchors, though only one
  is needed.
- **The other 17** have an anchor that runs past the part shared with the
  canonical: the shared stretch at that end is shorter than 10 residues, or
  the peptide itself is.

## 5. Rules compared (anchor length 10)

| Rule | Correct | Wrong | Missed | Unjoined (correct) |
|---|---|---|---|---|
| first occurrence (current) | 880 | 12 | 42 | 64 |
| unique: each anchor exactly once, start not after end (MED-3) | 880 | 10 | 44 | 64 |
| consistent: as unique, but a repeated anchor is resolved by the pair whose canonical span best fits the peptide length; a tie is not joined | 882 | 10 | 42 | 64 |
| first + complete ends | 915 | 2 | 17 | 64 |
| unique + complete ends | 915 | 0 | 19 | 64 |
| **consistent + complete ends** | **917** | **0** | **17** | 64 |

"Complete ends": a peptide ending at a stop codon is not joined at its
C-terminus, and one starting at the gene's start codon is not joined at its
N-terminus. A peptide complete at both ends is the whole protein.

**Changes from current to proposed:** 25 missed → correct, and
12 wrong → correct. Nothing gets worse.

## 6. Anchor length

### 6.1 chr15

Full table: [sweep.tsv](stitch_audit/sweep.tsv).

| k | Current: correct / wrong / missed | Proposed: correct / wrong / missed |
|---|---|---|
| 4 | 891 / 21 / 22 | 926 / 0 / 8 |
| 5 | 899 / 13 / 22 | 933 / 0 / 1 |
| 6 | 898 / 13 / 23 | 933 / 0 / 1 |
| 7 | 895 / 13 / 26 | 933 / 0 / 1 |
| 8 | 889 / 12 / 33 | 927 / 0 / 7 |
| 10 (current) | 880 / 12 / 42 | 917 / 0 / 17 |
| 15 | 858 / 11 / 65 | 898 / 0 / 36 |
| 20 | 826 / 10 / 98 | 867 / 0 / 67 |

Longer anchors only lose joins. Under the proposed rule, the one miss left
at k = 5–7 is A3SS 1317 NPTN slice 2: its shared C-terminal stretch is
shorter than 5 residues.

### 6.2 Proteome (GENCODE v50 MANE Select)

[survey_anchors.py](stitch_audit/survey_anchors.py) uses the 19,168 MANE
Select proteins, translated from `gencode.v50.transcripts.fa.gz` with the
CDS from `gencode.v50.annotation.gtf`. 61 that do not translate cleanly
(non-ATG start, internal stop codon, selenoproteins) are excluded.
Anchors are taken where a slice's would be: the first k residues after,
and the last k before, each of the 172,210 internal exon boundaries.

- **Repeated:** the anchor occurs more than once in its own protein. The
  proposed rule then resolves it by the other anchor, or does not join.
- **Chance match:** an anchor from a randomly paired other protein (fixed
  seed) occurs in this protein. This is the rate at which an anchor not in
  the canonical would be found anyway, i.e. a wrong join.
- An earlier version paired each protein with the next one in the file.
  Neighbours are often paralogs, which inflated chance matches about
  100-fold, so it was replaced.

| k | Anchors | Repeated in own protein | Chance match |
|---|---|---|---|
| 4 | 343,351 | 3.94% | 0.78% |
| 5 | 343,056 | 1.36% | 0.055% |
| 6 | 342,704 | 0.88% | 0.011% |
| **7** | 342,402 | **0.71%** | **0.007%** |
| 8 | 342,102 | 0.60% | 0.007% |
| 10 (current) | 341,426 | 0.51% | 0.005% |
| 15 | 339,638 | 0.43% | 0.004% |
| 20 | 337,433 | 0.40% | 0.003% |

Both rates fall steeply up to k = 6–7 and then level off. Beyond 7, the
repeats are low-complexity or internal duplications that no practical
length removes. Longer anchors only miss slices with short shared ends.
**k = 7** has the lowest error rates of the lengths that lose nothing on
chr15.

### 6.3 UniProt canonicals

A second run takes the canonical from the local UniProt cache
(`cache/uniprot-cache.db`; 146 of the chr15 genes, 776 distinct join
attempts), as a default run would. UniProt sequences have no coordinates,
so only anchor multiplicity is measured. No anchor occurs more than once.
The rates of missing anchors are similar to the GTF canonical (start
anchor not found 58, end anchor not found 235). See
[summary.txt](stitch_audit/summary.txt).

## 7. Proposed change (PR 5)

1. **Complete ends are not joined.**
   - A peptide ending at a stop codon is joined at its N-terminus only
     (T4 force-translated slices, and T1 slices ending at a stop in their
     last codon).
   - A peptide starting at the gene's start codon is joined at its
     C-terminus only.
   - A peptide complete at both ends is written as it is.
   - `_join` (PR 4c) already chooses `n`/`c` by anchors and needs no
     change here.
2. **Anchor choice (MED-3).** Each anchor is used only if it occurs exactly
   once.
   - If either occurs more than once, the pair of occurrences (start not
     after end) whose canonical span is closest to the peptide length is
     used. A tie is not joined.
   - The start anchor must not be after the end anchor.
   - The same rule applies in `stitch_to_canonical_aa` and `_join`.
3. **Anchor length 7** (`params.aa_stitch_length`). The existing assertion
   (≥ 6) still holds.
4. **Tests:**
   - the 8 T4 cases stay truncated (e.g. UBE2Q2 SE 19972);
   - COMMD4 SE 5917 is not given the canonical N-terminus;
   - CD276 A3SS 1316 is joined at the second occurrence;
   - a unit test on a synthetic repeated anchor;
   - a slice shorter than 10 but at least 7 residues is joined (GTF2A2
     SE 11950).
5. **Verification:** re-run this audit on the working tree. It must give
   the proposed counts. Every changed FASTA record is listed with its exact
   annotated-protein matches. The golden test is updated with the changes
   explained.

**Expected output change:**
- **No junction changes fate.** Tiers do not depend on joining, except the
  file chosen for an alternative-termini slice; see decision 2.
- **Records move or change** between the joined files, the `_orphan` files
  and `psq_canonical`.

## 8. Decisions requested

**Approved 2026-11-11:**
1. "Unique + complete ends". The stricter rule is preferred: leaving two
   slices as orphans costs less than writing a wrong joined record. CD276
   A3SS 1316 therefore stays unjoined.
2. Anchor length 7, with a separate terminus-comparison length kept at 10.
3. No command-line option.

1. **Joining rule:** "consistent + complete ends" (recommended), or
   "unique + complete ends"? The latter leaves CD276's two slices unjoined
   (orphans) instead of joining them correctly.
2. **Anchor length:**
   - (a) **Recommended:** 7.
   - (b) Keep 10 and fix only the rule (17 missed joins remain on chr15).
   - **Side effect of (a):** PR 4c also uses `aa_stitch_length` to decide
     whether an alternative-termini product differs from the canonical at
     its N- or C-terminus (comparing the first and last k residues), and so
     which `psq_T1_alt_*` file it goes to.
   - I propose a separate parameter for that comparison, kept at 10, so
     that PR 5 does not move any alternative-termini record between files.
     The other choice is to let it follow the anchor length, with the
     effect measured in verification.
3. **Exposing the anchor length:** should it be a command-line option?
   Today it is only in `params.py`. My recommendation is no: keep it in
   `params.py` and document the measured choice.

## 9. Implementation and verification (PR 5)

Implemented as approved in §8:

- **One method for every join.** `Sequence._join` now serves both
  `stitch_to_canonical_aa` (T1, T2, T4) and the alternative-termini path.
  - It uses an anchor only if the anchor occurs once.
  - A join at both ends needs the start anchor not to be after the end
    anchor.
  - Given the slice and its translation start, it does not join an end
    that is already complete (`Sequence.ends_at_stop`; the gene's start
    codon).
- **Anchor length:** `params.aa_stitch_length` is 7.
- **Terminus comparison:** the alternative-termini path compares termini
  over the new `params.aa_terminus_length`, kept at 10.

**Verification** ([verify_5.py](stitch_audit/verify_5.py); evidence in
[stitch_audit/final/](stitch_audit/final/)):

| Check | Result |
|---|---|
| Pipeline joins differing from the audited rule (unique + complete ends, k = 7) | 0 of 994 distinct |
| Pipeline joins against the genomic truth | 927 correct, **0 wrong**, 3 missed (CD276 A3SS 1316 × 2 by decision, NPTN A3SS 1317), 64 unjoined (correct) |
| Junctions whose fate changed | 1: GTF2A2 SE 11945, T1 ALT TERMINUS → T1 |
| Annotated proteins in some output file before but in none now | 1: CD276-201 |
| Annotated proteins in some output file now but in none before | 7 |
| Golden test | updated (below) |

- **The 994 distinct joins (998 in the audit)** differ because some
  products now differ, which changes which calls are distinct.
- **GTF2A2 SE 11945:** its slice-2 peptide is 7 residues long.
  - With a 10-residue anchor it could not be joined, so it was written to
    `psq_T1_alt_cterm_orphan` and the junction counted as ALT TERMINUS.
  - Joined now, it gives exactly the canonical (109 aa), as does slice 1:
    an internal event ending at the canonical's stop codon is plain T1.
- **CD276-201** (316 aa, the 2Ig isoform) was produced only by the wrong
  join of §4. The first occurrence of the repeated anchor happened to give
  an annotated isoform. The junction does not support that isoform, and
  the genomic truth places the slice in the other copy of the domain.
  Under the approved rule the slice is not joined (§8, decision 1).
- **The 7 annotated proteins now written include RGMA-201, GALK2-201,
  BLOC1S6-201, EMC4-201 and RPS27L-201.** These are T4 products that are
  annotated protein-coding isoforms ending at their own stop codon. Before,
  they were joined C-terminally into the canonical.

**FASTA records** (baseline `f1db483` → PR 5):

| File | Before | After |
|---|---|---|
| `psq_T1` | 173 | 179 |
| `psq_T1_alt_cterm` | 18 | 18 |
| `psq_T1_alt_cterm_orphan` | 2 | 0 |
| `psq_T1_alt_nterm` | 4 | 5 |
| `psq_T1_orphan` | 35 | 28 |
| `psq_T2` | 23 | 32 |
| `psq_T2_orphan` | 14 | 4 |
| `psq_T4` | 226 | 251 |
| `psq_T4_orphan` | 76 | 61 |
| `psq_canonical` | 268 | 268 |

Records move from the orphan files to the joined files. Every changed
record, with its annotated matches, is in
[final/fasta_diff.tsv](stitch_audit/final/fasta_diff.tsv).

**Regression tests** ([tests/test_stitch.py](../tests/test_stitch.py)). All
five fail at `f1db483`.

| Test | Asserts |
|---|---|
| `test_peptide_ending_at_a_stop_is_not_joined_c_terminally` | UBE2Q2 SE 19972: the T4 record is the 197-aa truncated product, not the canonical |
| `test_peptide_starting_at_the_gene_start_is_not_joined_n_terminally` | COMMD4 SE 5917: the product starts MSSVKL (175 aa) |
| `test_repeated_anchor_is_not_joined` | CD276 A3SS 1316: slice 1 is an orphan |
| `test_join_rule_on_synthetic_canonical` | unique anchors join; a repeated anchor or a start anchor after the end anchor does not |
| `test_short_peptide_is_joined` | GTF2A2 SE 11950: the 7-residue peptide is joined (giving the canonical), not orphaned |
