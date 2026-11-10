# Audit: translation-bound (`tx0`/`tx1`) defect in `Junction._get_translated_region`

| | |
|---|---|
| **Date** | 2026-11-03 |
| **Status** | Pre-fix audit. No production code changed. Input to PR 4 of the [v0.4.0 roadmap](../REVIEW.md#9-v040-hardening-roadmap-revised). |
| **Code audited** | JCAST v0.3.5a, commit `3ed17b0` (branch `v0.4.0-hardening`), [jcast/junctions.py:274-357](../jcast/junctions.py#L274-L357) |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions, Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 (corrected 2026-11-10; the versions first listed here were wrong. Re-running in this environment reproduces all three evidence files byte for byte.) |
| **Reproduce** | `python reports/tx1_audit/audit_tx1.py` (≈ 2 min, no network). See [§9](#9-reproducing-this-audit). |
| **Evidence** | [tx1_audit/summary.txt](tx1_audit/summary.txt) · [tx1_audit/junctions.tsv](tx1_audit/junctions.tsv) (one row per junction) · [tx1_audit/fasta_diff.tsv](tx1_audit/fasta_diff.tsv) (every changed output sequence) |

Every number in this report is printed by the audit script or can be read
directly from `junctions.tsv`.

---

## 1. Summary

**The defect.** JCAST trims each splice event to the annotated coding region
before translating it. The lower genomic bound is `tx0` and the upper bound is
`tx1`. The bound that should mark the **end of translation (the stop codon)**
is instead computed from **start codons**:

- On the `+` strand, `tx1` is the most downstream start codon end + 2.
- On the `−` strand, both bounds end up on the same start codon.

Across all 1,105 junctions, the C-terminal bound coincides with an annotated
stop codon **0 times**.

**Consequences, measured on chr15:**

1. **C-terminal splice events are effectively untranslatable.** 63 junctions
   span the stop-codon exon. Excluding the 19 skipped for read count, the
   shipped code translates **2 of 44** (1 T2, 1 T4) and fails the other 42.
   Their slices run into the 3′ UTR, so the real stop codon appears mid-slice
   and is rejected as premature. With correct bounds, **36 of 44** translate
   (7 T1, 4 T2, 25 T4).
2. **Spurious trims at start codons create artefactual sequences.** In 129
   junctions the misplaced bound fires inside a start-codon exon. It either
   cuts that exon a few nucleotides after the ATG or, on the `−` strand,
   keeps the **5′-UTR side** of the exon and discards the coding side. Some
   of these produce "successful" T1/T2 output that is really translated
   5′ UTR, or an "isoform" from which the splice event itself was deleted.
3. **The other 917 junctions are unaffected.** Their fates are identical under
   both versions, so the fix is surgical.

**Effect of a corrected implementation on chr15:**

| | Shipped | Corrected | Δ |
|---|---|---|---|
| T1 / T2 / T3 / T4 / FAIL / SKIP | 185 / 25 / 2 / 236 / 360 / 297 | 189 / 26 / 2 / 261 / 330 / 297 | 48 junctions change fate |
| `psq_T1.fasta` records | 152 | 164 | **+12 gained, 0 lost** |
| T1 records exactly matching an annotated Ensembl protein | 58 | 66 | +8 (all from the 12 gained) |
| Sequences that match an annotated protein and would be removed | — | — | **1** (a 43-aa BLOC1S6 orphan fragment, §5.5) |

**For users of ≤ v0.3.5 outputs.** On this data, the fix removes **no**
sequence from `psq_T1.fasta` or `psq_canonical.fasta`; it only adds to them.
The artefactual sequences it removes sit in the lower-confidence files (T2,
T4 and orphans; §5.4). Published T1 + canonical databases are therefore
**incomplete for C-terminal events but not contaminated by this defect**.
Lower-tier files may contain a small number of 5′-UTR-derived artefacts.

**Recommendation.** Proceed with PR 4 using the bound definition in §7, which
is the "corrected" variant audited here. Land the regression tests and the
golden-test update in the same commit, plus the CHANGELOG text in §7.4. Track
the residual issues in §8 separately; they exist under both versions.

---

## 2. Background: how the bounds are used

[`Junction.trim_cds`](../jcast/junctions.py#L427) calls
`_get_translated_region` to set `tx0`/`tx1`, then
[`_trim`](../jcast/junctions.py#L454) shaves exon coordinates in genomic
orientation, independent of strand:

- An exon that *contains* `tx0` has its start moved up to `tx0`. If that exon
  is an alternative exon, the exon below it is dropped (set to −1).
- An exon that *contains* `tx1` has its end moved down to `tx1`. If that exon
  is an alternative exon, the exon above it is dropped.

`_trim` therefore expects `tx0` and `tx1` to be the **lowest and highest
genomic coordinates of the coding sequence**:

| Strand | `tx0` (lower genomic bound) should be | `tx1` (upper genomic bound) should be |
|---|---|---|
| `+` | start codon (5′ end of CDS) | **stop codon end** (3′ end of CDS) |
| `−` | **stop codon start** (3′ end of CDS) | start codon end (5′ end of CDS) |

The bound that marks the C-terminus is `tx1` on `+` and `tx0` on `−`.

Trimming at the stop codon matters for translation:
[`make_pep(..., terminate=True)`](../jcast/helpers.py#L157) returns an empty
peptide if it meets a stop codon anywhere except the last codon. A slice that
runs past the stop into the 3′ UTR is therefore classified as containing a
premature termination codon (PTC).

---

## 3. The defect

### 3.1 Code

[jcast/junctions.py:296-348](../jcast/junctions.py#L296-L348), abridged:

```python
# start codons → tx0
gtf0_start = gtf0.query('feature == "start_codon" & ...').loc[:, 'start']
self.tx0 = sorted(gtf0_start, reverse=self.strand == '-')[startsite_index]

# "2020-07-25 now getting end codons ..."  ← comment says end codons
gtf0_end = gtf0.query('feature == "start_codon" & ...').loc[:, 'end']   # ← queries start codons
self.tx1 = sorted(gtf0_end, reverse=self.strand == '+')[0]

if self.strand == '-' and self.tx1 > 0 and self.tx0 > 0:
    self.tx0, self.tx1 = self.tx1, self.tx0
    self.tx1 += 2
elif self.strand == '+':
    self.tx1 += 2
```

The comment and the variable name both say "end codon", but the query selects
`start_codon` rows. It looks like a copy-paste of the block above it.

### 3.2 What it actually computes

Let *S* be the set of annotated start codons (protein-coding, TSL ≤ 2) for the
gene.

| Strand | `tx0` computed | `tx1` computed | Intended C-terminal bound |
|---|---|---|---|
| `+` | min(*S*.start) ✔ | max(*S*.end) + 2, i.e. **2 nt past the most downstream start codon** ✘ | max(stop_codon.end) |
| `−` | min(*S*.end), i.e. **a start codon** ✘ | max(*S*.start) + 2 = start codon end ✔ | min(stop_codon.start) |

When a gene has a single start codon on the `−` strand, both bounds collapse
onto it (`tx0 == tx1`). This is visible throughout `junctions.tsv`. For
example SEC11A gets `current_tx0 = current_tx1 = 84691617`, while its stop
codon starts at 84,670,019.

**Measured** (`junctions.tsv`; `+` n = 558, `−` n = 547):

- The C-terminal bound equals the corrected (stop-codon) bound in **0 %** of
  junctions on either strand.
- The 5′ bound matches in 77 % (`+`) and 79 % (`−`). The remainder differ
  only because `trim_cds` retries with the *second* start codon when the
  trimmed anchor is not a CDS exon, and that retry is triggered by the
  spurious trims (§8.3). The 5′-bound formula itself is the same in both
  versions.
- The misplaced C-terminal bound sits a median **28 kb** (range 431 nt to
  42.7 Mb) from the real one. It lies near the start of the gene, not its end.

### 3.3 Two failure modes

**Mode A — no trimming at the stop codon** (junction spans the C-terminal
exon). Example: **MFGE8, SE 5739, `−` strand.** The anchor exon
88,898,683–88,899,532 contains the stop codon (starts at 88,899,395).

| | `tx0` | Anchor kept | Slice 1 length | First in-frame stop | Fate |
|---|---|---|---|---|---|
| Shipped | 88,913,319 (start codon) | 88,898,683–88,899,532 (incl. 3′ UTR) | 1,191 nt | nt 476 → PTC | **FAILURE** |
| Corrected | 88,899,395 (stop codon) | 88,899,395–88,899,532 | 479 nt | nt 476 = final codon | **T1** |

**Mode B — a spurious trim at a start codon** (junction contains the
start-codon exon). Example: **FAM219B, SE 5946, `−` strand.** The skipped exon
74,906,278–74,906,365 contains the start codon (end 74,906,321). Upstream in
transcript terms means higher coordinates, so the region *above* 74,906,321 is
5′ UTR.

| | Bounds | Exons kept | What was translated | Fate |
|---|---|---|---|---|
| Shipped | `tx0 = tx1 = 74,906,321` | alt1 → **74,906,321–74,906,365 (5′-UTR side)**; downstream (5′-UTR exon) kept; anchor (coding) **dropped** | 5′ UTR, 3-frame fallback | **T1** (artefact) |
| Corrected | 74,902,619 – 74,906,321 | alt1 → 74,906,278–74,906,321 (coding side); 5′-UTR exon dropped; anchor kept | Coding sequence | FAILURE (§5.3) |

On the `+` strand, Mode B instead cuts the start-codon exon 4 nt after the
ATG. For example EMC4 (start codon 34,225,115) gets `tx1` = 34,225,119.

---

## 4. Method

**Harness.** [audit_tx1.py](tx1_audit/audit_tx1.py) runs every junction in
`tests/data/rmats` through the shipped per-junction pipeline
([`main._translate_one`](../jcast/main.py#L140)) twice:

1. **Shipped:** `_get_translated_region` as committed.
2. **Corrected:** `_get_translated_region` replaced at runtime by
   `corrected_get_translated_region` (defined in the script, specification in
   §7.1). Nothing else is changed.

Settings match the golden test (`-r 1`, `-c`, `--mask NnXx`, `-q 0 1`).
Canonical proteins come from the GTF (`params.use_gtf_only = True`) and
UniProt is stubbed out. The run is therefore deterministic and offline. Tier
assignment does not depend on the canonical source; only the routing between
stitched and orphan files does.

The shipped run reproduces the committed golden test exactly
([tests/test_pipeline.py](../tests/test_pipeline.py): FAIL 360, T4 236,
SKIP 297, T1 185, T2 25, T3 2; identical FASTA digests). This confirms the
harness matches production behaviour.

**Independent truth set.** Every transcript with CDS rows in the GTF (any
biotype, any TSL) is translated directly from the genome with Biopython,
giving 2,979 distinct proteins. A JCAST output record "matches an annotated
protein" when its full sequence equals one of them. As a sanity check, all
268/268 canonical records match.

**Audit-trail note.** An earlier prototype of the corrected logic, used for
the figures quoted in REVIEW.md on the same date, silently dropped the
TSL ≤ 2 filter. pandas resolves `@tsl` in `DataFrame.query` from the calling
frame's locals, and the prototype evaluated it inside a nested function. The
script here keeps `tsl` local to the querying frame. REVIEW.md §5.1 has been
updated to these final figures.

**Limitations.** This is one chromosome and one tissue's rMATS output, and
these are relative effects. The truth set measures agreement with Ensembl
annotation, not biological validity: a novel isoform will legitimately fail
to match. UniProt-canonical runs (the default) will route some sequences
differently between stitched and orphan files.

---

## 5. Results

### 5.1 Trimming activity

From `summary.txt`: junctions in which the lower or upper bound changed any
coordinate.

| Type | n | Shipped lower / upper | Corrected lower / upper |
|---|---|---|---|
| MXE | 90 | 28 / 10 | 25 / 12 |
| SE | 720 | 133 / 109 | 106 / 84 |
| RI | 88 | 2 / 2 | 6 / 4 |
| A5SS | 71 | 9 / 4 | 5 / 1 |
| A3SS | 136 | 13 / 2 | 8 / 2 |

188 junctions end up with different trimmed coordinates.

### 5.2 Fate changes

48 of 1,105 junctions (4.3 %) change fate. No junction moves into or out of
SKIPPED (297 under both).

| Shipped → Corrected | n | Direction |
|---|---|---|
| FAILURE → T4 | 26 | better |
| FAILURE → T1 | 8 | better |
| FAILURE → T2 | 3 | better |
| T4 → T2 | 1 | better |
| T4 → FAILURE | 3 | worse |
| T2 → FAILURE | 2 | worse |
| T1 → FAILURE | 2 | worse |
| T1 → T4 | 2 | worse |
| T2 → T4 | 1 | worse |

### 5.3 Every change is explained by the two failure modes

Each junction is classed by whether the corrected C-terminal trim fires
(**A**: the event spans the stop-codon exon), and whether the shipped code's
misplaced C-terminal bound fires (**B**: a spurious trim at a start codon).

| Class | n | Shipped T1/T2/T4/FAIL/SKIP | Corrected T1/T2/T4/FAIL/SKIP | Fate changes |
|---|---|---|---|---|
| A only | 59 | 0 / 0 / 1 / 40 / 18 | 7 / 4 / 25 / 5 / 18 | 36 |
| B only | 125 | 11 / 3 / 15 / 65 / 31 | 8 / 1 / 16 / 69 / 31 | 11 |
| A and B | 4 | 0 / 1 / 0 / 2 / 1 | 0 / 0 / 0 / 3 / 1 | 1 |
| **Neither** | **917** | 174 / 21 / 220 / 253 / 247 | **identical** | **0** |

(T3 = 2 in "Neither" under both versions; omitted for width.)

**Improvements (38).** 36 are class A. The slice previously ran past the stop
codon into the 3′ UTR. With the stop codon as the final codon, the junction
now translates. The other 2 are class B: DET1 A5SS 748 and MEIS2 A3SS 3449,
where the spurious start-codon trim had broken an otherwise translatable
junction. In the T4 cases, one slice still carries a genuine
in-frame PTC: the junction's alternative form introduces a stop, the usual
signature of an NMD-targeted isoform. Representative cases:

| Junction | Strand | Shipped: slice length, first stop | Corrected: slice length, first stop | Fate |
|---|---|---|---|---|
| MFGE8 SE 5739 | − | 1,191 nt, stop at 476 | 479 nt, stop at 476 (final) | FAIL → T1 |
| CASC4 SE 6426 | + | 2,762 nt, stop at 407 | 410 nt, stop at 407 (final) | FAIL → T1 |
| COMMD4 SE 5921 | + | 367 nt, stop at 296 | 299 nt, stop at 296 (final) | FAIL → T1 |
| CHSY1 SE 5606 | − | 4,074 nt, stop at 2,406 | 2,409 nt (final); slice 2 PTC at 462 | FAIL → T4 |

**Worse outcomes (10).** 9 are class B and 1 (OAZ2) is class A+B. In every case the
shipped "success" was an artefact. Traced at nucleotide level:

| Junction | Strand | Shipped → Corrected | What the shipped code had actually translated |
|---|---|---|---|
| FAM219B SE 5946 | − | T1 → FAIL | 5′ UTR (worked example, §3.3) |
| MTMR10 SE 6602 | − | T1 → T4 | 5′ UTR (anchor dropped, start exon trimmed to UTR side). Corrected slice 2 has a genuine PTC. |
| OAZ2 SE 24775 | − | T2 → FAIL | 5′ UTR. OAZ2 (antizyme 2) needs a programmed +1 ribosomal frameshift, so a frame-based translator *should* fail. |
| SPG21 MXE 4040 | − | T2 → T4 | An 11-nt sliver of the start exon plus 5′-UTR exons. Its 282-aa annotated protein is still output, now in T4 (§5.5). |
| WDR73 A5SS 3645 | − | T1 → FAIL | The alternative exons were cut to 7 nt, so **both slices were identical** (75 nt): a "T1 isoform pair" with no splice difference. |
| BLOC1S6 SE 6398, 6399, 6401 | + | T4 → FAIL | Downstream exon dropped at an internal start codon, so the skip form became "anchor only". Corrected: both forms carry in-frame stops. |
| BLOC1S6 SE 6400 | + | T1 → T4 | Same as above: slice 2 was the anchor exon alone. |
| ZNF592 SE 5769 | + | T2 → FAIL | Upstream exons are entirely 5′ UTR (the start codon is in the downstream exon). Shipped code translated them. The corrected code fails because `_trim` never removes whole UTR exons (§8.1). |

### 5.4 Output files

From `summary.txt` and `fasta_diff.tsv`. "Lost" and "gained" are sequences
present under only one version.

| File | Shipped | Corrected | Lost | Gained | Gained that match an annotated protein |
|---|---|---|---|---|---|
| `psq_canonical.fasta` | 268 | 268 | 0 | 0 | — |
| `psq_T1.fasta` | 152 | 164 | **0** | 12 | 8 |
| `psq_T1_orphan.fasta` | 48 | 46 | 3 | 1 | 0 |
| `psq_T2.fasta` | 20 | 23 | 2 | 5 | 1 |
| `psq_T2_orphan.fasta` | 13 | 11 | 6 | 4 | 0 |
| `psq_T3_orphan.fasta` | 4 | 4 | 0 | 0 | — |
| `psq_T4.fasta` | 193 | 216 | 2 | 25 | 12 |
| `psq_T4_orphan.fasta` | 72 | 81 | 7 | 16 | 1 |

Records exactly matching an annotated Ensembl protein:

| File | Shipped | Corrected |
|---|---|---|
| `psq_T1.fasta` | 58 / 152 (38 %) | 66 / 164 (40 %) |
| `psq_T2.fasta` | 3 / 20 | 3 / 23 |
| `psq_T4.fasta` | 59 / 193 | 70 / 216 |

The 4 gained T1 sequences that do *not* match an annotation are candidates
for unannotated isoforms. They stitch to the canonical protein at both ends,
the same evidence standard as every other T1 entry. They are listed in
`fasta_diff.tsv` for manual review.

### 5.5 Sequences matching an annotated protein that the fix removes

Three "lost" records match an annotated protein. Two are only reclassified:

- **APH1B MXE 3368:** the same 217-aa sequence moves from T4 to T2. Both
  slices now translate; the length difference is not a multiple of 3, hence
  T2.
- **SPG21 MXE 4040:** the same 282-aa sequence moves from T2 to T4.

One is genuinely no longer produced:

- **BLOC1S6 SE 6401:** `MSVPGPSSPDGALTRPPYCLEAGEPTPEPSCIVRHTGTRDFKI`
  (43 aa), formerly in `psq_T4_orphan.fasta`. It matches a short annotated
  BLOC1S6 protein, but the shipped code produced it from an artefactually
  truncated slice (the downstream exon was dropped at an internal start
  codon), not from this junction's splice event. Under the corrected bounds,
  both forms of the junction carry in-frame stops (FAIL), which is the
  correct reading of this event.

---

## 6. Impact assessment

| Who | Effect | Action |
|---|---|---|
| Users of `psq_T1.fasta` + `psq_canonical.fasta` (the recommended search database) | Missing isoforms for splice events at the C-terminus. On chr15 this is 12 additional T1 sequences (+8 %). No contamination found from this defect. | Re-run with v0.4.0 to recover C-terminal isoforms. Earlier identifications remain valid. |
| Users who also searched T2–T4 or orphan files | As above, plus a small number of artefactual sequences: 5′-UTR translations, and "isoforms" missing the splice event. On chr15, 20 records leave the T2/T4/orphan files (2 of them, APH1B and SPG21, reappear in another tier; §5.5). | Re-run. Treat earlier hits that map only to removed records (listed in `fasta_diff.tsv`, `current_only`) as unsupported. |
| `-m` read-count model | Unaffected. Junction counts are computed before trimming. | None. |

The share of junctions in class A (5.7 %) and B (11.7 %) on chr15 is a rough
guide for genome-wide runs. The true proportion depends on how many
detected events fall in terminal exons.

---

## 7. Fix specification for PR 4

### 7.1 Bound definition

With *S* and *T* = annotated start and stop codon rows (protein-coding,
TSL ≤ `params.tsl_threshold` when that column exists):

| Strand | `tx0` | `tx1` |
|---|---|---|
| `+` | the `startsite_index`-th smallest *S*.start | max(*T*.end) |
| `−` | min(*T*.start) | the `startsite_index`-th largest *S*.end |

- No `+2` offsets. GTF codon coordinates are inclusive, so these values
  already include the full codon. The old `+2` would leave 2 nt past the stop,
  which `make_pep` treats as a PTC (REVIEW.md §4, LOW).
- Each bound independently becomes `None` (no trimming on that side) when it
  cannot be determined. On chr15 that is 51 `+` and 23 `−` junctions in
  genes with no stop codon at TSL ≤ 2.
- `num_start_codons` keeps its meaning (distinct start positions), so the
  `trim_cds` second-start-site retry behaves as before.
- TSL handling: check `'transcript_support_level' in gtf.annot.columns`
  instead of catching exceptions. This also fixes REVIEW.md MED-7 (the `None`
  → `len()` path).
- Choosing the most upstream start and most downstream stop across transcripts
  gives the widest annotated CDS envelope for the gene. This is deliberate,
  and matches the existing "longest transcript" intent for the start codon.
  Its limitation is in §8.2.

### 7.2 Regression tests (land before the code change; they must fail first)

| Test | Junction | Asserts |
|---|---|---|
| `+` strand C-terminal bound | CASC4, SE 6426 | `tx0 == 44289030` (start codon), `tx1 == 44413406` (stop codon end), and after trimming `down_ee == 44413406` |
| `−` strand C-terminal bound | MFGE8, SE 5739 | `tx0 == 88899395` (stop codon start), `tx1 == 88913319` (start codon end), and after trimming `anc_es == 88899395` |
| C-terminal event translates | MFGE8, SE 5739 | fate is T1 |
| Existing `test_trim_start` | SEMA4B, A5SS 742 | `tx0 == 90201579`. Unchanged under the corrected bounds (verified). |

### 7.3 Golden test update

`tests/test_pipeline.py` must change in the same commit:

- `EXPECTED_FATES`: FAILURE 330, PARTIAL TIER 4 261, SKIPPED 297,
  SUCCESS TIER 1 189, SUCCESS TIER 2 26, SUCCESS TIER 3 2.
- `EXPECTED_FASTA`: regenerate. Record counts must equal the "Corrected"
  column of §5.4.

If the implementation in PR 4 reproduces these numbers exactly, it matches
the variant audited here.

### 7.4 Proposed CHANGELOG entry

> **Changed (affects output).** The translation-end bound used to trim splice
> events to the coding region was computed from start codons instead of stop
> codons (since 2020-07). As a result, splice events spanning a protein's
> C-terminal exon were almost never translated, and some events containing a
> start-codon exon yielded 5′-UTR-derived or event-less sequences in lower
> tiers. The bound now uses annotated stop codons. On the bundled chr15
> data, 48 of 1,105 junctions change tier: `psq_T1.fasta` gains 12 sequences
> and loses none, and `psq_canonical.fasta` is unchanged. Re-running is
> recommended to recover C-terminal isoforms. See
> `reports/2026-11-03-tx1-translation-bounds-audit.md`.

---

## 8. Residual issues (out of scope for PR 4; present under both versions)

1. **Whole UTR exons are never removed.** `_trim` only adjusts an exon that
   *contains* a bound. An exon lying entirely outside [`tx0`, `tx1`] is kept.
   Under the corrected bounds, 172 junctions keep at least one pure-UTR exon
   (121 of them FAIL). ZNF592 (§5.3) is an example. Suggested follow-up:
   drop exons that fall entirely outside the bounds, with its own audit and
   golden update.
2. **One bound for all transcripts.** The bounds are the gene's widest CDS
   envelope. A junction belonging to a transcript with an earlier stop, or a
   later start, is trimmed to the envelope rather than to its own CDS, which
   produces PTCs (T4/FAIL). The proper fix is per-transcript assignment
   (REVIEW.md, v1.0 splice-graph item).
3. **Second-start-site retry.** `trim_cds` retries with the second start
   codon when the trimmed anchor is not an annotated CDS exon. Spurious trims
   in the shipped code trigger this retry. It should be revisited after
   PR 4, which makes it more meaningful.
4. **Phase lookup after trimming.** `get_translated_phase` takes the frame
   from an exon chosen by type and strand (e.g. `down` for `−`-strand SE). If
   trimming has correctly dropped that exon as UTR, the phase is `None` and
   translation falls back to the 3-frame search (FAM219B, OAZ2). The phase
   should come from the first *retained* coding exon.
5. **The three-frame fallback keeps the last frame** (REVIEW.md MED-2). It
   interacts with item 4.

---

## 9. Reproducing this audit

From the repository root, in an environment with JCAST's dependencies:

```bash
python reports/tx1_audit/audit_tx1.py
```

The script regenerates `summary.txt`, `junctions.tsv` and `fasta_diff.tsv`
in `reports/tx1_audit/`. It creates and deletes its own temp directory, needs
no network, and leaves the working tree unchanged apart from those three
files.

**Input checksums** (SHA-256, first 16 hex digits):

| File | SHA-256 |
|---|---|
| `tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf` | `4d77058667240293` |
| `tests/data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz` | `7f17face433949c4` |
| `tests/data/rmats/A3SS.MATS.JC.txt` | `28274cedf4939dd5` |
| `tests/data/rmats/A5SS.MATS.JC.txt` | `f6c1980e6ed62b2f` |
| `tests/data/rmats/MXE.MATS.JC.txt` | `ea3877af20b96208` |
| `tests/data/rmats/RI.MATS.JC.txt` | `fc35b607c9389aee` |
| `tests/data/rmats/SE.MATS.JC.txt` | `5f5e826f717653f7` |

**`junctions.tsv` columns.** `input_*` are rMATS coordinates converted to
1-based inclusive. `current_*` and `corrected_*` are the bounds, trimming
flags, GTF phase, fate and trimmed coordinates under each version. `cds_min`
and `cds_max` are the gene's annotated CDS envelope (incl. stop codon,
TSL ≤ 2). `current_tx_*_inside_cds` flags bounds strictly inside an annotated
CDS exon. `*_exon_dropped` flags a trim that removed a whole exon; this is
an artefact under the shipped bounds but can be legitimate UTR removal under
the corrected ones. `corrected_utr_exon_kept` marks residual issue §8.1.

**`fasta_diff.tsv` columns.** `file`, `present_in` (`current_only` /
`corrected_only`), `matches_annotated_protein`, `header`, `sequence`.
