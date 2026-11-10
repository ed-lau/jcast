# Audit: trimming splice events to the coding region (PR 4a)

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Audit of the PR 4a implementation, run before it is committed. Baseline is PR 4 (commit `da72a0a`); candidate is the PR 4a working tree. |
| **Code audited** | `Junction._trim`, `Junction.trim_cds`, and the new `Junction.event_outside_cds` ([jcast/junctions.py](../jcast/junctions.py)); the new skip in `main._translate_one` |
| **Data** | Bundled chr15 test set: 1,105 rMATS junctions, Ensembl GRCh38.89 GTF, GRCh38 chr15 DNA |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/utr_audit/audit_utr.py` (≈ 2.5 min, no network). See [§8](#8-reproducing-this-audit). |
| **Evidence** | [utr_audit/summary.txt](utr_audit/summary.txt) · [utr_audit/junctions.tsv](utr_audit/junctions.tsv) (one row per junction) · [utr_audit/fasta_diff.tsv](utr_audit/fasta_diff.tsv) (every changed output sequence) |

Every number in this report is printed in `summary.txt` or can be read
directly from the two TSV files. The only exception is §4.3, which describes
an intermediate candidate that kept the retry; it was run with the same
script before the retry was removed.

---

## 1. Summary

**Decision being implemented.** Exons lying entirely in the UTR are removed
before translation, because translating UTR is outside JCAST's purpose. This
was item 1 of the
[tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions).

**Defects found in `_trim`.** It used a separate hand-written rule for each
junction type. Besides never removing whole UTR exons, this had two more
defects:

1. **172 junctions kept a whole exon outside the coding region** after PR 4.
   Most of these failed. A few produced "successes" translated from 5′ UTR.
2. **A5SS and A3SS on the `−` strand were trimmed at the wrong bound.** The
   rules assumed `+`-strand exon order. On `−`, a bound inside an alternative
   exon or the anchor was never applied. 11 junctions are now trimmed
   differently.
3. **The retained region of RI events was never trimmed.** Only the flanking
   exons were (10 junctions).

**Change.**

- **One strand-agnostic rule for every exon.** An exon containing a bound is
  cut at that bound, and an exon lying entirely outside [`tx0`, `tx1`] is
  removed.
- **New skip for events outside the CDS.** If trimming leaves both slices
  identical, the splice event lies outside the coding region. The junction
  is skipped with the new fate `SKIPPED. Splice event lies outside the
  annotated coding region.` The canonical is still written with `-c`.

**The second-start-site retry had to go.** You had deferred this item.

- When the trimmed anchor was not exactly a CDS exon, the retry re-trimmed
  the junction to the gene's *second* start codon.
- With whole-exon removal this deleted real coding exons. In 20 junctions,
  exons coding in a TSL 1–2 protein-coding transcript were removed (§4.3).
- With the retry removed, that drops to zero, and no exon coding in any
  TSL ≤ 2 protein-coding transcript is ever removed.

**Output on chr15.**

- **Fates:** 103 junctions change.
  - 90 become `SKIPPED (outside coding region)`: 83 were FAILURE, 3 T4 and
    4 T1.
  - 10 move up: 6 FAILURE → T1, 2 FAILURE → T2, 2 FAILURE → T4.
  - 3 move down to FAILURE. All three are explained in §5.3.
- **FASTA:** `psq_T1.fasta` gains 5 records (164 → 169) and loses none.
  `psq_canonical.fasta` is unchanged.
- **Annotated proteins:** no sequence matching an annotated protein is lost.
  Of the 23 sequences removed from any file, none matches one. Of the 15
  sequences gained, 3 match one.

**New defect found (not fixed here).** For RI events, slice 1 contains the
upstream and downstream exons twice. This affects every RI junction (§6,
REVIEW.md HIGH-11).

---

## 2. Background: the old trimming rules

`_trim` received the gene-wide coding bounds from PR 4 (`tx0` lower, `tx1`
upper, genomic). It then applied hand-written rules per junction type, for
example for SE:

- if `tx0` is inside `anc`, cut `anc` at `tx0`;
- if `tx0` is inside `alt1`, remove `anc` and cut `alt1`;
- the mirror rules apply for `tx1` with `down` and `alt1`.

Three properties follow:

- **An exon is only touched if a bound falls inside it, or if it lies
  beyond an exon that contains one.** A UTR exon whose neighbour is also UTR
  survives. In the tx1 audit, ZNF592 MXE 433 kept all three 5′-UTR exons.
- **A5SS/A3SS rules assume `+`-strand order.** The rules always trim the
  alternative exons with `tx0` and the anchor with `tx1` (A5SS), or the
  reverse (A3SS). On `−` the exons are in the opposite genomic order, so the
  bound that actually falls in an exon is never applied.
- **RI trims only `anc_es` and `down_ee`.** The retained region (`alt1`),
  which spans both flanks, keeps its UTR ends.

---

## 3. The change

### 3.1 Trimming

```python
for exon in ['anc', 'alt1', 'alt2', 'down']:
    if absent or already removed: continue
    if ee < tx0 or es > tx1:        remove (es = ee = -1)
    else:
        if es < tx0 <= ee:          es = tx0
        if es <= tx1 < ee:          ee = tx1
```

- A bound of `None` (undeterminable, from PR 4) leaves that side untrimmed.
- The bounds are genomic and the slices are concatenated in genomic order,
  so one rule is correct for every type on either strand.
- **Boundary handling differs slightly from the old rules.** A bound on an
  exon's last base (`tx0 == ee`) now trims the exon to 1 nt; the old strict
  comparisons kept the whole exon. On chr15 this case did not occur.

### 3.2 Events outside the coding region

`Junction.event_outside_cds` compares the exons of the two slices after
trimming, as `make_slice_localgenome` composes them. If they are identical,
the event cannot change the protein. Examples:

- SE: the skipped exon was removed;
- MXE: both alternative exons were removed;
- A5SS/A3SS: the alternative part was trimmed away.

`_translate_one` returns `SKIPPED_OUTSIDE_CDS` after the read-count and
q-value filters, so existing skip reasons are unchanged.

### 3.3 Second-start-site retry removed

`trim_cds` used to retry with the second-most-upstream start codon whenever
the trimmed anchor did not match a CDS exon exactly. The retry is gene-wide
and not tied to the junction's transcript. It fired for 89 junctions under
PR 4. Why it must go once whole exons are removed is shown in §4.3.

---

## 4. Method and validation

### 4.1 Runs

[utr_audit/audit_utr.py](utr_audit/audit_utr.py) runs `main._translate_one`
on all 1,105 junctions twice, each in its own Python process:

- **Baseline:** a git worktree of `da72a0a` (PR 4).
- **Candidate:** the working tree (PR 4a).

Settings are as in `tests/test_pipeline.py`: `-r 1`, `--mask NnXx`,
`-q 0 1`, `-c`, GTF canonical, UniProt stubbed. The script records, per
junction:

- bounds, trimmed coordinates, phase and fate;
- whether the second-start retry fired.

It then compares the two versions' FASTA output and checks every
changed sequence against all proteins translated from the GTF.

### 4.2 Is anything coding removed?

Every exon the candidate removes was intersected with every `CDS` row of the
gene, from any transcript, of any biotype and any TSL.

- **Overlap with a TSL 1–2 protein-coding CDS: none.**
- **Overlap with other CDS rows: 24 junctions.** All of them involve
  transcripts that the TSL ≤ 2 filter excludes by design:

  | Best protein-coding TSL among the overlapping transcripts | Junctions |
  |---|---|
  | TSL 5 | 14 |
  | TSL 4 | 5 |
  | TSL 3 | 3 |
  | none (overlap is with NMD/other biotypes only) | 2 |

- **Parts trimmed off an exon** (rather than whole exons) lie outside
  [`tx0`, `tx1`]. Without the retry, that range is the full envelope of TSL
  ≤ 2 protein-coding start and stop codons. They can therefore only be coding
  in transcripts with no start or stop codon annotation (`cds_start_NF` /
  `cds_end_NF`) or above the TSL threshold.

### 4.3 Why the retry had to be removed

An intermediate candidate had the new trimming but kept the retry. With it:

- **49 junctions** removed an exon overlapping some CDS.
- **20 of them** removed an exon coding in a TSL 1–2 protein-coding
  transcript.
- **The retry had fired in all 20.**

Examples:

- MCTP2 SE 5653: the bound moved from the first start codon (CDS from
  94,298,266) to 94,358,548.
- MGA SE 31207: it moved past most of the gene.
- ZNF280D A3SS 1355: an exon of ZNF280D-201 (TSL 1) was removed.

Compared with that intermediate, removing the retry changes 16 fates:

- **14 junctions** are no longer wrongly skipped as outside the CDS. They
  now FAIL, as they did under PR 4.
- **FAH RI 943** goes from that wrong skip to T1.
- **ZFYVE19 SE 6502** goes from T4 to FAILURE (§5.3).

---

## 5. Results

### 5.1 Fates

| Fate | PR 4 | PR 4a |
|---|---|---|
| FAILURE | 330 | 240 |
| PARTIAL TIER 4 | 261 | 259 |
| SKIPPED (low coverage) | 297 | 297 |
| SKIPPED (outside coding region) | – | 90 |
| SUCCESS TIER 1 | 189 | 190 |
| SUCCESS TIER 2 | 26 | 27 |
| SUCCESS TIER 3 | 2 | 2 |

| Transition | Junctions |
|---|---|
| FAILURE → SKIPPED (outside coding region) | 83 |
| PARTIAL TIER 4 → SKIPPED (outside coding region) | 3 |
| SUCCESS TIER 1 → SKIPPED (outside coding region) | 4 |
| FAILURE → SUCCESS TIER 1 | 6 |
| FAILURE → SUCCESS TIER 2 | 2 |
| FAILURE → PARTIAL TIER 4 | 2 |
| SUCCESS TIER 1 → FAILURE | 1 |
| SUCCESS TIER 2 → FAILURE | 1 |
| PARTIAL TIER 4 → FAILURE | 1 |
| **Total** | **103** |

- **No exon outside the coding region remains.** 172 junctions kept one
  under PR 4; none do under PR 4a.
- **The 83 FAILURE → SKIPPED junctions** are 68 SE, 6 A5SS, 6 A3SS and 3 MXE.
  Their events lie in the 5′ or 3′ UTR of every TSL ≤ 2 protein-coding
  transcript. 9 of them involve exons coding only in TSL 3–5 transcripts.

### 5.2 Improvements

| Junction | Strand | PR 4 → PR 4a | Cause |
|---|---|---|---|
| FAM219B A3SS 1304 | − | FAIL → T1 | `−`-strand A3SS: alternative exons now cut at the stop codon |
| ADPGK A3SS 1551 | − | FAIL → T1 | same |
| GTF2A2 A3SS 2352 | − | FAIL → T1 | `−`-strand A3SS: anchor now cut at the start codon; the GTF phase becomes available |
| HEXA A3SS 5130 | − | FAIL → T1 | same |
| ANXA2 MXE 929 | − | FAIL → T1 | 5′-UTR downstream exon removed |
| FAH RI 943 | + | FAIL → T1 | 5′-UTR anchor removed and retained region trimmed (no retry) |
| CLK3 SE 5971 | + | FAIL → T2 | 3′-UTR downstream exon removed (coding only in CLK3-212, TSL 3) |
| TTBK2 A3SS 2960 | − | FAIL → T2 | `−`-strand A3SS anchor trimmed; phase becomes available |
| HDDC3 A3SS 1238 | − | FAIL → T4 | `−`-strand A3SS alternative exons cut at the stop codon |
| NMB A3SS 6113 | − | FAIL → T4 | same |

HDDC3 and NMB are two of the alternative C-terminus cases in the
[C-terminal audit](2026-11-10-alt-cterm-tiering-audit.md). PR 4c would make
them `SUCCESS TIER 1 ALT CTERM`.

### 5.3 Outcomes that get worse, each explained

| Junction | PR 4 → PR 4a | Explanation |
|---|---|---|
| ZNF280D A3SS 1355 (−) | T1 → SKIPPED | After PR 4 trimming, both alternative exons were already identical, so the PR 4 "T1" was an event-less pair: the 23-aa orphan `RKTEITCWLVWSHMMGDNPFQPK`, translated from a 5′-UTR anchor. Removing that anchor leaves identical slices, so skipping is correct. |
| TLN2 SE 34106, SE 34107, A5SS 3688 (+) | T1 → SKIPPED | The events lie outside every TSL ≤ 2 CDS. They overlap only TLN2-205/-209/-211, all TSL 5. All their PR 4 output was in `psq_T1_orphan.fasta`, and none matches an annotated protein. |
| ZNF609 SE 6146, LRRK1 SE 19275, TLN2 RI 4110 (+) | T4 → SKIPPED | Same reason: they are coding only in TSL 5 transcripts. |
| ANXA2 SE 12446 (−) | T1 → FAIL | The PR 4 T1 orphans (`ALGARPRLSC…`) were translated from a 5′-UTR exon above the start codon, which is now removed. The remaining slices have no GTF phase, and the three-frame fallback fails. Correct translation depends on PR 4b (phase from the first retained coding exon). |
| GTF2A2 SE 11945 (−) | T2 → FAIL | Same pattern: the PR 4 T2 orphans (`GWRGGRRSRN…`) were 5′-UTR-derived, and the remaining slices have no phase (PR 4b). |
| ZFYVE19 SE 6502 (+) | T4 → FAIL | The PR 4 T4 sequence began at the second start codon, which the retry chose. With the retry removed, the anchor is cut at the first start codon. That no longer matches a CDS exon exactly, so the phase lookup finds nothing (PR 4b). |

None of the sequences these junctions produced under PR 4 matches an
annotated protein.

### 5.4 Output files

| File | PR 4 | PR 4a | Lost | Gained (matching an annotated protein) |
|---|---|---|---|---|
| `psq_canonical.fasta` | 268 | 268 | 0 | 0 |
| `psq_T1.fasta` | 164 | 169 | **0** | 5 (1) |
| `psq_T1_orphan.fasta` | 46 | 37 | 11 | 2 (0) |
| `psq_T2.fasta` | 23 | 22 | 1 | 0 |
| `psq_T2_orphan.fasta` | 11 | 13 | 2 | 4 (0) |
| `psq_T3_orphan.fasta` | 4 | 4 | 0 | 0 |
| `psq_T4.fasta` | 216 | 217 | 3 | 4 (2) |
| `psq_T4_orphan.fasta` | 81 | 75 | 6 | 0 |

None of the 23 lost sequences matches an annotated protein. Every changed
sequence is listed in `fasta_diff.tsv`.

### 5.5 Phase availability

Removing exons changes which exon `get_translated_phase` inspects:

- 13 junctions lose their GTF phase;
- 5 gain one.

This is the issue PR 4b addresses: take the phase from the first retained
coding exon and enforce one-frame translation. ANXA2 SE 12446, GTF2A2 SE 11945
and ZFYVE19 SE 6502 should be re-checked there.

---

## 6. New finding: RI slice 1 duplicates its flanking exons (REVIEW.md HIGH-11)

**The duplication.**

- In the rMATS RI file, `riExonStart_0base..riExonEnd` spans the whole
  retained exon, from the upstream exon's start to the downstream exon's end.
  This is so in all 88 chr15 rows.
- JCAST maps that span to `alt1`.
- `make_slice_localgenome` then builds slice 1 as `anc + alt1 + down`.

The upstream and downstream exons therefore appear twice. For example, in
UNC45A RI 916, slice 1 is 37 + 515 + 324 = 876 nt, where the true retained
transcript is 515 nt.

**Consequences.**

- A duplication whose length isn't a multiple of 3 usually shifts the frame
  into a stop codon. That fits the RI fates under PR 4a: 57 of 88 are T4 and
  only 1 is T1.
- The forced T4 translation of slice 1 begins with the anchor exon, then
  runs into the duplicated sequence.

**Fix.** For RI, slice 1 should be `alt1` alone (after trimming).

This is outside PR 4a's scope and changes output, so it gets its own PR and
audit. It is proposed as **PR 4a2**, before the phase work in PR 4b: phase
lookup for RI depends on the slice composition.

---

## 7. Follow-ups

1. **PR 4a2: RI slice composition** (§6).
2. **PR 4b: phase from the first retained coding exon** (§5.5). Re-check the
   three cases in §5.3.
3. **One gene-wide envelope**
   ([tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions)
   item 2, deferred). Transcripts with a later start codon (ANXA2 SE 12446)
   need per-transcript bounds. The retry was a crude substitute.
4. **PR 4c: alternative C-termini.** Re-run
   [that audit](2026-11-10-alt-cterm-tiering-audit.md) after 4a2 and 4b,
   since it was measured on PR 4 alone.
5. **The TSL threshold now matters more.** Exons coding only in TSL 3–5
   transcripts are removed (§4.2), so events confined to them are skipped.
   Users who want those isoforms can raise `params.tsl_threshold`. It is
   not a CLI option yet.

---

## 8. Reproducing this audit

From the repository root, with the package dependencies installed:

```bash
python reports/utr_audit/audit_utr.py               # baseline da72a0a vs the working tree
python reports/utr_audit/audit_utr.py --baseline X  # any other baseline commit
```

The script:

- creates and removes its own git worktree for the baseline;
- copies the GTF to a temporary directory, so no `.cached` file is left in
  the repository;
- stubs out UniProt;
- rewrites the three evidence files in `reports/utr_audit/`.

The inputs and their checksums are the same as for the
[tx1 audit, §9](2026-11-03-tx1-translation-bounds-audit.md#9-reproducing-this-audit).

After PR 4a is committed, `--baseline da72a0a` still compares against PR 4.
The candidate is whatever the working tree holds, so to reproduce the
evidence files here exactly, run the script at commit `89e1a82` (PR 4a).
Later PRs, starting with PR 4a2's RI change, alter the candidate.
