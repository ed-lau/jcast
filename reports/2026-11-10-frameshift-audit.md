# Audit: the T1/T2 frameshift test and per-slice frames (PR 4b2)

| | |
|---|---|
| **Date** | 2026-11-10 |
| **Status** | Implemented. A regression from PR 4b, found in the [alternative-termini audit §6](2026-11-10-alt-termini-audit.md#6-regression-found-the-frameshift-test-ignores-per-slice-frames-proposed-pr-4b2). Baseline is commit `5c0a697`. |
| **Code** | `Sequence.frameshift` ([jcast/sequences.py](../jcast/sequences.py)) |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | `python reports/frameshift_audit/audit_frameshift.py` at the PR 4b2 commit (≈ 2 min, no network) |
| **Evidence** | [frameshift_audit/summary.txt](frameshift_audit/summary.txt) · [frameshift_audit/junctions.tsv](frameshift_audit/junctions.tsv) · [frameshift_audit/fasta_diff.tsv](frameshift_audit/fasta_diff.tsv) |

## 1. Defect and fix

T1 and T2 differ only in whether the splice event shifts the reading frame
of the sequence downstream of it. `Sequence.frameshift` tested this as:

```python
(len(slice1_nt) - len(slice2_nt)) % 3 != 0
```

That is correct only if both slices are read from the same starting frame.
Since PR 4b, each slice has its own frame (`Junction.slice_phases`). That
happens when trimming leaves the slices starting in different exons. For
such slices, what matters is the *coding* length, the slice length minus
its frame offset:

```python
((len(slice1_nt) - phase1) - (len(slice2_nt) - phase2)) % 3 != 0
```

The PR 4b audit validated frames and products but not the T1/T2 split, so
it did not catch this.

## 2. Method

[frameshift_audit/audit_frameshift.py](frameshift_audit/audit_frameshift.py):

1. **Prediction.** From the pipeline's own slices and frames, it lists every
   junction that reaches the T1/T2 test (both slices translate) and whose
   slices have different frames. It evaluates both tests on each one.
2. **Comparison.** It runs every junction through `main._translate_one` at
   `5c0a697` and in the working tree, using the
   [phase audit](phase_audit/audit_phase.py) harness.
3. **Check.** Exactly the predicted junctions must change fate, in the
   predicted direction. Every changed output sequence is listed with its
   exact annotation matches.

## 3. Results

- **17 junctions** reach the T1/T2 test with slices in different frames.
  The two tests disagree on all 17.
- **Exactly those 17 change fate, each in the predicted direction:** 16 go
  from T2 to T1, and 1 (GTF2A2 SE 11950) from T1 to T2.
- **No other junction changes.**
- **No sequence is lost or created.** Records only move between the T1 and
  T2 files:

  | Junction | Before | After |
  |---|---|---|
  | LINS1 MXE 449 (exactly the annotated LINS1-207) | `psq_T2` | `psq_T1` |
  | PKM SE 36967, both slices | `psq_T2_orphan` | `psq_T1_orphan` |
  | GTF2A2 SE 11950 | `psq_T1_orphan` | `psq_T2`, joined at the N-terminus as T2 allows |

| Junction | Frames | Slice nt | Before | After |
|---|---|---|---|---|
| ADAL MXE 2115 | 1, 0 | 61, 192 | T2 | T1 |
| LINS1 MXE 449 | 2, 0 | 383, 132 | T2 | T1 |
| GTF2A2 SE 11945 | 0, 2 | 153, 26 | T2 | T1 |
| GTF2A2 SE 11950 | 0, 2 | 98, 26 | **T1** | **T2** |
| HYKK SE 17536 | 0, 2 | 477, 140 | T2 | T1 |
| WDR61 SE 22528 | 2, 0 | 243, 94 | T2 | T1 |
| OAZ2 SE 24775 | 0, 2 | 183, 131 | T2 | T1 |
| MGA SE 31207 | 0, 1 | 2013, 949 | T2 | T1 |
| PKM SE 36967 | 0, 2 | 201, 167 | T2 | T1 |
| DET1 SE 5744 | 0, 1 | 1116, 1093 | T2 | T1 |
| CD276 SE 5993 | 0, 2 | 155, 76 | T2 | T1 |
| DPP8 SE 6124 | 0, 2 | 307, 270 | T2 | T1 |
| ZNF280D SE 6275 | 0, 2 | 175, 147 | T2 | T1 |
| ZNF280D SE 6276 | 0, 2 | 175, 147 | T2 | T1 |
| TMEM62 SE 6443 | 0, 2 | 250, 138 | T2 | T1 |
| UBE3A SE 6665 | 0, 1 | 62, 42 | T2 | T1 |
| UBE3A SE 6666 | 0, 1 | 62, 42 | T2 | T1 |

**Fates after PR 4b2:**

| Fate | Junctions |
|---|---|
| FAILURE | 157 |
| PARTIAL TIER 4 | 302 |
| SKIPPED (low coverage) | 297 |
| SKIPPED (outside coding region) | 92 |
| SUCCESS TIER 1 | 227 |
| SUCCESS TIER 2 | 30 |

`psq_T1.fasta` grows from 170 to 171 records.

**Erratum in the termini audit.** That audit's §6 listed PKM SE 36967 as
T1 → T2. Its query used a column that, for N-anchored slices, is measured
from the alternative start codon. The erratum note in that report records
the correction.

## 4. Regression tests

| Test | Asserts |
|---|---|
| `test_per_slice_phase` | GTF2A2 SE 11945, frames (0, 2): T1 (it was T2) |
| `test_per_slice_phase_reconstructs_annotated_isoform` | LINS1 MXE 449: T1, and `psq_T1.fasta` contains the LINS1-207 protein (it was in `psq_T2`) |
| `test_frameshift_counts_each_slice_frame` | GTF2A2 SE 11950, frames (0, 2): T2 (it was T1) |

All three fail at `5c0a697`. The first two had been written in PR 4b and
asserted the defective behaviour; they were corrected here.
