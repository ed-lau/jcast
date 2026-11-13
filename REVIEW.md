# JCAST — Critical Review & v0.4.0 Hardening Roadmap

*Review date: 2026-11-03. Reviewed version: v0.3.5a (HEAD `89a9dea` plus the
uncommitted `jcast/gui/` skeleton). Supersedes the 2026-05-23 review.*

> **Status, 2026-11-13.** The sections below record the state at the review
> date, with each item's fix marked; §3 describes the current code. Progress
> since then is tracked in the roadmap
> ([§9](#9-v040-hardening-roadmap-revised)): PRs 0–10 and 14–17 are done.
> - JCAST installs and runs on current Python, pandas and numpy.
> - The tiering issues of §5 are fixed, each with an audit under
>   [reports/](reports/README.md).
> - With `--uniprot` and `-p 16`, a whole-genome run takes about 4 minutes
>   with `-m`, and about 15 minutes translating every junction. It used to
>   take days.
> - Each run records every junction's outcome in `junctions.tsv`, and the
>   log holds the run's settings and a summary.
>
> Open: PRs 11–13 (CI, argparse hygiene, GUI).

This document is the developer-facing companion to [README.md](README.md) and
[CHANGELOG.md](CHANGELOG.md). It catalogs the current state of the codebase,
correctness and scientific-integrity issues, and the hardening roadmap to
v0.4.0. The audience is JCAST maintainers and new contributors.

**How this revision was verified.** Unlike the first review, every claim below
was checked against running code, not just by reading it:

- A fresh environment resolved from today's indexes (Python 3.11, pandas 3.1,
  numpy 2.4, matplotlib 3.11, scikit-learn 1.9, gtfparse 1.3) plus the
  maintainer's local `.conda` env (pomegranate 1.1.2, torch 2.6).
- The test suite, run as-is.
- A full end-to-end tiering run over the bundled chr15 data (1,105 junctions,
  GTF-only canonical so it needs no network). It was run twice: as shipped,
  and with the `tx1` correction from §5.1 so the two outputs could be diffed.
- Cutoff selection compared across four mixture-model implementations on the
  13 bundled `tests/data/model/` datasets.

## Contents

0. [Summary](#0-summary)
1. [Corrections to the 2026-05-23 review](#1-corrections-to-the-2026-05-23-review)
2. [The pomegranate break and the immediate fix](#2-the-pomegranate-break-and-the-immediate-fix)
3. [Architecture overview](#3-architecture-overview)
4. [Bugs & correctness issues](#4-bugs--correctness-issues)
5. [Scientific integrity assessment](#5-scientific-integrity-assessment)
6. [Performance](#6-performance)
7. [Dependencies, packaging & environment](#7-dependencies-packaging--environment)
8. [Test coverage](#8-test-coverage)
9. [v0.4.0 hardening roadmap (revised)](#9-v040-hardening-roadmap-revised)
10. [GUI plan — notes](#10-gui-plan--notes)
11. [Appendix A — proposed `model.py` replacement](#appendix-a--proposed-modelpy-replacement)

---

## 0. Summary

**JCAST does not run on any freshly installed environment today.** Three
independent breakages stack up:

1. `jcast/model.py` imports `pomegranate` at module top level, and
   `jcast/main.py` imports `model` unconditionally. Any pomegranate ≥ 1.0
   therefore breaks *every* CLI invocation, not just `-m`
   ([§2](#2-the-pomegranate-break-and-the-immediate-fix)). pomegranate 0.x
   no longer builds on Python 3.11.
2. A gzipped genome (`.fa.gz`, which the README advertises) crashes on every
   run with `UnicodeDecodeError` (HIGH-8, upgraded).
3. `requirements.txt` has no upper bounds. It now resolves pomegranate 1.x,
   gtfparse 3.x (polars-based, incompatible), and pandas 3.x.

Once those are fixed, the core tiering logic runs cleanly on pandas 3 /
numpy 2. The full chr15 run raised no exceptions.

**The most important scientific finding** is that the `tx1` (translation-end)
bug is **not silent**, contrary to the first review. Splice events at the
C-terminal exon are essentially untranslatable (2 of 44 on chr15), and 12 %
of junctions get spurious trims at start codons. Correcting it changes the
fate of 48 of 1,105 junctions and recovers 12 additional T1 sequences
(none lost). Fixed in PR 4
([§5.1](#51--tx1-is-derived-from-start_codon-not-stop_codon--measured)).

---

## 1. Corrections to the 2026-05-23 review

| Item | First review said | Verified status |
|---|---|---|
| HIGH-6 `np.int_` | Removed in numpy 2 | **Wrong.** `np.int_` still exists in numpy 2.x (the removed aliases were `np.int`, `np.float`, …). Withdrawn. |
| HIGH-8 `FileType('r')` genome | "Undefined behavior" on `.gz` | **Worse than stated.** It is a deterministic `UnicodeDecodeError` on every `.fa.gz` run, because `gzip.open()` receives a text-mode handle. |
| HIGH-1 `except TypeError or ValueError` | HIGH | Construct is wrong, but the comparisons it guards can only raise `TypeError`, so there is no observed impact. **Reclassified LOW** (fix in passing). |
| HIGH-7 test suite | `test_model.py` "tests `gaussian_mixture` only" | `test_model.py` is *also* broken: it calls `RmatsResults` without `logger`, and it can't even be imported without pomegranate. Adding `logger=` to the other tests is not sufficient either, because they pass `str` paths while the code calls `path.name` (see HIGH-7). |
| §4.1 `tx1` | "Most of the time the bug is silent" | **Wrong.** The end-trim fires on 127 of 1,105 chr15 junctions; see §5.1 for measured output differences. |
| §2 `__init__.py` | "`from jcast import *` yields the empty set" | **Wrong.** Names in `__all__` cause the submodules to be imported, so it imports every module, which today triggers the pomegranate error. Still worth removing. |
| §2 rMATS reader | "explicit column lists … not positional assumptions" | **Wrong.** `df.columns = [...]` renames *by position* and ignores the rMATS header. It works for current rMATS-turbo output, but a format change would silently mislabel columns. |
| §5 pomegranate options | Recommended scikit-learn `GaussianMixture` | Revised after measurement. A dependency-free scipy EM reproduces the Gamma+Normal model; the sklearn model produces systematically different cutoffs. See §2. |
| §5 gtfparse | "2.x exists, migration separate" | gtfparse is now at **3.0.5**. `requirements.txt` (`>=1.2.1`) resolves it today, and `read_gtf` returns a polars frame, so `.to_csv` fails. Needs an immediate `<2` bound. |
| §1 Quarto site | "The Quarto site under `jcast_website/`" | `jcast_website/` is in `.gitignore`. The docs *source* is not version-controlled, only the rendered `docs/`. |

Everything else in the first review was confirmed, and is carried forward below
with its original ID so existing references (README → HIGH-2/HIGH-3) still
resolve.

---

## 2. The pomegranate break and the immediate fix

### Symptom

```
ImportError: cannot import name 'GeneralMixtureModel' from 'pomegranate'
```

pomegranate 1.0 (2023) was a ground-up rewrite on PyTorch. The top-level
`GeneralMixtureModel`, `GammaDistribution` and `NormalDistribution`, and the
`from_samples()` constructor, are gone. `setup.py` pins `<1.0`, but
`requirements.txt` (`>=0.14.5`) does not. Note also that the import sits at
the top of [jcast/model.py:13](jcast/model.py#L13), and
[jcast/main.py:10](jcast/main.py#L10) imports `model` unconditionally.
**The CLI fails even when `-m` is not used.**

### Options evaluated

Each option was run on the 12 non-empty `tests/data/model/` datasets. Cutoffs
are `min_count` as currently defined (count + 1 scale; see MED-11).

| Dataset | (a) scipy EM, Gamma+Normal | (b) pomegranate 1.1.2 port | (c) sklearn Box-Cox GMM (existing `gaussian_mixture`) |
|---|---|---|---|
| encode_human_adrenalgland | 13 | 12 | 6 |
| encode_human_heart | 15 | **error** (variances must be positive) | 6 |
| encode_human_lung | 8 – 11¹ | 12 | 7 |
| encode_human_spleen | 6 | **error** | 5 |
| encode_human_stomach | degenerate → fallback | 5 | 9 |
| encode_human_testis | 9 | **error** | 8 |
| encode_human_uterus | 6 | 6 | 9 |
| encode_mouse_heart | 7 | 3 (collapsed, sd = 0.00) | 7 |
| heart_mouse_old | 30 | 32 | 17 |
| ipsc_human_d0 | 15 | 16 | 9 |
| ipsc_human_d14 | 16 | **error** | 10 |
| kidney_old_totalRNA_gencodeM26 | 5 | **error** | 1 |

¹ Sensitive to the EM convergence tolerance (11 at `tol=1e-8`, 8 at `1e-6`).
See §5.6.

The original pomegranate 0.14.9 baseline could not be rebuilt on Python 3.11
(build failure), so it can't serve as a reference going forward either.

- **(a) Hand-rolled EM (recommended).** Same model family as v0.3.5
  (Gamma low component + Normal high component on ln(count)). It needs only
  numpy/scipy, which are already dependencies. It is deterministic, so the
  50-restart + KDE-of-parameters machinery is no longer needed, and it runs in
  about 1 s per dataset. Where pomegranate 1.x converges, the two agree within
  ±1–2 reads. A degeneracy guard returns `None`, so the caller falls back to
  `-r` with a warning instead of crashing. Code is in
  [Appendix A](#appendix-a--proposed-modelpy-replacement), about 70 lines and
  tested on all 12 datasets with matplotlib 3.11.
- **(b) Port to pomegranate 1.x.** **Not viable.** It fails on 5 of 12
  datasets, degenerates on a sixth, and pulls in PyTorch (~1 GB).
- **(c) Switch to the existing sklearn path.** It works on all 12, but it is
  a different model (Box-Cox Gaussian). Its cutoffs are systematically lower
  (heart 6 vs 15, kidney 1), which would silently change results relative to
  v0.3.5. If adopted, it needs its own changelog entry.
- **(d) Pin `pomegranate<1`.** Not installable on current Python. Dead end.

### Immediate patch (PR 0) — ✅ applied 2026-11-03

Verified with no pomegranate installed: the model reproduces the table above
(stomach falls back), the sklearn path still runs, and `python -m jcast … -m`
completes end to end on chr15 (pandas 3.1, numpy 2.4, matplotlib 3.11). In
addition to the list below, PR 0 replaced the removed
`pd.core.computation.ops.UndefinedVariableError` with
`pandas.errors.UndefinedVariableError` (first half of MED-7), so the pandas
bound could be widened to `>=1.5,<4` instead of forcing `numpy<2`.

1. Replace `general_mixture_model` / `plot_general_mixture_model` with
   Appendix A. Drop the `pomegranate` import, and move the `sklearn` import
   inside `gaussian_mixture` so it is only needed when that path is used.
2. In `run_jcast`, handle `min_count is None`: log a warning and use
   `args.read`.
3. Fix `get_junction_count_array` to use
   `pd.concat(self.rmats_list, ignore_index=True)` (HIGH-4). Without this,
   `-m` crashes on pandas ≥ 2 right after the import is fixed.
4. Remove `plt.style.use('seaborn-white')` (HIGH-5). It raises `OSError` on
   matplotlib ≥ 3.6, and it has no effect after plotting anyway.
5. Remove `pomegranate` from `setup.py` / `requirements.txt`.

---

## 3. Architecture overview

JCAST is a small flat package (~3,400 lines of Python, 14 modules plus the GUI
skeleton) that runs as a CLI, in one process or, with `-p`, several.

```
jcast/
├── __main__.py        Entry — dispatches to gui.launch_gui (--gui) or main.main()
├── main.py            CLI parsing, run loop (one process or -p workers), tier dispatch (_translate_one)
├── junctions.py       RmatsResults reader; Junction (trim_cds, get_translated_phase, min_jc)
├── sequences.py       Sequence: nucleotide retrieval, translation, canonical lookup
│                       (UniProt FASTA, GTF or UniProt API), stitching, FASTA writing
├── uniprot.py         --uniprot: UniProtFasta (canonical by sequence, then unique gene name),
│                       uniprot_canonicals.tsv, psq_isoforms.tsv, jcast_uniprot_merged.fasta
├── annots.py          ReadAnnotations (GTF, indexed by gene), ReadGenome (FASTA), AnnotatedTranscript
├── helpers.py         get_local_nuc, get_complementary, make_pep, write_seqrecord_to_fasta
├── model.py           Mixture-model read-count cutoff (numpy/scipy EM; vestigial sklearn path)
├── outcomes.py        junctions.tsv and the end-of-run summary
├── cache.py           The user cache directory (platformdirs, or --cache_dir)
├── params.py          Module-level run parameters (TSL threshold, anchor length, …)
├── constants.py       Genetic code table
├── fates.py           Per-junction outcome strings and their short names
├── logger.py          Logger setup
└── gui/               Optional PySide6 skeleton
```

### Data flow

1. **`RmatsResults`** reads all five `*.MATS.JC.txt` files (MED-8) and
   renames their columns positionally.
2. **`ReadAnnotations`** parses the GTF via `gtfparse`, caches the table as
   a pickle in the user cache directory (PR 7), and indexes its rows by
   gene (PR 9).
3. **`ReadGenome`** loads the genome FASTA (optionally gzip) into a `SeqIO`
   dict held in memory.
4. **`UniProtFasta`**, with `--uniprot`, reads the supplied UniProt FASTA
   (PR 14).
5. **Optional `-m`**: fit a mixture to ln(count + 1) and use the decision
   boundary as the read cutoff, in reads (PR 6).
6. **Per junction**, in rMATS file order. With `-p`, workers translate and
   the main process writes their results in that order (PR 15).
   1. The read-count and FDR filters (PR 6). A skipped junction only
      writes its gene's canonical, with `-c`.
   2. `trim_cds` trims every exon to the gene's coding region and removes
      exons outside it (PR 4/4a).
   3. `get_translated_phase` reads each slice's frame at its 5′ end (PR 4b).
   4. Slice nucleotides are built (RI per PR 4a2). A slice with a base
      outside `ACGT` is skipped (PR 3).
   5. Canonical lookup: the UniProt FASTA, then the GTF (with the
      `Ensembl_canonical` fallback, PR 17); or the UniProt API without
      `--uniprot`.
   6. The junction is skipped if its event lies outside the coding region.
   7. Tier T1, T1 ALT TERMINUS, T2 or T4, or fail. The records are written
      to `psq_*.fasta`, and a row to `junctions.tsv` (PR 16).

   An exception fails only its junction (PR 7).
7. **At the end:** the `--uniprot` outputs, then the log's summary of fates
   and records per file (PR 16).

### Tier classification (current, after PR 4c)

Before tiering, each junction is processed in four steps:

1. **Trimming.** The junction is trimmed to the gene's annotated coding
   region, [`tx0`, `tx1`]. This runs from the most upstream start codon to
   the most distal stop codon of protein-coding transcripts at TSL ≤
   `params.tsl_threshold`. Exons entirely outside it are removed.
2. **Slices.** The two slices are built from the retained exons. For RI,
   slice 1 is the retained exon alone.
3. **Frames.** Each slice gets its own reading frame. It is read at the
   slice's 5′ end from the protein-coding CDS rows covering it: first at the
   TSL threshold, then any TSL, taking the majority frame if transcripts
   disagree.
4. **Translation.** Each slice is translated once, in that frame. A slice
   with no annotated frame is not translated, and no frame is ever guessed.

| Tier / fate | Meaning |
|---|---|
| **T1** | Both slices translate in their annotated frames with no stop codon before the last codon. Their coding lengths (slice length minus that slice's frame; PR 4b2) differ by a multiple of 3. |
| **T1 ALT TERMINUS** | New in PR 4c. At least one slice is anchored at an annotated terminus (TSL ≤ 2 protein-coding), and both reach an acceptable end. *C-anchored*: its first stop is an annotated stop codon; it is joined at its N-terminus only. *N-anchored*: it lacks the gene's start codon that the other slice has, and is translated from another annotated start codon in it; it is joined at its C-terminus only. The output file depends on which terminus differs from the canonical: `psq_T1_alt_cterm`, `psq_T1_alt_nterm`, `psq_T1_alt_both`, or `psq_T1` if neither (`_orphan` if unjoinable). The frameshift test is waived only when a slice ends at an internal annotated stop. |
| **T2** | Same as T1, but the coding lengths differ by a non-multiple of 3 (frameshift). An unstitchable slice may be joined at its N-terminus only. |
| **T4** | Exactly one slice translates. The other is translated in its frame **up to** its first stop codon (a truncated product; nothing if it has no frame). The product is joined at its N-terminus only if it is ≥ `params.ptc_threshold` × the translated slice's length. Otherwise, or if it cannot be joined, it goes to `psq_T4_orphan.fasta`. |
| ~~T3~~ | **Removed in PR 4b.** It re-translated in a guessed frame when the annotated frame failed. |
| Skip | One of four reasons: low count; FDR outside `-q`; masked nucleotides; or, new in PR 4a, the event lies outside the annotated coding region (both slices are identical after trimming). The canonical is written if `-c` is set. |
| Fail | Neither slice translates: no annotated frame, or a stop codon in both. |

**Output files.**

- T1/T2/T4 slices that stitch to the canonical go to `psq_T#.fasta`, and
  T1 ALT TERMINUS slices to `psq_T1_alt_*.fasta`.
- Slices that do not stitch go to the matching `_orphan` file.
- A stitched slice identical to the canonical is written to
  `psq_canonical.fasta`. With `-c`, every gene's canonical is (PR 8).
- `junctions.tsv` has one row per junction; `logfile.log` has the run's
  settings, warnings and summary (PR 16).
- With `--uniprot`: `uniprot_canonicals.tsv`, `psq_isoforms.tsv` and
  `jcast_uniprot_merged.fasta` (PR 14).

**Known limitations of the current tiers**, with the audit for each:

- **Alternative termini need annotation.** An alternative stop or start
  codon that the annotation lacks, or has only above the TSL threshold, is
  still read as a premature stop (T4/FAIL) or a mid-CDS fragment.
  [Termini audit](reports/2026-11-10-alt-termini-audit.md).
- **Start-less slices without an annotated alternative start** (34 on chr15)
  are still read from mid-CDS. They inflate T1/T2 counts without adding
  wrong records.
  [Phase audit §5.4](reports/2026-11-10-phase-audit.md).
- **One coding region per gene.** Transcripts with a later start codon or
  an earlier stop are trimmed to the gene-wide range (ZFYVE19 SE 6502).

---

## 4. Bugs & correctness issues

Tags: **HIGH** blocks running on current environments or biases scientific
output. **MED** reduces robustness or gives wrong output in edge cases.
**LOW** is UX or style. ✔ = reproduced by execution; ◎ = confirmed by
reading.

### HIGH-2 — `make_pep` hard-fails on non-ACGT bases ✔ — **fixed in PR 3**

[jcast/helpers.py:226](jcast/helpers.py#L226). `code[nt[i:i+3]]` raises
`KeyError` for lowercase (soft-masked) and IUPAC-ambiguous bases. Reproduced
for `atg…`, `NNN`, and `RAA` on both strands. Fix: `nt = str(nt).upper()`
and `code.get(codon)`, with an explicit policy for unknown codons: either
emit `X`, or treat it like a stop and skip.

### HIGH-3 — `get_complementary` raises before its fallback ✔ — **fixed in PR 3**

[jcast/helpers.py:148-152](jcast/helpers.py#L148-L152). `raise` makes the
`comp += 'X'` line unreachable. Fix together with HIGH-2, ideally using
`Bio.Seq.reverse_complement`, which already handles IUPAC codes and case.

**PR 3** ([audit](reports/2026-11-11-masked-bases-audit.md)):
- `get_local_nuc` uppercases genome sequence when it is read.
- A slice with any base outside `ACGT` is skipped (`SKIPPED_MASKED`).
- `make_pep` returns `''` for a codon it cannot translate.
- `get_complementary` is `Bio.Seq.reverse_complement`.
- `--mask` is deprecated and ignored.
- On chr15, a soft-masked genome now reproduces the unmasked output exactly.

### HIGH-4 — `DataFrame.append` removed in pandas 2 ✔ — **fixed in PR 0**

[jcast/junctions.py:158](jcast/junctions.py#L158). Only the `-m` path hits
this. Fix: `pd.concat(self.rmats_list, ignore_index=True)`.

### HIGH-5 — `seaborn-white` style removed in matplotlib 3.6 ✔ — **fixed in PR 0**

[jcast/model.py:101](jcast/model.py#L101),
[jcast/model.py:196](jcast/model.py#L196). Raises `OSError`. Delete the
calls; they run after plotting, so they never had an effect.

### HIGH-7 — Test suite is broken ✔ — **fixed in PR 1**

6 of 7 tests fail, and `test_model.py` fails at import. Two separate problems
need fixing:

- Every constructor call omits the required `logger` argument.
- The tests pass `str` paths, but `ReadAnnotations`/`ReadGenome` call
  `self.path.name`. The classes are coupled to `argparse.FileType`.

Fix HIGH-8 first so the classes accept plain paths. Then add `logger=` to the
tests.

### HIGH-8 — `argparse.FileType('r')` for GTF and genome; `.fa.gz` always crashes ✔ *(upgraded)* — **fixed in PR 1**

[jcast/main.py:362-369](jcast/main.py#L362-L369),
[jcast/annots.py:88-91](jcast/annots.py#L88-L91). `ReadGenome` calls
`gzip.open(self.path, 'rt')` on the text-mode handle argparse opened, and the
gzip header read fails with `UnicodeDecodeError`. **Every run with a
gzipped genome crashes**, even though the README says "optionally .gz". Fix:
`type=str` plus an `os.path.isfile` check, and make both classes take a path.
Until then, add a Known-issues entry to the README.

### HIGH-9 — Top-level optional imports break the whole CLI ✔ *(new)* — **fixed in PR 0**

[jcast/model.py:11-13](jcast/model.py#L11-L13) imports `pomegranate` and
`sklearn` at module scope, and [jcast/main.py:10](jcast/main.py#L10) imports
`model` unconditionally. A missing or incompatible package for an optional
feature (`-m`) therefore kills every invocation. `scikit-learn` is not in
`install_requires`, so a plain `pip install jcast` hits this too. Fix: the
§2 patch, plus `scikit-learn` in `install_requires` or lazy imports.

### HIGH-10 — `requirements.txt` resolves an incompatible stack ✔ *(new)* — **fixed in PR 0**

No upper bounds, so it currently installs pomegranate 1.1.2 (breaks the
import), gtfparse 3.0.5 (polars frame; `.to_csv` fails in `read_gtf`), and
pandas 3.1. It also disagrees with `setup.py` on every bound. See §7.

### HIGH-11 — RI slice 1 contains its flanking exons twice ✔ *(new, 2026-11-10)* — **fixed in PR 4a2**

[jcast/sequences.py:93-95](jcast/sequences.py#L93-L95). The rMATS RI
`riExonStart_0base..riExonEnd` span (mapped to `alt1`) already runs from the
upstream exon's start to the downstream exon's end. This is so in all 88
chr15 rows. Yet slice 1 is built as `anc + alt1 + down`, so both flanks
appear twice. UNC45A RI 916: 37 + 515 + 324 = 876 nt instead of 515.

Most RI junctions therefore hit a frame-shifted stop codon and land in T4
(57 of 88 under PR 4a; 1 is T1). Slice 1 should be `alt1` alone. Found
during the PR 4a audit
([report §6](reports/2026-11-10-utr-exon-trimming-audit.md)). Fixed in PR 4a2
([audit](reports/2026-11-10-ri-slice-audit.md)): 5 RI fates change, and
58 duplicated-sequence RI records are replaced, mostly in T4.

### MED-1 — `sum_sjc` returns `min(SJC, IJC)` ◎ — **fixed in PR 6**

[jcast/junctions.py:245-274](jcast/junctions.py#L245-L274). The name, the
README header field 11 ("skipped-junction read count"), and the plot axis
label all say SJC. See §5.2.

*2026-11-12, PR 6:* the minimum is kept (both isoforms need reads) and the property is renamed `min_jc`. Help text, README header field 11 and the plot now describe it ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §2.1).

### MED-2 — Three-frame fallback keeps the *last* frame; can discard a good one ◎ — **closed in PR 4b (fallback removed)**

[jcast/sequences.py:143-151](jcast/sequences.py#L143-L151). There's no
`break`. If frame 0 succeeds and a later frame fails, both slices are reset to
`''` and the junction drops to T4/fail. Add `break` after the first success,
or pick the longest translation deliberately.

*2026-11-10:* closed by removing the three-frame fallback in PR 4b. Once
frames were read correctly, the fallback changed only 15 chr15 junctions,
and no annotated sequence depended on it
([phase audit §5.3](reports/2026-11-10-phase-audit.md)).

### MED-3 — Stitch anchor is not checked for uniqueness or order ◎ — **fixed in PR 5**

[jcast/sequences.py:506-510](jcast/sequences.py#L506-L510). `str.find`
takes the first occurrence of the 10-mer, and nothing checks
`merge0 <= merge1`. Repeats can mis-stitch, and an end anchor found
*before* the start anchor duplicates canonical sequence. Require
`count == 1` for each anchor and `merge0 <= merge1`, or else send the slice to
orphan.

**PR 5** ([audit](reports/2026-11-11-stitch-anchor-audit.md)):
- The uniqueness and order checks were added.
- The audit found that most wrong joins had a different cause: a peptide
  was joined at an end that was already complete (a stop codon, or the
  gene's start codon). Those ends are no longer joined.
- The anchor length is now 7 (was 10).
- On chr15 against genomic coordinates: 0 wrong joins (was 12 of 998).

### MED-4 — Output directory never created before the log file opens ◎ — **fixed in PR 7**

[jcast/logger.py:16](jcast/logger.py#L16). Add
`os.makedirs(args.out, exist_ok=True)` before `get_logger`.

*2026-11-12, PR 7:* `run_jcast` creates it.

### MED-5 — UniProt SQLite cache written to CWD ◎ — **fixed in PR 7**

[jcast/sequences.py:397-400](jcast/sequences.py#L397-L400). Use
`platformdirs.user_cache_dir('jcast')` (machine-wide) or `args.out`. Wrap the
`INSERT` in a try/except for read-only volumes.

*2026-11-12, PR 7:* the cache is in the user cache directory (`platformdirs`), or `--cache_dir`. If it cannot be opened or written, JCAST warns and runs without it ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §3.5).

### MED-6 — GTF `.cached` written next to the input, keyed only by filename ◎ — **fixed in PR 7**

[jcast/annots.py:45-49](jcast/annots.py#L45-L49). This fails on read-only
reference volumes. It also serves a **stale cache** if the GTF is replaced
under the same name, and `to_csv` writes the index as an unnamed extra column.
Fix: cache in a user dir keyed on (path, size, mtime), with `index=False`.

*2026-11-12, PR 7:* the cache is a pickle in the cache directory, keyed by path, size, mtime and the pandas version, and renamed into place once written. The audit found a worse failure: a second run read a half-written cache without error ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §3.5).

### MED-7 — `_get_translated_region` exception handling is broken on pandas ≥ 2 ✔ *(expanded)* — **fixed in PR 4**

[jcast/junctions.py:300-338](jcast/junctions.py#L300-L338).

- `pd.core.computation.ops.UndefinedVariableError` **no longer exists**
  (verified on pandas 3.1). For a GTF without `transcript_support_level`, the
  `except` clause itself raises `AttributeError`. Use
  `pandas.errors.UndefinedVariableError` (available since pandas 1.5), or
  better, check `'transcript_support_level' in gtf.annot.columns` once.
- On the `KeyError` path, `gtf0_start` / `gtf0_end` are set to `None` and then
  passed to `len()`, which raises `TypeError`. Use empty Series instead.

### MED-8 — `RmatsResults` reads all five files regardless of `-s` ◎ *(new)*

[jcast/junctions.py:37-41](jcast/junctions.py#L37-L41) versus the existence
check at [jcast/main.py:755-757](jcast/main.py#L755-L757), which only checks
the selected types. A folder with only `SE.MATS.JC.txt` passes validation
and then raises `FileNotFoundError` (this is the existing `TODO` at
[main.py:108](jcast/main.py#L108)). Open: PR 12.
The check also uses `assert`, which is stripped under `python -O`; use
`parser.error`.

### MED-9 — `write_canonical`'s UniProt fallback is a no-op ◎ *(new)* — **fixed in PR 8**

[jcast/sequences.py:537-543](jcast/sequences.py#L537-L543). It calls
`get_canonical_aa_uniprot()` but discards the return value instead of
assigning `self.canonical_aa`. The fallback never writes anything, yet still
costs a network request, and it does so even when `params.use_gtf_only` is
set.

*2026-11-12, PR 8:* the fallback is removed. It never wrote anything: the canonical is already looked up, with `-c`, on every path ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §7 item 9).

### MED-10 — A network failure aborts the entire run ◎ *(new)* — **fixed in PR 8**

*Avoidable since PR 14: with `--uniprot`, JCAST makes no network calls
(MED-5 and MED-9 do not apply either). The API path is unchanged.*

[jcast/sequences.py:429-431](jcast/sequences.py#L429-L431). Once its retries
are exhausted, `rqs.get(...)` raises `requests.exceptions.ConnectionError` /
`RetryError`, which nothing catches. A transient EBI outage hours into a run
kills it. Catch the exception, log it, return an empty record, and let the
junction fall through to orphan/fail.

*2026-11-12, PR 8:* a failed or refused request raises `UniProtAPIError`, and that junction fails (`FAILURE (exception)`, PR 7). The run continues. The audit found a worse case than the abort: a refusal such as HTTP 429 was taken as "no canonical" and silently orphaned every record. Requests now time out after 30 s, and 429 is retried ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §4.3).

### MED-11 — Read-count thresholds are off by one, twice ◎ *(new)* — **fixed in PR 6**

- **`-r`**: the help text and README say "minimum read count for a junction
  to be translated", but [jcast/main.py:196](jcast/main.py#L196) skips when
  `sum_sjc <= args.read`. The default `-r 1` therefore requires ≥ 2 reads.
- **`-m`**: `min_count` is computed on the `count + 1` scale
  ([jcast/model.py:63](jcast/model.py#L63)) and then compared to raw counts
  with `<=`. The effective threshold is `min_count + 1` raw reads.

Neither is catastrophic, but both make reported cutoffs disagree with what
was applied. Pick one convention (`keep if count >= threshold`), document it,
and log the effective threshold.

*2026-11-12, PR 6:* a junction is kept if its count ≥ the cutoff, for `-r` and `-m`. The model returns its cutoff in reads. Measured before the fix: `-m` cut two reads above the boundary, dropping 3,410 of 131,233 whole-genome junctions. Default `-r 1` now keeps every junction ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md) §2, §7).

### MED-12 — Output is not idempotent ◎ *(new)* — **fixed in PR 7**

`psq_*.fasta` files are opened in append mode
([jcast/helpers.py:272](jcast/helpers.py#L272)), as is `logfile.log`.
Re-running into the same `-o` with different parameters silently merges
results across runs. Refuse a non-empty output dir without `--overwrite`, or
write to a temp dir and move it into place.

*2026-11-12, PR 7:* a directory holding JCAST output is refused unless `--overwrite`, which removes that output first. The log is written fresh.

### MED-14 — `-c` does not guarantee the canonical for genes whose junctions all succeed ◎ *(new, 2026-11-10)* — **fixed in PR 8**

With `-c`, the skip, FAILURE and (from PR 4c) alternative-terminus paths
write the gene's canonical. T1, T2 and T4 write it only when a slice
stitches into it exactly. On chr15, 2 protein-coding genes get no canonical
record despite `-c`. The documented intent is that the canonical is written
even when isoform translation fails. Write it on every path when `-c` is set
(PR 8, with MED-9).

*2026-11-12, PR 8:* `_translate_one` writes the canonical after every path, looking it up first where the path had not, e.g. masked junctions. Coverage is now complete: chr15 276 of 276 genes with `--uniprot`, and the whole genome (§10 of the audit).

### MED-13 — `get_logger` stacks handlers on every call ◎ *(new; matters for the GUI)* — **fixed in PR 7**

[jcast/logger.py:30-31](jcast/logger.py#L30-L31). Each `run_jcast` call adds
another `FileHandler` and `StreamHandler` to the same `'jcast'` logger. In a
long-lived process (the GUI, tests, notebooks), run *n* logs every line *n*
times and keeps writing into earlier runs' log files. Remove existing
handlers first, and close the file handler when the run ends.

*2026-11-12, PR 7:* `get_logger` closes and removes earlier handlers; `run_jcast` closes its own when it ends.

### LOW

- ~~`except TypeError or ValueError` should be `(TypeError, ValueError)`
  (was HIGH-1). Only `TypeError` can occur here, so there is no observed
  impact.~~ **Gone in PR 4a**: `_trim` was rewritten and tests `None`
  bounds explicitly.
- [jcast/helpers.py:234](jcast/helpers.py#L234): the terminal-stop tolerance
  accepts a stop followed by 0–1 trailing nt but rejects 2 trailing nt
  (`'ATGAAATAAG'` → `'MK'`, `'ATGAAATAAGG'` → `''`). Both cases are the same
  biologically.
- ~~[jcast/junctions.py:218](jcast/junctions.py#L218): `gene_symbol is None`
  never fires, because pandas parses `NA` as `NaN`. Use `pd.isna`. The README
  "defensive substitution" for `NA` gene names is therefore ineffective.~~
  **Fixed in PR 2**: `pd.isna` or blank gives `NoGeneName_<gene id>`.
  Before, headers carried `nan`, and `str()`/`repr()` of the junction raised
  `TypeError`. Nothing in the package calls them, so translation was
  unaffected.
- [jcast/helpers.py:30-35](jcast/helpers.py#L30-L35): chromosome-name
  reconciliation only strips `chr`. It never adds it (Ensembl rMATS +
  GENCODE genome → `KeyError`).
- ~~`jcast/__init__.py`: the glob-built `__all__` makes `from jcast import *`
  import every submodule. Remove it.~~ **Removed in PR 10.**
- [jcast/main.py:724-727](jcast/main.py#L724-L727): `--g_or_ln` is read by
  nothing. Remove it (from the README too). Open: PR 12.
- ~~`sequences.py`: stray `print()`.~~ **Gone.**
- [jcast/helpers.py:45-125](jcast/helpers.py#L45-L125): `get_rest_nuc` is
  unused, uses bare `except:`, and has `print`-based logging. Remove it.
- ~~`model.py`: the KDE grid `linspace(-5, 20, 500000)` silently clips any
  parameter outside that range, and `gaussian_kde` raises if all restarts
  agree.~~ **Gone in PR 0**: Appendix A replaced the model.
- `tests/data/model/encode_mouse_fetalheart/*.MATS.JC.txt` are header-only
  (0 junctions).
- Date-stamped comments (`# 2020-07-25 …`) throughout. Move them to git
  history.

---

## 5. Scientific integrity assessment

JCAST output feeds MS database searches, so silent sequence errors propagate
into identifications.

### 5.1 — `tx1` is derived from `start_codon`, not `stop_codon` — measured

> **Status: fixed in PR 4 (2026-11-10).** The implementation reproduces the audited corrected variant
> exactly: the same fate counts and the same per-file record counts (audit §5.4, §7.3).

> **Full audit:** [reports/2026-11-03-tx1-translation-bounds-audit.md](reports/2026-11-03-tx1-translation-bounds-audit.md),
> reproducible with `python reports/tx1_audit/audit_tx1.py`. The figures below
> are from that audit. They supersede the earlier prototype numbers, which
> had silently dropped the TSL filter.

[jcast/junctions.py:320-348](jcast/junctions.py#L320-L348). The block
commented "now getting end codons" queries `feature == "start_codon"`. The
C-terminal bound (`tx1` on `+`, `tx0` on `−`) therefore sits at a start
codon: 2 nt past it on `+`, and *on* it on `−`, where both bounds collapse.
It matches an annotated stop codon in 0 of 1,105 chr15 junctions.

**Two failure modes:**

- **A — C-terminal events are untranslatable.** 63 junctions span the
  stop-codon exon. The shipped code translates 2 of the 44 that pass the read
  filter; corrected bounds translate 36 (7 T1, 4 T2, 25 T4). The slice
  otherwise runs into the 3′ UTR and the true stop is read as a PTC.
- **B — spurious trims at start codons** (129 junctions). These cut the
  start-codon exon a few nt after the ATG, or on `−` keep its 5′-UTR side.
  Some T1/T2 "successes" are 5′-UTR translations, or slice pairs from which
  the splice event was deleted.

**Output impact** (corrected bounds as specified in the audit, §7):

- 48 of 1,105 junctions change fate, all in classes A/B; the other 917 are
  identical. 38 improve. The 10 that get worse were all traced to artefactual
  shipped "successes".
- `psq_T1.fasta`: 152 → 164. **12 gained (8 match an annotated Ensembl
  protein), 0 lost.** `psq_canonical.fasta` is unchanged. One annotated-match
  sequence is no longer produced: a 43-aa BLOC1S6 orphan fragment that arose
  from an artefactual trim.

**Recommendations:** implement the bound definition in audit §7.1 (no `+2`
offsets; each bound falls back to `None` independently). Land the regression
tests in audit §7.2 first, then update the golden test (audit §7.3) and the
CHANGELOG (audit §7.4) in the same commit. Residual issues that exist under
both versions are listed in audit §8.

**Decisions on the residual issues (2026-11-10):**

| Audit §8 item | Decision |
|---|---|
| 1. Pure-UTR exons are never removed | Drop them. Translating UTR is outside JCAST's purpose. **Done in PR 4a.** |
| 2. One gene-wide CDS envelope | Defer. Revisit after PR 4. |
| 3. Second-start-site retry | Deferred, but **removed in PR 4a**: once UTR exons are removed, it deleted exons coding in TSL 1–2 transcripts in 20 junctions ([PR 4a audit §4.3](reports/2026-11-10-utr-exon-trimming-audit.md)). |
| 4. Phase lookup after trimming | Each slice's frame is read at its 5′ end, and translation uses one frame only. **Done in PR 4b.** |
| 5. Three-frame fallback (MED-2) | **Removed in PR 4b**, together with Tier 3 and the `psq_T3*.fasta` outputs. The audit showed only 15 junctions depended on it: 3 T3 guesses, 10 FAILUREs that hid a translatable slice, and 2 pseudogene T1s. No annotated sequence was lost ([phase audit §5.3, §7](reports/2026-11-10-phase-audit.md)). |
| *(PR 4b finding)* Slices lacking the start codon | **Keep reading a frame for them.** Withholding it removed real isoforms (GTF2A2-208, LINS1-207). An annotated-start rule, mirroring the annotated-stop rule, is to be added in PR 4c ([phase audit §5.4](reports/2026-11-10-phase-audit.md)). |

Alternative C-terminal isoforms, which today are classed as PTC (T4/FAIL), are
evaluated in
[reports/2026-11-10-alt-cterm-tiering-audit.md](reports/2026-11-10-alt-cterm-tiering-audit.md).
That work, extended to alternative N-termini, is PR 4c. All audits and
decisions are indexed in [reports/README.md](reports/README.md).

### 5.2 — `sum_sjc` is `min(SJC, IJC)` ◎ — **resolved in PR 6**

The cutoff model, the FASTA header (`r<count>`), the plot axis, and the README
all describe this as the skipped-junction count. The v0.3.2 changelog says
"sum of all replicates", and the 2022-05 switch to `min` is undocumented.
Decide on one semantics, rename the property (`min_jc`), and update the
header documentation.

*2026-11-12, PR 6:* see MED-1.

### 5.3 — Soft-masked / ambiguous bases ✔ — **fixed in PR 3**

Covered by HIGH-2/3. The README's workaround also suggested hard-masked
genomes, which silently skip every junction touching a repeat (19% on
chr15). The PR 3 audit found this, and the README now advises unmasked or
soft-masked genomes.

### 5.4 — Stitch anchors can mis-place sequence ◎ — **fixed in PR 5**

See MED-3.

### 5.5 — Tier 4 caveats ◎

The README now surfaces the warning prominently. No further change needed.

### 5.6 — The modeled read cutoff is a weakly identified heuristic ✔ *(new)*

The input is a discrete lattice (ln 2, ln 3, …) dominated by low counts. The
"expressed" Normal component is usually a broad shoulder, not a separate mode.
Consequences, observed across the 12 bundled datasets:

- One dataset (`encode_human_stomach`) has no stable two-component fit under
  any implementation. EM slowly absorbs the Gamma component into the ln 2
  spike.
- The cutoff for `encode_human_lung` moves from 11 to 8 when only the EM
  convergence tolerance changes.
- pomegranate 1.x and the sklearn model disagree with each other by up to 2×.

**Recommendation.** Keep `-m`, but present it as a heuristic:

- Log the fitted parameters and effective cutoff.
- Fall back to `-r` with a warning when the fit is degenerate.
- Write the cutoff into a run manifest in the output dir.
- Don't describe it as a statistical decision boundary in the docs.

**Status:** the parameters and the cutoff (in reads, PR 6) are logged, and a
degenerate fit falls back to `-r` with a warning (PR 0). Since PR 16 the
log is the run's record: version, arguments, inputs, model and cutoff.

---

## 6. Performance *(new section)*

None of this matters for the chr15 fixture, but all of it matters genome-wide
(~10⁵ junctions, ~3 M GTF rows).

> **Status, 2026-11-12:** the first three are fixed: the per-gene index and
> the in-memory de-duplication in PR 9, filters first in PR 6. Together
> with `-p` (PR 15) a whole-genome run went from days to minutes
> ([parallel audit](reports/2026-11-12-parallel-audit.md)). The fourth
> applies only to the API path; `--uniprot` makes no API calls.

- **Repeated full-GTF scans.** `gtf.annot.query('gene_id == …')` scans the
  entire GTF 3–5 times per junction (`_get_translated_region` ×1–2,
  `get_translated_phase`, `trim_cds`, `get_annotated_transcripts_gtf`). Build
  `dict(tuple(annot.groupby('gene_id')))` once, or index by `gene_id`. This is
  likely the dominant runtime cost.
- **Quadratic FASTA writing.**
  [`write_seqrecord_to_fasta`](jcast/helpers.py#L246-L275) re-parses the
  entire output FASTA on every write to de-duplicate, which is O(N²). Keep an
  in-memory `set` of written sequences per output file.
- **Canonical lookup before filters.**
  [jcast/main.py:190](jcast/main.py#L190) fetches the canonical protein
  (UniProt network call on cache miss) before the read-count and FDR filters,
  so low-count junctions pay for it even without `-c`. Filter first; fetch
  lazily.
- **Per-call connections.** Each UniProt lookup opens a new SQLite connection
  and a new `requests.Session`. Hold one of each per run.

---

## 7. Dependencies, packaging & environment

### Dependency table (as of 2026-11-03)

| Dep | setup.py | requirements.txt | Resolves today to | Action |
|---|---|---|---|---|
| biopython | `>=1.78,<2` | `>=1.78` | 1.88 | OK |
| gtfparse | `>=1,<2` | `>=1.2.1` | **3.0.5 (breaks)** | ✅ `>=1.3,<2` (PR 0); polars migration later |
| pandas | `>=1.0,<2` | `>=1.3.0` | 3.1 | ✅ `>=1.5,<4` (PR 0) |
| numpy | *(implicit)* | *(implicit)* | 2.4 | ✅ `>=1.20,<3` declared (PR 0) |
| scipy | *(missing)* | `>=1.7.0` | 1.17 | ✅ added to setup.py (PR 0) |
| scikit-learn | *(commented out)* | `>=1.0` | 1.9 | ✅ lazy-imported; only needed for the vestigial `gaussian_mixture` path (PR 0) |
| matplotlib | `>=3.4,<4` | `==3.4.2` | 3.11 | ✅ `>=3.4,<4` in both; HIGH-5 fixed (PR 0) |
| requests | `>=2,<3` | `>=2.24.0` | — | OK |
| tqdm | `>=4,<5` | `>=4.61.2` | — | OK |
| pomegranate | `>=0.13,<1.0` | `>=0.14.5` | **1.1.2 (breaks)** | ✅ removed (PR 0) |
| platformdirs | — | — | — | ✅ `>=2.5,<5` added for the cache directory (PR 7) |
| PySide6 | extra `gui` | — | 6.11 | OK |

The tiering pipeline itself (everything except `-m`) ran without error on
pandas 3.1 / numpy 2.4 / biopython 1.88 once HIGH-8 was worked around.
Modernizing is a small change.

> **Update — PR 10 (2026-11-13):** every item below is fixed. The metadata
> is in `pyproject.toml`, the Python floor is 3.11, and the dependency
> bounds are tested at both ends on Python 3.11 to 3.15
> ([packaging check](reports/2026-11-13-packaging-check.md)): numpy
> ≥ 1.23.2, scipy ≥ 1.9.2, matplotlib ≥ 3.6, requests ≥ 2.20.
> `requirements.txt` is removed; scikit-learn is in the `test` extra.

### Packaging issues *(new)* — ✅ fixed in PR 10

- `find_packages()` picks up `tests/` (it has an `__init__.py`), so a
  top-level package named `tests` is installed into site-packages. Use
  `find_packages(exclude=['tests*'])`.
- `data_files=[('tests', …)]` installs the ~65 MB chr15 genome/GTF into
  `<sys.prefix>/tests/`. Remove it; `MANIFEST.in` already ships them in the
  sdist.
- `python_requires='>=3.7'` and the 3.6–3.9 classifiers are stale. The GUI
  uses `QWidget | None` (PEP 604), which raises `TypeError` at import on
  < 3.10.
- Version is duplicated in `setup.py` and `jcast/__init__.py`. Move metadata
  to `pyproject.toml` (`[project]`) with a single version source.
- `tox.ini` targets py37–39 only.

### Python floor

Bump to **3.10** for v0.4.0, as previously recommended. This is now also
required by the GUI's annotations. **Done in PR 10, at 3.11:** 3.10 reached
its end of life in October 2026.

---

## 8. Test coverage

> **Update — PR 1 (2026-11-03):** all 18 tests pass. Tests are pytest-style
> with session fixtures in [tests/conftest.py](tests/conftest.py), which load
> the chr15 data once. The GTF cache and UniProt cache go to temp dirs, and
> data paths are anchored at the test file. Network tests carry
> `@pytest.mark.network` and the golden run carries `@pytest.mark.slow`; run
> `pytest -m "not network"` offline. Gaps 1, 2 (marking only; network tests
> still run by default), 5 and 6 below are done. The table and hygiene notes
> describe the state before PR 1.
>
> **Update — 2026-11-12:** 148 tests in 16 files, about 2.5 minutes. Each
> output-changing PR added regression tests that fail on the code before it.
> Gap 4 was done in PR 3 (`tests/test_masking.py`). Gap 7 is PR 11.

| File | Tests | Status (2026-11-03 run) |
|---|---|---|
| [tests/test_translate.py](tests/test_translate.py) | 1 | Passes |
| [tests/test_gtf.py](tests/test_gtf.py) | 2 | **Fail** — missing `logger`; then `str` path vs `.name` |
| [tests/test_junction.py](tests/test_junction.py) | 2 | **Fail** — same; `test_pkm_manual` also needs network (UniProt) |
| [tests/test_annotated.py](tests/test_annotated.py) | 2 | **Fail** — same; also needs network |
| [tests/test_model.py](tests/test_model.py) | 3 | **Import error** (pomegranate), and also missing `logger` |

Other test hygiene problems: tests depend on CWD (`tests/data` relative
paths), `tearDown` deletes `.cached` files from the repo's data directory,
`test_model` writes plots into `tests/out/`, and nothing asserts on
`min_count`.

### Gaps to fill

1. Repair the constructors (after HIGH-8). Use `tmp_path` for outputs and
   caches, and anchor paths at `Path(__file__).parent`.
2. Mark network tests (`@pytest.mark.network`) and default to
   `params.use_gtf_only` elsewhere.
3. ✅ Regression tests pinning `tx0`/`tx1` for one `+` and one `−` transcript
   (the gate for §5.1). Added in PR 4.
4. ✅ `make_pep` / `get_complementary` with lowercase and IUPAC input.
   Added in PR 3.
5. A model test asserting `min_count` within a tolerance band per dataset,
   plus the degenerate-fit fallback (the stomach dataset is a natural fixture).
6. **End-to-end golden test.** Run the chr15 fixture in GTF-only mode and
   assert per-fate counts. Baseline as shipped: FAIL 360, T4 236, SKIP 297,
   T1 185, T2 25, T3 2. This catches silent regressions during refactors.
   ✅ Added in PR 1. Its values are updated with each output-changing PR; the
   history is in [reports/README.md](reports/README.md).
7. GitHub Actions running pytest on 3.10–3.13.

---

## 9. v0.4.0 hardening roadmap (revised)

The order changed. PR 0 is new and must come first, because the package
currently can't be installed and run from scratch.

| # | PR | Theme | Gate |
|---|---|---|---|
| 0 | ✅ **Unblock**: Appendix A model, lazy optional imports, `pd.concat`, drop `seaborn-white`, remove pomegranate, bound `gtfparse<2`, `pandas.errors.UndefinedVariableError` | Env | Done 2026-11-03 |
| 1 | ✅ Paths not file handles (HIGH-8), then repair tests (HIGH-7) + golden e2e test | Tests | Done 2026-11-03 |
| 2 | ✅ `NaN` gene symbols → `NoGeneName_<gene id>` (MED-7 was fixed in PR 4) | Compat | Done 2026-11-11 |
| 3 | ✅ Soft-masking / IUPAC (HIGH-2/3): uppercase at retrieval, skip slices with non-ACGT bases, `--mask` deprecated ([audit](reports/2026-11-11-masked-bases-audit.md)) | Correctness | Done 2026-11-11 |
| 4 | ✅ **`tx0`/`tx1` from start/stop codons** (§5.1) — regression test first; changelog entry | Sci-integrity | Done 2026-11-10 |
| 4a | ✅ **Drop pure-UTR exons**: strand-agnostic `_trim`, `SKIPPED_OUTSIDE_CDS`, second-start retry removed ([audit](reports/2026-11-10-utr-exon-trimming-audit.md)) | Sci-integrity | Done 2026-11-10 |
| 4a2 | ✅ **RI slice composition** (HIGH-11): slice 1 is the retained exon alone ([audit](reports/2026-11-10-ri-slice-audit.md)) | Sci-integrity | Done 2026-11-10 |
| 4b | ✅ **Per-slice frame at the slice's 5′ end**, one-frame translation only; three-frame fallback, Tier 3 and `psq_T3*.fasta` removed (MED-2 closed) ([audit](reports/2026-11-10-phase-audit.md)) | Sci-integrity | Done 2026-11-10 |
| 4b2 | ✅ **Frameshift test from coding lengths.** `Sequence.frameshift` compared raw slice lengths, which is wrong once slices have their own frames (PR 4b). Fixed: 17 chr15 junctions reclassified, no sequence lost or created ([audit](reports/2026-11-10-frameshift-audit.md)) | Sci-integrity (regression) | Done 2026-11-10 |
| 4c | ✅ **Alternative termini** (done 2026-11-10; [audit and verification](reports/2026-11-10-alt-termini-audit.md)). C-terminal: a stop at an annotated protein-coding stop codon is not a PTC; join N-terminally only; write to `psq_T1_alt_cterm.fasta` ([report](reports/2026-11-10-alt-cterm-tiering-audit.md)). N-terminal: a slice lacking the gene's start codon that contains another annotated start is translated from it and joined C-terminally only ([phase audit §5.4](reports/2026-11-10-phase-audit.md)). Re-evaluated on the current pipeline in [the termini audit](reports/2026-11-10-alt-termini-audit.md); decisions approved 2026-11-10 (TSL ≤ 2 sites, one file per terminus, fate `SUCCESS TIER 1 ALT TERMINUS`). | Sci-integrity | Done 2026-11-10 |
| 5 | ✅ **Joining to the canonical** (MED-3): complete ends not joined, unique ordered anchors, anchor length 7 ([audit](reports/2026-11-11-stitch-anchor-audit.md)) | Correctness | Done 2026-11-11 |
| 6 | ✅ Threshold semantics (MED-11) and `sum_sjc` naming (MED-1 / §5.2). Decided 2026-11-12 ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md)): keep if count ≥ threshold; default `-r 1` (no filter); `-m` cut at the model's boundary, reported in reads; the count stays min(IJC, SJC), renamed `min_jc`; the read-count and FDR filters applied first ([log audit](reports/2026-11-12-log-audit.md)) | Sci-integrity | Done 2026-11-12 |
| 7 | ✅ Output dir, caches, idempotent output, logger handlers (MED-4/5/6/12/13); catch an exception per junction, record it as a fate and continue instead of aborting the run ([masked-bases audit §7](reports/2026-11-11-masked-bases-audit.md#7-decisions-requested)). Decided 2026-11-12 ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md)): refuse a directory holding JCAST output unless `--overwrite`; a failing junction writes nothing; caches in the user cache directory (`platformdirs`, `--cache_dir`), GTF cache keyed by path, size and mtime, written atomically, as a pickle | Robustness | Done 2026-11-12 |
| 8 | ✅ Network resilience + `write_canonical` fallback (MED-9/10); `-c` writes the canonical on every path (MED-14). Decided 2026-11-12 ([audit](reports/2026-11-12-thresholds-runs-canonical-audit.md)): `-c` on every path, masked junctions included; a failed API lookup fails the junction (timeout, 429 retried); `write_canonical`'s fallback removed; gene names used for `--uniprot` only when unique in the GTF ([mapping audit](reports/2026-11-12-uniprot-mapping-audit.md)) | Robustness | Done 2026-11-12 |
| 9 | ✅ Per-gene GTF index, in-memory FASTA de-dup ([parallel audit](reports/2026-11-12-parallel-audit.md)); filtering before the canonical lookup done in PR 6 | Performance | Done 2026-11-12 |
| 10 | ✅ Packaging: `pyproject.toml` metadata with the version set once (0.4.0), Python ≥ 3.11, only `jcast` installed, the tests and their data in the sdist, dependency bounds tested at both ends on Python 3.11–3.15 ([packaging check](reports/2026-11-13-packaging-check.md); decided 2026-11-13: Python ≥ 3.11, CI on Linux and macOS, version 0.4.0; `4a039a9`) | Env | Done 2026-11-13 |
| 11 | CI (GitHub Actions) | Infra | PRs 1, 10 |
| 12 | argparse hygiene: drop `--g_or_ln`, `-s` reads only selected files (MED-8), `parser.error` not `assert` | UX | PR 1 |
| 13 | GUI v0 → wizard (§10) | New surface | PRs 7, 12 |
| 14 | ✅ `--uniprot FASTA`: canonical from a supplied UniProt FASTA (no API calls), GTF fallback, built-in Swiss-Prot merge ([evaluation](reports/2026-11-11-uniprot-fasta-evaluation.md)) | Reproducibility | Done 2026-11-11 |
| 15 | ✅ Parallel junctions: `-p/--processes` (fork; workers return records, the main process writes in order); one process without `--uniprot` ([parallel audit](reports/2026-11-12-parallel-audit.md)) | Performance | Done 2026-11-12 |
| 16 | ✅ Run record: per-junction `junctions.tsv`, a small log with an end-of-run summary, sequences only with `--verbose`, no `~module~function` suffix ([log audit](reports/2026-11-12-log-audit.md) §7; whole-genome log 346 MB → 5 kB, output unchanged) | UX | Done 2026-11-12 |
| 17 | ✅ GTF canonical from GENCODE's `Ensembl_canonical` tag when no transcript passes the TSL ≤ 2 / start-and-stop filter: complete protein_coding transcripts, readthrough transcripts excluded ([`Ensembl_canonical` audit](reports/2026-11-12-ensembl-canonical-audit.md) §10, `8ddf8c3`; 17 genes in the whole-genome run, fates and merged database unchanged; 370 of the 371 coding genes still without a canonical are readthrough genes) | Sci-integrity | Done 2026-11-12 |

### Deferred to v0.5.x

- gtfparse 2.x/3.x (polars) migration.
- `Config` dataclass replacing module-state `params.py` (prerequisite for a
  clean GUI parameter panel).
- Split `sequences.py` into translation / canonical / stitch modules.
- gff3 input; APPRIS-tag canonical selection (the `Ensembl_canonical`
  fallback was done in v0.4.0, PR 17); "longest sequence" fallback when
  UniProt fails.
- ~~Optional concurrency.~~ Done in PR 15 (`-p`, with `--uniprot`).
- StringTie input support.

### Deferred to v1.0

- Stable Python API beyond the CLI/GUI.
- Splice-graph-based canonical mapping
  ([jcast/main.py:520](jcast/main.py#L520) `TODO`).
- Reproducibility manifest: input hashes, tool version, effective thresholds
  and model parameters, written to the output dir. **Mostly done in PR 16:**
  the log records the version, the arguments, the inputs (sha256 for the
  rMATS files and the UniProt FASTA; size and modification time for the GTF
  and genome), the model and the cutoff. A machine-readable manifest is
  still open.

---

## 10. GUI plan — notes

The PySide6 plan from the first review (wizard → model preview → run panel,
`QThread` worker calling `run_jcast`, optional `jcast[gui]` extra) still
stands.

**`-p` from the GUI needs a different start method (added 2026-11-12).**
- `-p` forks its workers (PR 15). Forking a process that runs Qt threads
  is unsafe: the child copies locks that other threads may hold.
- So the GUI must keep `processes=1` (the default) when it calls
  `run_jcast` in-process. Alternatively, it can run `python -m jcast` as a
  subprocess, or the pool can use `spawn` or `forkserver` for GUI runs.
- `--uniprot` has no such constraint and should be exposed, with
  `--overwrite`, `--cache_dir` and `--verbose`.

The skeleton in [jcast/gui/](jcast/gui/) imports and launches. Points
the first review missed:

- ~~**Logger handler stacking (MED-13)** must be fixed before the GUI can
  run more than one job per session.~~ Fixed in PR 7.
- **Model-plot embedding** needs the plot functions refactored to draw onto a
  caller-supplied `Axes` instead of calling `plt.subplots` / `savefig`
  internally. Appendix A already returns plain parameters, which makes this
  easy.
- **The draggable cutoff:** MED-11 is settled (PR 6). The cutoff is in
  reads, and the model is fitted on ln(count + 1).
- **Cancellation.** `run_jcast` has no cancellation hook. The worker needs a
  per-junction check of a flag/event.
- `--gui` is only recognized via `python -m jcast`, not by the `jcast`
  console script (`main()` rejects it). Either accept it in argparse or
  document `jcast-gui` as the only entry.

Standalone bundles (PyInstaller / briefcase) stay in v0.5 scope. Dropping
pomegranate (and thus torch) makes those bundles far smaller.

---

## Appendix A — proposed `model.py` replacement

This is a drop-in replacement for `general_mixture_model` and
`plot_general_mixture_model`: the same Gamma+Normal model on ln(count + 1),
fitted by deterministic EM with numpy/scipy only. Tested on the 12 bundled
datasets (results in §2) with numpy 2.4, scipy 1.17 and matplotlib 3.11.
`min_count` keeps its current (count + 1) semantics, so behavior only changes
by the model fit itself; MED-11 is a separate decision.

```python
import math
import os

import numpy as np
import scipy.special as special
import scipy.stats as stats


def general_mixture_model(sum_sjc_array,
                          max_iter=500,
                          tol=1e-6,
                          min_weight=0.05):
    """Model skip junction histogram with a Gamma (trace) + Normal (expressed) mixture, fitted by EM.

        Args:
            sum_sjc_array(numpy array): rMATS skipped-junction count results matrix (counts + 1).
            max_iter(int): Maximum number of EM iterations.
            tol(float): Relative change in log-likelihood at which EM is considered converged.
            min_weight(float): Fit is treated as degenerate if either component weight falls below this.
        Returns:
            ln_sjc(numpy array): Natural log of read counts from rMATS.
            params(tuple): (gamma shape, gamma rate, normal mean, normal sd, gamma weight, normal weight),
                or None if the fit is degenerate.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not
                trace, or None if the fit is degenerate.
    """
    ln_sjc = np.log(np.asarray(sum_sjc_array, dtype=float).ravel())
    x = np.clip(ln_sjc, 1e-6, None)

    # Initialize the gamma on the lower half of the data and the normal on the upper half.
    lo, hi = x[x <= np.median(x)], x[x > np.median(x)]
    shape, rate = lo.mean() ** 2 / lo.var(), lo.mean() / lo.var()
    mu, sd = hi.mean(), hi.std()
    weights = np.array([0.5, 0.5])

    ll_old = -np.inf
    for _ in range(max_iter):
        # E-step
        log_p = np.vstack([np.log(weights[0]) + stats.gamma.logpdf(x, shape, scale=1 / rate),
                           np.log(weights[1]) + stats.norm.logpdf(x, mu, sd)])
        ll = special.logsumexp(log_p, axis=0)
        resp = np.exp(log_p - ll)

        # M-step
        weights = resp.mean(axis=1)
        if weights.min() < min_weight:
            return ln_sjc, None, None
        mu = np.average(x, weights=resp[1])
        sd = math.sqrt(np.average((x - mu) ** 2, weights=resp[1]))
        # Gamma MLE: solve log(k) - digamma(k) = log(mean) - mean(log x), Minka (2002) start + Newton steps
        mean = np.average(x, weights=resp[0])
        s = math.log(mean) - np.average(np.log(x), weights=resp[0])
        if not s > 0:
            return ln_sjc, None, None
        shape = (3 - s + math.sqrt((s - 3) ** 2 + 24 * s)) / (12 * s)
        for _ in range(5):
            shape -= (math.log(shape) - special.digamma(shape) - s) / (1 / shape - special.polygamma(1, shape))
        rate = shape / mean
        if not all(np.isfinite([shape, rate, mu, sd])) or sd <= 0:
            return ln_sjc, None, None

        if ll.sum() - ll_old < tol * abs(ll.sum()):
            break
        ll_old = ll.sum()

    params = (shape, rate, mu, sd, weights[0], weights[1])

    # Decision boundary: where the weighted densities cross between the gamma mode and the normal mean.
    grid = np.linspace(max((shape - 1) / rate, 1e-3), mu, 2000)
    diff = (weights[0] * stats.gamma.pdf(grid, shape, scale=1 / rate)
            - weights[1] * stats.norm.pdf(grid, mu, sd))
    crossings = np.flatnonzero(np.diff(np.sign(diff)))
    if len(crossings) == 0:
        return ln_sjc, None, None
    min_count = math.ceil(np.exp(grid[crossings[0] + 1]))

    return ln_sjc, params, min_count


def plot_general_mixture_model(ln_sjc,
                               params,
                               min_count,
                               write_dir: str,
                               filename: str):
    """Plot the Gamma + Normal mixture fitted to skip junction read counts.

        Args:
            ln_sjc(numpy array): Natural log of read counts from rMATS.
            params(tuple): Output of general_mixture_model.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not trace.
            write_dir(str): Path to directory for plot output.
            filename(str): Name of plot output file.
        Returns:
            None
    """
    import matplotlib.pyplot as plt

    shape, rate, mu, sd, w0, w1 = params
    x = np.linspace(0, 10, 2000)
    eval_dist0 = stats.gamma.pdf(x, shape, scale=1 / rate) * w0
    eval_dist1 = stats.norm.pdf(x, mu, sd) * w1

    fig, ax = plt.subplots(constrained_layout=True, figsize=(10, 4))
    ax.hist(ln_sjc, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    ax.set_xlabel('Natural logarithm of total skipped junction counts')
    ax.set_ylabel('Frequency')
    ax.set_title('Skipped junction count distributions')
    ax.set_ylim(0, 1)

    ax.plot(x, eval_dist0, color="C0")
    ax.plot(x, eval_dist1, color="C3")
    ax.plot(x, eval_dist0 + eval_dist1, color="black")
    ax.fill_between(x, 0, eval_dist0, color="C0", alpha=0.2)
    ax.fill_between(x, 0, eval_dist1, color="C3", alpha=0.2)
    ax.axvline(linewidth=2, ls='--', color='black', x=np.log(min_count))
    ax.text(np.log(min_count) + math.log(1.25), 0.75, "min_counts = " + str(min_count))

    # Second axis for untransformed reads
    tcks = ax.get_xticks()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Untransformed total splice junction read counts')
    ax2.set_xticklabels(np.exp(tcks[1:-1]).astype(int))

    fig.savefig(fname=os.path.join(write_dir, filename + '.pdf'))
    plt.close(fig)
```

Matching change in `run_jcast` ([jcast/main.py:58-65](jcast/main.py#L58-L65)):

```python
ln_sjc, mix_params, min_count = model.general_mixture_model(sum_sjc_array=rmats_results.sum_sjc_array)
if min_count is None:
    logger.warning(f'Mixture model fit was degenerate; falling back to -r {args.read}.')
    min_count = args.read
else:
    logger.info(f'Mixture model parameters {mix_params}; modeled min_count {min_count}.')
    model.plot_general_mixture_model(ln_sjc, mix_params, min_count, write_dir=args.out, filename='model')
```
