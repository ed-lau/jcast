# JCAST — Critical Review & v0.4.0 Hardening Roadmap

*Review date: 2026-11-03. Reviewed version: v0.3.5a (HEAD `89a9dea` plus the
uncommitted `jcast/gui/` skeleton). Supersedes the 2026-05-23 review.*

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
bug is **not silent**, contrary to the first review. It truncates 15 % of SE
and 11 % of MXE junctions on chr15. Correcting it changes the fate of 45 of
1,105 junctions and recovers 12 additional T1 sequences
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

JCAST is a small flat package (~2,500 lines of Python, 12 modules plus the GUI
skeleton) that runs as a single-process CLI.

```
jcast/
├── __main__.py        Entry — dispatches to gui.launch_gui (--gui) or main.main()
├── main.py            CLI parsing, top-level run loop, tier dispatch (_translate_one)
├── junctions.py       RmatsResults reader; Junction (trim_cds, get_translated_phase, sum_sjc)
├── sequences.py       Sequence: nucleotide retrieval, translation, canonical lookup
│                       (GTF or UniProt), stitching, FASTA writing
├── annots.py          ReadAnnotations (GTF), ReadGenome (FASTA), AnnotatedTranscript
├── helpers.py         get_local_nuc, get_complementary, make_pep, write_seqrecord_to_fasta
├── model.py           Mixture-model read-count cutoff (pomegranate + vestigial sklearn path)
├── params.py          Module-level run parameters (TSL threshold, stitch length, …)
├── constants.py       Genetic code table
├── fates.py           Per-junction outcome strings
├── logger.py          Logger setup
└── gui/               Optional PySide6 skeleton (uncommitted)
```

### Data flow

1. **`RmatsResults`** reads all five `*.MATS.JC.txt` files and renames their
   columns positionally.
2. **`ReadAnnotations`** parses the GTF via `gtfparse` and caches a TSV copy
   next to the input.
3. **`ReadGenome`** loads the genome FASTA (optionally gzip) into a `SeqIO`
   dict held in memory.
4. **Optional `-m`**: fit a mixture to ln(count + 1) and use the decision
   boundary as the read cutoff.
5. **Per junction**: `trim_cds` → `get_translated_phase` → slice nucleotides →
   mask check → canonical lookup → read/FDR filters → tier T1–T4 or fail →
   append to `psq_*.fasta`.

### Tier classification

| Tier | Meaning |
|---|---|
| **T1** | Both slices translate without stop codons in the GTF phase; no frameshift. |
| **T2** | Same, but the slice-length difference is not divisible by 3 (frameshift). |
| **T3** | The GTF phase fails, but another frame translates both slices cleanly. |
| **T4** | One slice hits a PTC. The other is force-translated, and kept only if ≥ `ptc_threshold` × the good slice. |
| Skip | Low count, FDR outside range, or masked nucleotides. |
| Fail | None of the above. The canonical is written if `-c` is set. |

Slices that translate but fail to stitch to the canonical go to
`psq_T#_orphan.fasta`.

---

## 4. Bugs & correctness issues

Tags: **HIGH** blocks running on current environments or biases scientific
output. **MED** reduces robustness or gives wrong output in edge cases.
**LOW** is UX or style. ✔ = reproduced by execution; ◎ = confirmed by
reading.

### HIGH-2 — `make_pep` hard-fails on non-ACGT bases ✔

[jcast/helpers.py:226](jcast/helpers.py#L226). `code[nt[i:i+3]]` raises
`KeyError` for lowercase (soft-masked) and IUPAC-ambiguous bases. Reproduced
for `atg…`, `NNN`, and `RAA` on both strands. Fix: `nt = str(nt).upper()`
and `code.get(codon)`, with an explicit policy for unknown codons: either
emit `X`, or treat it like a stop and skip.

### HIGH-3 — `get_complementary` raises before its fallback ✔

[jcast/helpers.py:148-152](jcast/helpers.py#L148-L152). `raise` makes the
`comp += 'X'` line unreachable. Fix together with HIGH-2, ideally using
`Bio.Seq.reverse_complement`, which already handles IUPAC codes and case.

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

### MED-1 — `sum_sjc` returns `min(SJC, IJC)` ◎

[jcast/junctions.py:245-274](jcast/junctions.py#L245-L274). The name, the
README header field 11 ("skipped-junction read count"), and the plot axis
label all say SJC. See §5.2.

### MED-2 — Three-frame fallback keeps the *last* frame; can discard a good one ◎

[jcast/sequences.py:143-151](jcast/sequences.py#L143-L151). There's no
`break`. If frame 0 succeeds and a later frame fails, both slices are reset to
`''` and the junction drops to T4/fail. Add `break` after the first success,
or pick the longest translation deliberately.

### MED-3 — Stitch anchor is not checked for uniqueness or order ◎

[jcast/sequences.py:506-510](jcast/sequences.py#L506-L510). `str.find`
takes the first occurrence of the 10-mer, and nothing checks
`merge0 <= merge1`. Repeats can mis-stitch, and an end anchor found
*before* the start anchor duplicates canonical sequence. Require
`count == 1` for each anchor and `merge0 <= merge1`, or else send the slice to
orphan.

### MED-4 — Output directory never created before the log file opens ◎

[jcast/logger.py:16](jcast/logger.py#L16). Add
`os.makedirs(args.out, exist_ok=True)` before `get_logger`.

### MED-5 — UniProt SQLite cache written to CWD ◎

[jcast/sequences.py:397-400](jcast/sequences.py#L397-L400). Use
`platformdirs.user_cache_dir('jcast')` (machine-wide) or `args.out`. Wrap the
`INSERT` in a try/except for read-only volumes.

### MED-6 — GTF `.cached` written next to the input, keyed only by filename ◎

[jcast/annots.py:45-49](jcast/annots.py#L45-L49). This fails on read-only
reference volumes. It also serves a **stale cache** if the GTF is replaced
under the same name, and `to_csv` writes the index as an unnamed extra column.
Fix: cache in a user dir keyed on (path, size, mtime), with `index=False`.

### MED-7 — `_get_translated_region` exception handling is broken on pandas ≥ 2 ✔ *(expanded)*

[jcast/junctions.py:300-338](jcast/junctions.py#L300-L338).

- `pd.core.computation.ops.UndefinedVariableError` **no longer exists**
  (verified on pandas 3.1). For a GTF without `transcript_support_level`, the
  `except` clause itself raises `AttributeError`. Use
  `pandas.errors.UndefinedVariableError` (available since pandas 1.5), or
  better, check `'transcript_support_level' in gtf.annot.columns` once.
- On the `KeyError` path, `gtf0_start` / `gtf0_end` are set to `None` and then
  passed to `len()`, which raises `TypeError`. Use empty Series instead.

### MED-8 — `RmatsResults` reads all five files regardless of `-s` ◎ *(new)*

[jcast/junctions.py:36-40](jcast/junctions.py#L36-L40) versus the existence
check at [jcast/main.py:437-440](jcast/main.py#L437-L440), which only checks
the selected types. A folder with only `SE.MATS.JC.txt` passes validation
and then raises `FileNotFoundError` (this is the existing `TODO` at main.py:35).
The check also uses `assert`, which is stripped under `python -O`; use
`parser.error`.

### MED-9 — `write_canonical`'s UniProt fallback is a no-op ◎ *(new)*

[jcast/sequences.py:537-543](jcast/sequences.py#L537-L543). It calls
`get_canonical_aa_uniprot()` but discards the return value instead of
assigning `self.canonical_aa`. The fallback never writes anything, yet still
costs a network request, and it does so even when `params.use_gtf_only` is
set.

### MED-10 — A network failure aborts the entire run ◎ *(new)*

[jcast/sequences.py:429-431](jcast/sequences.py#L429-L431). Once its retries
are exhausted, `rqs.get(...)` raises `requests.exceptions.ConnectionError` /
`RetryError`, which nothing catches. A transient EBI outage hours into a run
kills it. Catch the exception, log it, return an empty record, and let the
junction fall through to orphan/fail.

### MED-11 — Read-count thresholds are off by one, twice ◎ *(new)*

- **`-r`**: the help text and README say "minimum read count for a junction
  to be translated", but [jcast/main.py:196](jcast/main.py#L196) skips when
  `sum_sjc <= args.read`. The default `-r 1` therefore requires ≥ 2 reads.
- **`-m`**: `min_count` is computed on the `count + 1` scale
  ([jcast/model.py:63](jcast/model.py#L63)) and then compared to raw counts
  with `<=`. The effective threshold is `min_count + 1` raw reads.

Neither is catastrophic, but both make reported cutoffs disagree with what
was applied. Pick one convention (`keep if count >= threshold`), document it,
and log the effective threshold.

### MED-12 — Output is not idempotent ◎ *(new)*

`psq_*.fasta` files are opened in append mode
([jcast/helpers.py:272](jcast/helpers.py#L272)), as is `logfile.log`.
Re-running into the same `-o` with different parameters silently merges
results across runs. Refuse a non-empty output dir without `--overwrite`, or
write to a temp dir and move it into place.

### MED-13 — `get_logger` stacks handlers on every call ◎ *(new; matters for the GUI)*

[jcast/logger.py:30-31](jcast/logger.py#L30-L31). Each `run_jcast` call adds
another `FileHandler` and `StreamHandler` to the same `'jcast'` logger. In a
long-lived process (the GUI, tests, notebooks), run *n* logs every line *n*
times and keeps writing into earlier runs' log files. Remove existing
handlers first, and close the file handler when the run ends.

### LOW

- [jcast/junctions.py:484-578](jcast/junctions.py#L484-L578)
  `except TypeError or ValueError` should be `(TypeError, ValueError)`
  (was HIGH-1). Only `TypeError` can occur here, so there is no observed
  impact.
- [jcast/helpers.py:232](jcast/helpers.py#L232): the terminal-stop tolerance
  accepts a stop followed by 0–1 trailing nt but rejects 2 trailing nt
  (`'ATGAAATAAG'` → `'MK'`, `'ATGAAATAAGG'` → `''`). Both cases are the same
  biologically.
- [jcast/junctions.py:218](jcast/junctions.py#L218): `gene_symbol is None`
  never fires, because pandas parses `NA` as `NaN`. Use `pd.isna`. The README
  "defensive substitution" for `NA` gene names is therefore ineffective.
- [jcast/helpers.py:29-34](jcast/helpers.py#L29-L34): chromosome-name
  reconciliation only strips `chr`. It never adds it (Ensembl rMATS +
  GENCODE genome → `KeyError`).
- [jcast/__init__.py:6-11](jcast/__init__.py#L6-L11): the glob-built `__all__`
  makes `from jcast import *` import every submodule. Remove it.
- [jcast/main.py:421-425](jcast/main.py#L421-L425): `--g_or_ln` is read by
  nothing. Remove it (from the README too).
- [jcast/sequences.py:153](jcast/sequences.py#L153): stray `print()`.
- [jcast/helpers.py:43-121](jcast/helpers.py#L43-L121): `get_rest_nuc` is
  unused, uses bare `except:`, and has `print`-based logging. Remove it.
- [jcast/model.py:44-50](jcast/model.py#L44-L50): the KDE grid
  `linspace(-5, 20, 500000)` silently clips any parameter outside that range,
  and `gaussian_kde` raises if all restarts agree. This is moot if Appendix A
  is adopted.
- `tests/data/model/encode_mouse_fetalheart/*.MATS.JC.txt` are header-only
  (0 junctions).
- Date-stamped comments (`# 2020-07-25 …`) throughout. Move them to git
  history.

---

## 5. Scientific integrity assessment

JCAST output feeds MS database searches, so silent sequence errors propagate
into identifications.

### 5.1 — `tx1` is derived from `start_codon`, not `stop_codon` — measured

[jcast/junctions.py:322-347](jcast/junctions.py#L322-L347). The block
commented "now getting end codons" queries `feature == "start_codon"`. On
`+` strand, `tx1` becomes the most downstream *start* codon (+2). On `−`
strand, after the swap, the lower bound `tx0` is likewise a start codon.
Alternative transcripts often start translation at internal ATGs, so this
coordinate frequently lies *inside* coding exons.

**Measured on chr15 (1,105 junctions processed):**

| Type | Junctions | End-side trim fired |
|---|---|---|
| SE | 720 | 109 (15 %) |
| MXE | 90 | 10 (11 %) |
| A5SS | 71 | 4 |
| RI | 88 | 2 |
| A3SS | 136 | 2 |

**Output impact.** The full pipeline was run with a prototype fix:
`+`: `[min start_codon.start, max stop_codon.end]`;
`−`: `[min stop_codon.start, max start_codon.end]`; no `+2` offset.

- 45 / 1,105 junctions change fate. 23 FAIL → T4, 8 FAIL → T1, 1 FAIL → T2,
  1 T4 → T2. In the other direction, 4 T4 → FAIL, 3 T1 → FAIL, and 4 T1/T2 →
  T4.
- `psq_T1.fasta`: 152 → 164 records. **12 sequences gained, 0 lost; no
  existing T1 sequence changed.** Stitching back to the canonical masks
  truncated C-termini, so the damage shows up as *lost isoforms* (e.g. a
  skip-form slice reduced to the anchor exon) rather than corrupted ones.
- T4 and orphan files change more (≈ 20 sequences each way).

**Recommendations:**

1. Don't simply swap the feature string. The `+2` offset was tuned for start
   codons, and with stop codons it leaves 2 nt past the stop, which
   `make_pep` then treats as a PTC (see the LOW item above).
2. Pin `tx0`/`tx1` for one `+` and one `−` transcript in a regression test
   *before* changing the code.
3. Investigate the 3 T1 → FAIL transitions in the prototype before shipping.
4. Call out the behavior change in the changelog, so users of ≤ v0.3.5 can
   decide whether to re-run.

### 5.2 — `sum_sjc` is `min(SJC, IJC)` ◎

The cutoff model, the FASTA header (`r<count>`), the plot axis, and the README
all describe this as the skipped-junction count. The v0.3.2 changelog says
"sum of all replicates", and the 2022-05 switch to `min` is undocumented.
Decide on one semantics, rename the property (`min_jc`), and update the
header documentation.

### 5.3 — Soft-masked / ambiguous bases ✔

Covered by HIGH-2/3. The README's workaround (use an unmasked genome) is
accurate.

### 5.4 — Stitch anchors can mis-place sequence ◎

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

---

## 6. Performance *(new section)*

None of this matters for the chr15 fixture, but all of it matters genome-wide
(~10⁵ junctions, ~3 M GTF rows).

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
| PySide6 | extra `gui` | — | 6.11 | OK |

The tiering pipeline itself (everything except `-m`) ran without error on
pandas 3.1 / numpy 2.4 / biopython 1.88 once HIGH-8 was worked around.
Modernizing is a small change.

### Packaging issues *(new)*

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
required by the GUI's annotations.

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
3. Regression tests pinning `tx0`/`tx1` for one `+` and one `−` transcript
   (the gate for §5.1).
4. `make_pep` / `get_complementary` with lowercase and IUPAC input.
5. A model test asserting `min_count` within a tolerance band per dataset,
   plus the degenerate-fit fallback (the stomach dataset is a natural fixture).
6. **End-to-end golden test.** Run the chr15 fixture in GTF-only mode and
   assert per-fate counts. Baseline as shipped: FAIL 360, T4 236, SKIP 297,
   T1 185, T2 25, T3 2. This catches silent regressions during refactors.
7. GitHub Actions running pytest on 3.10–3.13.

---

## 9. v0.4.0 hardening roadmap (revised)

The order changed. PR 0 is new and must come first, because the package
currently can't be installed and run from scratch.

| # | PR | Theme | Gate |
|---|---|---|---|
| 0 | ✅ **Unblock**: Appendix A model, lazy optional imports, `pd.concat`, drop `seaborn-white`, remove pomegranate, bound `gtfparse<2`, `pandas.errors.UndefinedVariableError` | Env | Done 2026-11-03 |
| 1 | ✅ Paths not file handles (HIGH-8), then repair tests (HIGH-7) + golden e2e test | Tests | Done 2026-11-03 |
| 2 | pandas ≥ 2 exception handling (MED-7), `NaN` gene symbols | Compat | PR 1 |
| 3 | Soft-masking / IUPAC (HIGH-2/3) | Correctness | PR 1 |
| 4 | **`tx0`/`tx1` from start/stop codons** (§5.1) — regression test first; changelog entry | Sci-integrity | PRs 1, 2 |
| 5 | Three-frame `break` (MED-2), stitch anchor checks (MED-3) | Correctness | PR 1 |
| 6 | Threshold semantics (MED-11) and `sum_sjc` naming (MED-1 / §5.2) | Sci-integrity | PR 1 |
| 7 | Output dir, caches, idempotent output, logger handlers (MED-4/5/6/12/13) | Robustness | PR 1 |
| 8 | Network resilience + `write_canonical` fallback (MED-9/10) | Robustness | PR 1 |
| 9 | Per-gene GTF index, in-memory FASTA de-dup, filter before canonical lookup | Performance | PR 1 (golden test guards it) |
| 10 | Packaging: pyproject metadata, Python ≥ 3.10, fix `find_packages`/`data_files`, bounds | Env | PR 0 |
| 11 | CI (GitHub Actions) | Infra | PRs 1, 10 |
| 12 | argparse hygiene: drop `--g_or_ln`, `-s` reads only selected files (MED-8), `parser.error` not `assert` | UX | PR 1 |
| 13 | GUI v0 → wizard (§10) | New surface | PRs 7, 12 |

### Deferred to v0.5.x

- gtfparse 2.x/3.x (polars) migration.
- `Config` dataclass replacing module-state `params.py` (prerequisite for a
  clean GUI parameter panel).
- Split `sequences.py` into translation / canonical / stitch modules.
- gff3 input; APPRIS-tag canonical selection; "longest sequence" fallback when
  UniProt fails.
- Optional concurrency. Note that FASTA appends and the shared SQLite cache
  are not thread-safe today.
- StringTie input support.

### Deferred to v1.0

- Stable Python API beyond the CLI/GUI.
- Splice-graph-based canonical mapping
  ([jcast/main.py:224](jcast/main.py#L224) `TODO`).
- Reproducibility manifest: input hashes, tool version, effective thresholds
  and model parameters, written to the output dir.

---

## 10. GUI plan — notes

The PySide6 plan from the first review (wizard → model preview → run panel,
`QThread` worker calling `run_jcast`, optional `jcast[gui]` extra) still
stands. The skeleton in [jcast/gui/](jcast/gui/) imports and launches. Points
the first review missed:

- **Logger handler stacking (MED-13)** must be fixed before the GUI can run
  more than one job per session.
- **Model-plot embedding** needs the plot functions refactored to draw onto a
  caller-supplied `Axes` instead of calling `plt.subplots` / `savefig`
  internally. Appendix A already returns plain parameters, which makes this
  easy.
- **The draggable cutoff** only makes sense once MED-11 settles which scale
  the cutoff is on.
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
