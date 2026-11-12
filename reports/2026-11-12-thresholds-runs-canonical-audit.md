# Audit: read-count thresholds, robustness of real runs, and canonicals under `-c` (PRs 6–8)

| | |
|---|---|
| **Date** | 2026-11-12 |
| **Status** | Proposals in §2.4, §3.6 and §4.4; decisions requested in §6 and taken 2026-11-12 (§7), item 9 pending. **PR 7 implemented and verified (§8).** PRs 6 and 8 to follow. The audit was measured at `8b845aa`. |
| **Question** | Do the read-count cutoffs apply what they say (PR 6)? What goes wrong in real runs: output directory, reruns, one failing junction, the logger, the caches (PR 7)? Does `-c` write every gene's canonical, and what does a UniProt API failure do to a run (PR 8)? |
| **Reproduce** | `python reports/runs_audit/audit_runs.py [--rmats <dir> --jcast_out <dir> --uniprot <fasta>]` (≈ 10 min on chr15; the whole-genome arguments add read counts and `-c` coverage from an existing run) · `python reports/runs_audit/bench_gtf_cache.py --cached <gtf>.cached --scratch <dir>` |
| **Inputs** | The chr15 test set (1,105 junctions); the 12 rMATS sets in `tests/data/model/`; the user's whole-genome run of 2026-11-12 (NMDI RNA-seq, 599,644 junctions; GENCODE v50; Swiss-Prot 2026_04; `-m -c --uniprot -p 8`), not bundled. |
| **Evidence** | [runs_audit/audit/summary.txt](runs_audit/audit/summary.txt) · [thresholds.tsv](runs_audit/audit/thresholds.tsv) · [threshold_fates.tsv](runs_audit/audit/threshold_fates.tsv) · [runs.txt](runs_audit/audit/runs.txt) · [canonical.tsv](runs_audit/audit/canonical.tsv) · [canonical.txt](runs_audit/audit/canonical.txt) · [gtf_cache.txt](runs_audit/audit/gtf_cache.txt) |

## 1. Summary

- **PR 6: both read-count cutoffs keep fewer junctions than they say.**
  - `-r N` is documented as the lowest count translated. The code
    translates only counts above N, so the default `-r 1` needs 2 reads.
  - `-m` keeps junctions two reads above the model's own boundary. On the
    user's whole-genome run, that dropped 3,410 junctions: 2.6% of those
    at or above the boundary. Across the 11 bundled data sets with a model
    fit, it dropped 2–20%.
  - The count is the smaller of the summed inclusion and skipping
    junction reads. The README, the option help, the plot and the
    property name all say skipping reads.
  - No junction in any data set has a count of 0. A literal `-r 1` would
    therefore filter nothing, so the default should become `-r 2`. The
    default output is then unchanged.
- **PR 7: real runs fail in six ways.** Each is reproduced on chr15 in
  [runs.txt](runs_audit/audit/runs.txt):
  - The output directory must already exist, even the default `out/`.
  - **A rerun into the same directory keeps the first run's output.** A
    `-r 30` run into a `-r 1` directory leaves files byte-identical to the
    `-r 1` run.
  - **One exception aborts the run.** Output is partial and the
    `--uniprot` outputs are missing, with or without `-p`.
  - A second run in the same process logs every line twice, once into the
    first run's log.
  - **A half-written GTF cache is read without error.** A second run
    started while the first writes the cache sees only part of the GTF
    (measured with half the file: half the genes). The
    cache also goes stale when a GTF is replaced under the same name.
  - Both caches fail on read-only directories.
- **PR 8: `-c` misses a few canonicals.**
  - 1 chr15 gene and 5 of 13,622 whole-genome genes get no canonical.
  - **Masked junctions never write the canonical:** it is looked up only
    after the masked-base check.
  - `write_canonical`'s fallback sends two API requests and ignores the
    answers: 346 requests on chr15.
- **PR 8: an API failure aborts the run, or silently changes it.**
  - **No network:** the run aborts at the first lookup (MED-10).
  - **Any refusal other than a server error** (for example HTTP 429, too
    many requests) is treated as "no canonical". The run completes, and
    every joined record moves to an orphan file: chr15 `psq_T1` 149 → 0,
    `psq_canonical` 183 → 0.
  - **Requests have no timeout,** so a stalled connection hangs the run.
  - `--uniprot` avoids all three.

## 2. PR 6: read-count thresholds (MED-11, MED-1 / §5.2)

### 2.1 The count

`Junction.sum_sjc` ([jcast/junctions.py:250](../jcast/junctions.py#L250))
returns the smaller of two sums, each over both samples and their
replicates: the inclusion junction counts (IJC) and the skipping junction
counts (SJC). It has done so since 2022-05.

Everything else describes it as the skipping count:
- the property name;
- the `-r` help text;
- README header field 11 ("Skipped-junction read count");
- the model plot's axis label.

The 2022-05 change is not in the changelog.

The minimum asks both isoforms to have read support. On chr15, 75
junctions with at least 2 skipping reads have fewer than 2 inclusion
reads. Whole genome: 18,859.

### 2.2 `-r`

[main.py:328](../jcast/main.py#L328) skips a junction when
`count <= args.read`. The help text and the README say `-r` is the lowest
count translated. `-r N` therefore needs N + 1 reads.

**No junction in any data set has a count of 0.** The minimum is 1 in the
chr15 set, in all 12 bundled sets and in the whole-genome set
([thresholds.tsv](runs_audit/audit/thresholds.tsv)). If `-r 1` did what it
says, the default would filter nothing. In practice, the default removes
every junction with exactly 1 read:

| | Junctions | Count 1 (skipped by the default) |
|---|---|---|
| chr15 | 1,105 | 297 (27%) |
| 12 bundled sets | 8,553 – 116,887 each | 14–28% each |
| Whole genome | 599,644 | 83,534 (14%) |

On chr15, translating those 297 junctions (`-r 1` read literally) adds:
- 81 T1, 19 T2 and 114 T4 fates;
- 65 `psq_T1` records.

It also removes 2 canonical records. Genes whose junctions were all
skipped wrote their canonical on the skip path; once translated, they hit
MED-14 (§4.1) ([threshold_fates.tsv](runs_audit/audit/threshold_fates.tsv)).

### 2.3 `-m`

The mixture model is fitted to ln(count + 1)
([junctions.py:163](../jcast/junctions.py#L163)):
- It returns `min_count` = ⌈e^b⌉, where b is the boundary between the two
  components ([model.py:78](../jcast/model.py#L78)).
- That value is on the count + 1 scale. The boundary in reads is
  `min_count − 1`.
- `_translate_one` then skips `count <= min_count`. So it keeps
  `count ≥ min_count + 1`: **two reads above the boundary.**
- The log ("modeled min_count 41") and the plot ("min_counts = 41", on an
  axis labelled as untransformed read counts) both show the count + 1
  value.

| Data set | Model `min_count` | Kept now (count ≥ min_count + 1) | Kept at the boundary (count ≥ min_count − 1) | Lost to the offset (share of those at the boundary) |
|---|---|---|---|---|
| chr15 | 27 | 114 | 122 | 8 |
| **Whole genome (user's run)** | **41** | **127,823** | **131,233** | **3,410 (2.6%)** |
| encode_human_spleen | 6 | 12,369 | 15,480 | 3,111 (20%) |
| kidney_old_totalRNA_gencodeM26 | 5 | 15,533 | 18,616 | 3,083 (17%) |
| encode_human_lung | 8 | 11,808 | 13,781 | 1,973 (14%) |
| heart_mouse_old | 30 | 31,457 | 32,108 | 651 (2%) |

All 12 sets are in [thresholds.tsv](runs_audit/audit/thresholds.tsv). The
other seven lose 5–14%. `encode_human_stomach` has no stable fit, and the
run falls back to `-r` (REVIEW.md §5.6). The lower the cutoff, the larger
the share lost.

The whole-genome kept count matches the user's run: 599,644 − 471,821
low-coverage skips = 127,823.

### 2.4 Proposal

1. **One rule for both:** keep a junction if its count ≥ the threshold.
   - The same applies to the `-r` fallback when the model fit is
     degenerate.
   - Log the effective threshold in reads at the start of every run.
2. **`-r` default 2.** Default output stays as it is, and `-r` means what
   it says. Scripts that pass `-r N` explicitly now keep junctions with
   exactly N reads: one read fewer than before (CHANGELOG).
3. **`-m` cut at the model's boundary.**
   - `general_mixture_model` returns the cutoff in reads, ⌈e^b⌉ − 1.
   - The log and the plot report that value.
   - Whole genome: 3,410 more junctions translated.
4. **The count stays the minimum** of summed IJC and SJC. It is the
   stricter rule, and both slices are written, so both isoforms need
   reads.
   - Rename `sum_sjc` to `min_jc` (REVIEW.md §5.2), and
     `sum_sjc_array` with it.
   - Correct the `-r` help text, README header field 11, the plot label
     and the CHANGELOG. The value in headers (`r<count>`) does not change.
5. **The golden test** calls `_translate_one` with `read=1`. That becomes
   `read=2`, and its values stay as they are.

## 3. PR 7: robustness of real runs (MED-4, 5, 6, 12, 13; per-junction exceptions)

Each item is measured on chr15 through `run_jcast` with `--uniprot`
([runs.txt](runs_audit/audit/runs.txt)).

### 3.1 The output directory must exist (MED-4)

`get_logger` opens `<out>/logfile.log` before anything creates `<out>`,
which raises `FileNotFoundError`. This includes the default, `-o out`.

### 3.2 A rerun mixes runs (MED-12)

Every `psq_*.fasta` is opened for appending
([helpers.py:291](../jcast/helpers.py#L291)). A sequence already in the
file is not written again.

Measured: chr15 with `-r 1`, then with `-r 30` into the same directory.
- **The output files afterwards are byte-identical to the `-r 1` run
  alone.**
- A `-r 30` run on its own writes 40 `psq_T1` records, not 186.
- The log holds both runs (2,210 junctions).
- The `--uniprot` outputs (`uniprot_canonicals.tsv`, `psq_isoforms.tsv`,
  `jcast_uniprot_merged.fasta`) are rewritten, but built from the mixed
  `psq_*` files.

Nothing warns. A rerun with different settings, a different GTF or a
different UniProt release looks complete.

The user's whole-genome output directory holds a single run (its log has
one start line), so it is not affected.

### 3.3 One exception aborts the run (masked-bases audit §7, item 3)

To simulate one failing junction, the 500th junction of 1,105 was made to
raise.
- **`-p 1`:** the run stops at junction 500, with 410 of 847 records
  written.
- **`-p 3`:** the run stops at junction 490, with 405 records written.
  Workers take junctions in batches of 20, and the junctions before the
  failing one in its batch are lost with it.
- In both, the `--uniprot` outputs are not written.

The other ways to abort a run (§3.5, §4.3) end the same way.

### 3.4 Logger handlers stack (MED-13)

A second `run_jcast` in the same process adds a second pair of handlers.
- Its 10,813 log lines go into its own log and into the first run's log
  (now 21,626 lines).
- Console warnings appear twice.

This matters for the GUI, notebooks and tests. The tests clear the
handlers by hand today.

### 3.5 The caches (MED-5, MED-6)

**UniProt API cache** (`cache/uniprot-cache.db`):
- It is written to the working directory
  ([sequences.py:558](../jcast/sequences.py#L558)).
- From a read-only working directory, the first API lookup raises
  `PermissionError`.
- With `--uniprot` it is never used.

**GTF cache** (`<gtf>.cached`, [annots.py:44-48](../jcast/annots.py#L44-L48)):
- **Written next to the GTF.** A GTF in a read-only directory with no
  cache yet raises `PermissionError`. This is common for shared
  references.
- **Keyed by file name only.** A GTF replaced under the same name is read
  from the old cache: 2,155 genes instead of the new file's 1,152.
- **Not written atomically.** A run that starts while another is writing
  the cache reads the partial file without error: 48,431 of 97,275 rows
  and 1,139 of 2,155 genes. **This is the worst of the six,** because it
  gives wrong output with no error. It can happen whenever two JCAST runs
  on a new GTF are started together, for example one per data set.
- `to_csv` writes the index as an extra `Unnamed: 0` column. It is
  harmless.

**The cache's format affects output.** JCAST always reads the GTF back
from the tab-separated cache, even on the run that writes it. That
round trip sets the column types:
- `exon_number` arrives from gtfparse as text ('2', '10') and becomes a
  number. JCAST sorts exons by it, and as text '10' sorts before '2'.
- `transcript_support_level` 'NA' becomes missing, so
  `transcript_support_level <= 2` can be evaluated.
- `frame` becomes an integer.

Any change to the cache must keep these types.

**Speed.** Reading the GENCODE v50 cache takes 49 s per run (2.7 GB). A
pickle of the same parsed table reads in 8.5 s (1.9 GB), and the table and
its column types come back identical
([gtf_cache.txt](runs_audit/audit/gtf_cache.txt)). An earlier run of the
same measurement gave 54 s and 6.6 s. The cache was read from the shared
volume and the pickle from the local disk.

### 3.6 Proposal

1. **Create the output directory** (`os.makedirs(out, exist_ok=True)`)
   before the logger opens.
2. **Fresh output per run.**
   - If the output directory already holds JCAST output, refuse to run
     unless `--overwrite` is given. JCAST output means `logfile.log`,
     `psq_*.fasta`, `model.pdf`, `jcast_uniprot_merged.fasta`,
     `psq_isoforms.tsv` and `uniprot_canonicals.tsv`.
   - `--overwrite` deletes those files only. Other files in the directory
     are left alone and do not block a run.
   - The in-memory de-duplication (PR 9) then starts empty each run and no
     longer needs to re-read files.
3. **One failing junction does not stop the run.**
   - The junction gets a new fate, `FAILURE (exception). <type>:
     <message>`.
   - Its traceback goes to the log at ERROR level, and the run continues.
   - **A failing junction writes nothing.** Every junction's records are
     collected and written when it completes. Serial runs then work like
     `-p` workers, and `-p` loses only the failing junction, not its
     batch.
   - The console shows the first failure. At the end, a warning gives the
     count of failed junctions and points to the log. The exit status
     stays 0.
   - `KeyboardInterrupt` still stops the run.
4. **Logger:** `get_logger` closes and removes the `jcast` logger's
   existing handlers. `run_jcast` closes its file handler when it ends,
   including on error. The tests stop clearing handlers by hand.
5. **Caches** go in a user cache directory: `~/Library/Caches/jcast` on
   macOS, `~/.cache/jcast` on Linux. They use `platformdirs`, a small
   pure-Python package and a new dependency.
   - `--cache_dir DIR` overrides the location, for clusters with small
     home quotas.
   - **GTF cache:**
     - Keyed by the GTF's absolute path, size and modification time.
     - Written under a temporary name and renamed into place, so a second
       run never reads a partial file.
     - Stored as a pickle of the table after the same tab-separated round
       trip, so the column types are unchanged.
     - The key includes the pandas version, because a pickle from another
       version is rebuilt.
     - Writing a GTF's new cache removes that GTF's older one.
     - The old `<gtf>.cached` files are no longer read. The first run
       after the change rebuilds the cache once.
   - **API cache:** the same directory.
   - **If a cache cannot be written:** log a warning and run without it.
6. **Golden tests:** unchanged.

## 4. PR 8: canonicals under `-c`, and UniProt API failures (MED-9, MED-10, MED-14)

### 4.1 `-c` misses canonicals (MED-14)

With `-c`, the skip, FAILURE and alternative-terminus paths write the
gene's canonical. T1, T2 and T4 write it only when a joined slice equals
it ([sequences.py:736](../jcast/sequences.py#L736)).

Genes with a canonical but no canonical record
([canonical.tsv](runs_audit/audit/canonical.tsv)):

| Run | Genes with a canonical | Without a record |
|---|---|---|
| chr15, GTF canonicals | 269 | 1 (TSPAN3: one T1 junction) |
| chr15, `--uniprot` | 276 | 1 (TSPAN3) |
| Whole genome, user's run (`-c --uniprot`) | 13,622 | 5 (4 Swiss-Prot canonicals, 1 from the GTF) |

MED-14 counted 2 chr15 genes when it was found, at PR 4c; at `8b845aa`
it is 1. The count grows as more junctions are translated: reading `-r 1`
literally would leave 2 more (§2.2).

**Masked junctions never write the canonical.**
- The masked-base check ([main.py:313](../jcast/main.py#L313)) comes
  before the canonical lookup ([main.py:322](../jcast/main.py#L322)), so
  `write_canonical` sees an empty canonical.
- Example: PIAS1's only chr15 junction, with N in its skipped exon. With
  `-c` it writes no canonical, either with GTF canonicals or with
  `--uniprot`.
- On the API path, it sends two requests and ignores the answers (§4.2).

### 4.2 `write_canonical`'s fallback does nothing (MED-9)

When the canonical is empty and `--uniprot` is not given,
`write_canonical` asks the API for it again
([sequences.py:672](../jcast/sequences.py#L672)):
- It asks twice, reviewed and then unreviewed.
- It ignores both answers, so it never writes anything.
- It does so in GTF-only mode too.

On the chr15 API path with `-c`, these calls are 346 of 1,110 requests.
`get_canonical_aa` has already asked the same two questions for that
junction. Removing the fallback changes no output.

### 4.3 API failures (MED-10)

The API path retries connection errors and HTTP 500/502/503/504 ten
times, about 100 s with backoff
([sequences.py:587](../jcast/sequences.py#L587)). Then:

- **Connection refused, or the retries run out:** the exception aborts the
  run at the first lookup (measured with a stub). Any junction can trigger
  this, hours into a run.
- **Any other non-200 answer** (429 too many requests, 403, 404, …): a
  warning is logged and the canonical is left empty. The junction is
  translated without it.
  - Measured with every request refused with HTTP 429: the run completes,
    but every joined record goes to an orphan file instead.

    | File | API answering | Every request refused (HTTP 429) |
    |---|---|---|
    | `psq_T1.fasta` | 149 | 0 |
    | `psq_T1_orphan.fasta` | 101 | 397 |
    | `psq_T4.fasta` | 186 | 0 |
    | `psq_T4_orphan.fasta` | 180 | 520 |
    | `psq_canonical.fasta` | 183 | 0 |

  - **A partial outage moves the affected genes' records the same way,
    with no error.** The output then depends on the network on the day.
- **No timeout** on `rqs.get` ([sequences.py:591](../jcast/sequences.py#L591)).
  A connection that stalls without closing blocks the run forever.

In these runs, the "API answering" baseline is a stub serving the local
API cache, with an empty answer for genes not in it. Only the differences
matter here, not the absolute counts.

### 4.4 Proposal

1. **`-c` writes the gene's canonical on every path.**
   - On T1, T2 and T4, whether or not a slice joins to it exactly.
   - On the masked path too: with `-c`, look up the canonical before the
     masked-base check.
   - Output: 1 more chr15 canonical record, and 5 whole-genome. The golden
     values change: `psq_canonical` 268 → 269, and the `--uniprot` test's
     275 → 276.
2. **Remove `write_canonical`'s API fallback.** No output change. With
   GTF canonicals, JCAST then makes no API requests.
3. **A failed lookup fails the junction, never the canonical.**
   - Any answer other than 200 is a failure. So is a connection error
     after the retries.
   - The lookup raises. With PR 7 (§3.6 item 3), the junction gets the
     exception fate, writes nothing, and the run continues.
   - The end-of-run warning counts those junctions and suggests
     `--uniprot`.
   - Add a 30 s timeout per request, and retry 429. urllib3 honours the
     server's `Retry-After`.
   - Later junctions of the same gene try again. A long outage slows the
     run (each lookup spends about 100 s retrying) but does not stop it.
4. **API path:** otherwise unchanged. `--uniprot` remains the recommended
   route for whole-genome runs.

## 5. Order and verification

- **PR 7 first,** because PR 8's failure handling uses its per-junction
  exceptions. Then PR 6, then PR 8. Each is its own commit, verified
  against this audit as before, with the code under test frozen in a
  worktree.

| PR | Expected change | Verification |
|---|---|---|
| 7 | No change to output | chr15 golden tests unchanged. A whole-genome sample (`verify_runs.py`) byte-identical to `8b845aa` with `-p 1` and `-p 16`. New tests: a missing output directory is created; a non-empty one is refused, or cleared with `--overwrite`; an injected exception gives one exception fate and the same output for every other junction, with and without `-p`; two runs in one process log separately; a stale or partial GTF cache is not used. |
| 6 | Default `-r` unchanged; `-m` keeps more | Golden values unchanged at `read=2`. Whole genome with `-m`: 131,233 junctions pass, exactly the 3,410 more listed by this audit. Tests: `-r N` keeps count N; the model's cutoff in reads. |
| 8 | 1 more chr15 canonical; 5 whole-genome | Golden `psq_canonical` 269 (GTF) and 276 (`--uniprot`); no other file changes. Whole-genome `-c` coverage 13,622 of 13,622. Tests: masked junction with `-c` writes the canonical; API stub refusing (429) or down gives exception fates, not orphans; no API request in GTF-only mode. |

## 6. Decisions requested

**PR 6**
1. Keep a junction if its count ≥ the threshold, for `-r` and `-m`
   (recommended)?
2. Default `-r`: **2** (recommended: default output unchanged), or 1 (no
   filter, since every junction has at least 1 read; chr15 + 297
   junctions, + 65 `psq_T1` records)?
3. `-m`: cut at the model's own boundary, and report the cutoff in reads
   (recommended)? Whole genome: + 3,410 junctions.
4. The count: keep the minimum of summed IJC and SJC, rename it `min_jc`,
   and correct the docs (recommended)? Or return to summed SJC as before
   2022-05 (chr15 + 75 junctions at `-r 1`, whole genome + 18,859)?

**PR 7**

5. Fresh output: (a) **refuse a directory holding JCAST output unless
   `--overwrite`** (recommended). (b) Always replace JCAST's files, with a
   warning (convenient, but a mistyped `-o` destroys a previous run).
   (c) Write to a temporary directory and move it into place at the end.
6. Per-junction exceptions (approved in principle, 2026-11-11). The
   details:
   - a failing junction writes nothing;
   - fate `FAILURE (exception)`;
   - traceback in the log;
   - the first failure on the console, then a count at the end;
   - exit status 0.

   Agree? Or should any failed junction give a non-zero exit status? That
   would make workflow managers such as Snakemake discard the whole
   output.
7. Caches:
   - in the user cache directory via `platformdirs` (new dependency),
     with `--cache_dir`;
   - the GTF cache keyed by path, size and modification time, written
     atomically, as a pickle (49 s → 8.5 s per whole-genome run)?

   Or keep the tab-separated format? It is identical in content and needs
   no pickle, but takes about 40 s longer per run.

**PR 8**

8. `-c` writes the canonical on every path, including masked junctions
   (recommended)?
9. Remove `write_canonical`'s API fallback (recommended; no output
   change)?
10. API failures:
    - (a) **fail the junction and continue**, with a timeout and 429
      retried (recommended);
    - (b) abort the run with a clear message;
    - (c) treat a failure as "no canonical", MED-10's original proposal.
      Not recommended: measured, it silently moves records to orphan
      files.

## 7. Decisions (2026-11-12)

1. **One rule for `-r` and `-m`:** keep a junction if its count ≥ the
   threshold. Approved as recommended.
2. **Default `-r 1`, so the default filters nothing.** This differs from
   the recommendation. `-m` is expected to be the usual choice.
   - The default output then gains the junctions with 1 read: on chr15,
     297 more junctions translated and 65 more `psq_T1` records
     ([threshold_fates.tsv](runs_audit/audit/threshold_fates.tsv),
     second column).
   - The golden values change in PR 6. It loses 2 canonical records
     (MED-14, §4.1), which PR 8 is expected to restore.
3. **`-m` is cut at the model's own boundary,** with the cutoff reported
   in reads. Approved.
4. **The count stays the minimum** of summed IJC and SJC, renamed
   `min_jc`, and the documentation is corrected. Approved.
5. **Fresh output:** a directory holding JCAST output is refused unless
   `--overwrite`, which deletes only JCAST's files. Approved.
6. **A failing junction** is skipped and recorded, and the run continues,
   as proposed:
   - it writes nothing;
   - its fate is `FAILURE (exception)`;
   - its traceback goes to the log;
   - the console shows the first failure and a count at the end;
   - the exit status stays 0.
7. **Caches** go in the user cache directory via `platformdirs`, with
   `--cache_dir`. The GTF cache is keyed by path, size and modification
   time, written atomically, and stored as a pickle. Approved.
8. **`-c` writes the canonical on every path,** including masked
   junctions. Approved.
9. **`write_canonical`'s API fallback:** the maintainer asked how it
   works. It was explained (see below), and the decision is pending.
10. **A failed API lookup fails the junction loudly:** it is recorded
    with the exception fate, and the run continues. A 30 s timeout is
    added and 429 is retried. Approved.

Two further questions were raised with these decisions and are audited
separately:
- the size and content of the log
  ([log audit](2026-11-12-log-audit.md));
- where the Ensembl-to-UniProt mapping under `--uniprot` comes from, and
  how complete it is next to the API
  ([mapping audit](2026-11-12-uniprot-mapping-audit.md)).

**Item 9, as explained.**
- `write_canonical` runs only with `-c`, on the paths that write the
  gene's canonical. When the junction's canonical is empty and
  `--uniprot` is not given, it asks the API again: first for a reviewed
  entry, then for an unreviewed one. It means to write whatever it gets.
- It never assigns the answer to `canonical_aa`, so the empty canonical
  is still empty and nothing is written.
- **When the canonical is empty at that point:**
  - **(i)** On the masked path, the canonical has not been looked up yet.
    Decision 8 moves the lookup before the masked-base check, which
    removes this case.
  - **(ii)** On the API path, `get_canonical_aa` has already asked the
    same two questions for this junction and got nothing. Asking again
    gives the same answer.
  - **(iii)** In GTF-only mode (`params.use_gtf_only`), the gene has no
    GTF canonical. The fallback calls the API although the mode
    promises not to.
- If the fallback were repaired instead of removed, it would change
  output only in case (iii). GTF-only runs would then mix in UniProt
  canonicals for those genes.
- Removing it changes no output and saves the requests: 346 of 1,110 on
  chr15 (§4.2).

## 8. Implementation and verification: PR 7

**What changed** (no change to output):

- **Output directory** ([jcast/main.py](../jcast/main.py)):
  - `run_jcast` creates `-o`.
  - It refuses a directory holding JCAST output (`logfile.log`,
    `psq_*.fasta`, `model.pdf`, the `--uniprot` outputs) unless
    `--overwrite` is given. That option removes those files first, and
    the log is opened for writing, not appending.
  - The CLI reports the refusal through `parser.error`.
- **Per-junction exceptions:**
  - `_translate_safely` wraps `_translate_one` for both the serial path
    and the `-p` workers.
  - It collects the junction's records, and the main loop writes them
    once the junction completes. Serial runs now write exactly as `-p`
    workers do.
  - On an exception:
    - the traceback is logged at ERROR to the log file only;
    - the fate is `FAILURE (exception). <type>: <message>`
      ([jcast/fates.py](../jcast/fates.py));
    - nothing of the junction is written.
  - The console shows the first failure, and a warning at the end counts
    them.
- **Logger** ([jcast/logger.py](../jcast/logger.py)):
  - `get_logger` closes and removes earlier handlers.
  - `close_logger` closes the run's handlers when `run_jcast` ends, also
    on error.
  - A record tagged `jcast_console=False` goes to the file only.
- **Caches** ([jcast/cache.py](../jcast/cache.py)):
  - The user cache directory via `platformdirs` (new dependency), or
    `--cache_dir`.
  - **GTF** ([jcast/annots.py](../jcast/annots.py)):
    - The table is still the GTF written as a tab-separated file and read
      back, now through a temporary file, and without the index column.
    - It is kept as a pickle named after the GTF's path digest, and its
      path, size, mtime and pandas version digest.
    - It is written to `<name>.<pid>.tmp` and renamed into place.
    - Older caches of the same GTF are removed.
    - A cache that cannot be read is rebuilt, with a warning.
  - **API** ([jcast/sequences.py](../jcast/sequences.py)): the same
    directory. A cache that cannot be opened or written gives one warning,
    and the lookup proceeds without it.
- **Tests:** 13 new in `tests/test_runs.py`, 95 in all.
  - The conftest points the caches at a session temp directory. Tests no
    longer clear logger handlers by hand.
  - A test checks that the cached table equals the old `<gtf>.cached`
    table without its index column.

**Verification.** The code under test was frozen in git worktrees: before
= `82a309b`; after = `82a309b` plus this change.

- **chr15** ([runs_audit/verify_7/verify_7.txt](runs_audit/verify_7/verify_7.txt),
  `verify_7.py`): with `-p 1` and `-p 3`, all 12 output files are byte
  for byte identical before and after, and so is the log from the first
  junction.

  The failure modes of §3, rerun on the new code:

  | §3 | Before (§3, runs.txt) | After |
  |---|---|---|
  | 3.1 missing `-o` | `FileNotFoundError` | created; 847 records written |
  | 3.2 `-r 30` into a `-r 1` directory | output identical to the `-r 1` run, 2 runs in one log | refused; with `--overwrite`, `psq` files identical to `-r 30` alone, 1 run in the log |
  | 3.3 exception in junction 500 | run aborted at junction 500 (`-p 1`) or 490 (`-p 3`); 410 or 405 of 847 records; no `--uniprot` outputs | run completed: 1,105 junctions, 1 exception fate, 846 of 847 records, `--uniprot` outputs written; `-p 1` and `-p 3` identical |
  | 3.4 two runs in one process | run 2 also wrote into run 1's log (10,813 → 21,626 lines) | run 1's log unchanged; no handlers left after the runs |
  | 3.5 API cache | under the working directory; `PermissionError` from a read-only one | in the cache directory; no error from a read-only working directory |
  | 3.5 GTF cache next to the GTF | yes, with an index column | nothing next to the GTF; no index column |
  | 3.5 GTF replaced under the same name | stale cache read (2,155 genes, not 1,152) | rebuilt: 1,152 genes; one cache per GTF |
  | 3.5 partial cache | read without error: 48,431 of 97,275 rows | another run's partial file is never read; a damaged cache is rebuilt (97,275 rows) |
  | 3.5 GTF in a read-only directory | `PermissionError` | read without error |

- **Whole genome**
  ([runs_audit/genome_verify_7/verify.txt](runs_audit/genome_verify_7/verify.txt),
  `reports/parallel_audit/verify_runs.py`):
  - 20,000 junctions sampled at random (seed 5) from the user's run, through
    the CLI with `-c -r 1 --uniprot`.
  - `HOME` pointed at a scratch directory, so the cache was built there.
  - **All 15 output files are byte for byte identical before and after,
    with `-p 1` and with `-p 16`.** So are the log messages from the first
    junction on.
  - The start-up messages differ only by the new cache lines.

  | Run | Wall time | Junction loop | Start-up (wall − loop) |
  |---|---|---|---|
  | before, 1 process | 425 s | 341 s | 84 s |
  | after, 1 process, building the GTF cache | 654 s | 308 s | 346 s (the cache took 4 min 56 s to build) |
  | after, 16 processes, reading the cache | 63 s | 27 s | 36 s |

  Building the cache from the 4.7 GB GENCODE GTF is a one-time cost per
  GTF. Later runs start about 50 s sooner. The machine's load average was
  16.8 at the start, so the times are approximate.

**Golden tests:** unchanged.

