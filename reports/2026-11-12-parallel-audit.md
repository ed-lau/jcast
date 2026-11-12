# Audit: run time on a whole genome, and parallelizing JCAST

| | |
|---|---|
| **Date** | 2026-11-12 |
| **Status** | Audit. Recommendations in §5; decisions requested in §6 and approved 2026-11-12 (§7). No code changed. Measured at `000deb1`. |
| **Question** | A deep RNA-seq run takes days on one thread. Where does the time go, and how much would running junctions in parallel gain now that `--uniprot` removes the API calls and the UniProt cache? |
| **Reproduce** | `python reports/parallel_audit/audit_parallel.py --rmats <dir> --gtf <gtf> --genome <fa> --uniprot <fasta> --write_probe <a large psq_T1.fasta>` (≈ 15 min) |
| **Inputs** | A real whole-genome run (NMDI RNA-seq, not bundled): rMATS 599,644 junctions in 17,068 genes; GENCODE v50 primary-assembly GTF (11.3 M rows, 78,941 genes); GRCh38 primary assembly; Swiss-Prot 2026_04 with isoforms. Mac Studio, 24 CPUs, 64 GB. The user's own JCAST run was using one core throughout. |
| **Evidence** | [parallel_audit/genome/summary.txt](parallel_audit/genome/summary.txt) · [profile.txt](parallel_audit/genome/profile.txt) · [variants.tsv](parallel_audit/genome/variants.tsv) · [parallel.tsv](parallel_audit/genome/parallel.tsv) · chr15 check: [parallel_audit/chr15_debug/](parallel_audit/chr15_debug/summary.txt) |

## 1. Summary

- **95% of the time goes to scanning the whole GTF.**
  - Every GTF lookup filters all 11.3 M rows for the junction's gene,
    about 9 times per junction, at 60 ms each.
  - Reading the 11.3 M rows once into a per-gene index takes 4.6 s.
    With it, each junction sees only its gene's rows.
  - **That alone is 19× faster on one core,** with byte-identical output.
- **The second cost grows during the run.** Before every write, the FASTA
  writer re-reads the whole output file to check for duplicates.
  - At the size `psq_T1.fasta` has now reached in the running job (15.8 MB),
    that check takes 154 ms per write.
  - Once the GTF scan is fixed, a junction takes about 30 ms, so these
    checks would dominate.
  - A set of the sequences already written, kept in memory, gives the same
    result for free.
- **Parallelizing then scales almost linearly,** with identical output.
  - 16 worker processes: 11.8× (454 junctions/s).
  - Total memory 10 GB, against 9.9 GB for one process: forked workers
    share the genome and GTF.
- **Projected for this data set** (599,644 junctions):

  | | Junctions/s | Time |
  |---|---|---|
  | Current code (the running job) | 2.8, slowing as the files grow | ≈ 60 h or more |
  | Per-gene index + in-memory de-duplication, 1 core | ≈ 33–38 | ≈ 5 h |
  | The same, 16 workers | ≈ 450 | ≈ 25 min, including ≈ 1.5 min to load the GTF and genome |

## 2. Why the old multithreading was disabled, and what changed

The `ThreadPoolExecutor` loop (`f84e10c`, July 2020) is still in
`run_jcast`, commented out. It had three problems:
- **Threads cannot run the work in parallel.** The time is spent in
  Python and pandas code that holds the interpreter lock (GIL). Only
  processes run it on several cores.
- **Shared state:** concurrent appends to the same FASTA files, and the
  SQLite UniProt cache.
- **API load and progress bars:** too many concurrent UniProt API calls,
  and no working progress bar.

With `--uniprot` (PR 14), the API calls and the cache are gone. The design
in §4 fixes the other two.

## 3. Where the time goes (current code)

300 junctions sampled at random (seed 1) from the five rMATS files, run
one at a time under cProfile:

| | Seconds | Share |
|---|---|---|
| Total | 201.6 | |
| `DataFrame.query` on the whole GTF (2,761 calls, 60–70 ms each, almost all in `isin` over 11.3 M gene ids) | 191.5 | 95% |
| … of which in `gene_proteins` (the `--uniprot` sequence mapping, once per gene) | 104.3 | |
| … in `trim_cds` and `get_translated_phase` | 92.5 | |
| `write_seqrecord_to_fasta` (outputs still small in this sample) | 0.5 | 0% |

`gene_proteins`, added in PR 14, scans the whole GTF once per gene. In
this sample nearly every junction is in a new gene, so it takes half the
time. In a full run, genes repeat (17,068 genes for 599,644 junctions), and
it adds ≈ 1.7 h to the ≈ 60 h. It becomes negligible with the per-gene
index.

**The write cost**, timed separately: one `write_seqrecord_to_fasta`
call against a file holding the first *n* records of the running job's
`psq_T1.fasta`:

| Records already in the file | File size | Per write |
|---|---|---|
| 1,000 | 0.8 MB | 5 ms |
| 10,000 | 12.2 MB | 61 ms |
| 12,215 | 15.8 MB | 154 ms |

With `-c`, every skipped or failed junction also checks
`psq_canonical.fasta`, which holds 10,048 records (7.7 MB) at this point in the run. The files grow
through the run, and so does the cost.

## 4. Candidate changes, measured

Each was patched in at runtime; jcast code is unchanged. The same 300
junctions as §3:

| Variant | Seconds | Junctions/s | Speedup | Output files | Fates |
|---|---|---|---|---|---|
| Current | 187.2 | 1.6 | 1.0 | – | – |
| Per-gene GTF index | 9.7 | 31.1 | 19.4 | identical | identical |
| In-memory de-duplication | 186.9 | 1.6 | 1.0 | identical | identical |
| Both | 9.0 | 33.2 | 20.7 | identical | identical |

The in-memory de-duplication gains nothing here only because a 300-junction
sample writes small files; §3 shows its cost at full-run sizes.

**The per-gene index** is built once with `groupby('gene_id')` (4.6 s).
Every GTF lookup in JCAST starts by filtering on the junction's gene id,
so the gene's rows give the same answers.

**Parallel prototype.** It uses both changes and a pool of forked worker
processes, each taking chunks of 20 consecutive junctions.
- **Workers:** each translates its chunk and returns, for each junction,
  the fate and the records it would have written. Workers write no
  files.
- **The main process:** receives results in junction order (`imap`), then
  de-duplicates and writes them. The output is therefore byte-identical to
  a serial run, whatever the number of workers.
- **Progress bar:** the main process can drive it, one tick per junction.

3,000 junctions sampled at random (seed 2):

| Workers | Seconds | Junctions/s | Speedup | Identical to serial | Peak memory (all processes) |
|---|---|---|---|---|---|
| 1 | 78.0 | 38.5 | 1.0 | yes | 8.6 GB |
| 4 | 21.4 | 140.1 | 3.6 | yes | 9.0 GB |
| 8 | 11.6 | 257.9 | 6.7 | yes | 9.4 GB |
| 16 | 6.6 | 453.6 | 11.8 | yes | 10.1 GB |

The bundled chr15 data give the same picture at smaller scale: the query
share is 85%, every variant's output is identical, and 4 workers give
3.6×. The index gains only 1.3× there, because the chr15 GTF has 97,275 rows.

**Not measured:** logging. The audit silenced the logger, while a real run
writes about 300 MB of log. Writing it from the main process, in order (§5),
keeps that cost on one core. It is small per line, but it would bound the
speedup at high worker counts.

## 5. Recommendations

1. **Per-gene GTF index** (roadmap PR 9).
   - Build `{gene id: rows}` once after reading the GTF, and give each
     junction its gene's rows.
   - Output unchanged; the golden tests guard it.
   - Do this first: on its own it turns days into hours.
2. **In-memory de-duplication of FASTA writes** (PR 9).
   - Keep the set of sequences written to each file.
   - On start-up, seed it from any existing file, as the current re-read
     effectively does.
   - Output unchanged.
3. **Parallel junctions, with `-n/--threads N` worker processes.**
   - Workers return fates, records and log lines. The main process writes
     all three in junction order, so the output and the log match a
     serial run's.
   - The progress bar runs in the main process, over results as they
     arrive.
   - **Start method:** fork, on macOS and Linux. The workers then share the
     genome and the GTF index without copying them; measured above, memory
     stays flat. Windows has no fork, so it runs serially.
   - **Default N:** 1, so a run never takes every core on a shared machine
     unless asked.
   - **Without `--uniprot`:** parallel runs would hit the API concurrently,
     and the SQLite cache is not safe across processes. Options: allow
     `-n` > 1 only with `--uniprot` or `use_gtf_only` (simplest), or fetch
     the needed canonicals serially first.
4. **Order:** 1 and 2 as one PR (no output change), then 3. Each is
   verified as here: identical files and fates, on chr15 and on this data
   set.

## 6. Decisions requested

1. Proceed with the per-gene index and in-memory de-duplication (§5.1–2) as
   the next PR?
2. Parallel junctions (§5.3):
   - Option name: `-n/--threads` (as the old commented-out option), or
     `-p/--processes`?
   - Default 1 (recommended), or all cores but one?
3. Without `--uniprot`, allow only one worker (recommended), or fetch
   canonicals first?
4. Fork with a serial fallback on Windows (recommended)? The alternative
   is spawn, which works on every OS. Each worker would then receive its
   genes' GTF rows and genome sequence instead of sharing them: more code,
   for a platform JCAST is rarely run on.

## 7. Decisions (2026-11-12)

1. **The per-gene GTF index and in-memory de-duplication go first, as one PR
   (PR 9).** The user stopped the running job, to rerun after this.
2. **Parallel junctions: `-p/--processes N`, default 1** (PR 15).
3. **Without `--uniprot`, only one process.** The CLI rejects `-p` > 1
   there, and the help text and README say why.
4. **Fork on macOS and Linux; serial on Windows.**

