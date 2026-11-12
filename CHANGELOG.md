# Changelog

All notable changes to JCAST are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/).

For the upcoming v0.4.0 release plan and detailed code review, see
[REVIEW.md](REVIEW.md).

## [Unreleased] — v0.4.0 hardening release

The v0.4.0 line is a hardening release: dependency modernization, correctness
fixes, scientific-integrity audits, and the introduction of an optional Qt
GUI. Implementation is staged as the patch series described in
[REVIEW.md → §9](REVIEW.md#9-v040-hardening-roadmap-revised).

### Tiers in v0.4.0 (summary of PRs 4–5)

How a junction is translated and tiered has changed substantially. Every
change has an audit under [reports/](reports/README.md).

| | v0.3.5a | v0.4.0 |
|---|---|---|
| Coding-region bounds | C-terminal bound taken from a **start** codon | Most upstream start codon to most distal stop codon (TSL ≤ 2 protein-coding) |
| UTR exons | Kept unless a bound fell inside them | Removed. Events entirely outside the coding region are **skipped** with a new fate |
| RI retained slice | Flanking exons included twice | Retained exon alone |
| Reading frame | One phase from a fixed exon matched exactly. Missing for 26 % of junctions after trimming | Each slice's own frame, read at its 5′ end |
| Frameshift test (T1 vs T2) | Raw slice lengths | Coding lengths, from each slice's frame (PR 4b2) |
| Alternative termini | A stop at an alternative annotated stop codon read as premature (T4/FAIL); slices without the start codon read from mid-CDS | **T1 ALT TERMINUS**: annotated alternative stop and start codons end and start the slice; written to `psq_T1_alt_cterm` / `psq_T1_alt_nterm` / `psq_T1_alt_both.fasta` (PR 4c) |
| No annotated frame | Three-frame guess, keeping the last frame tried | The slice is not translated |
| Tiers | T1, T2, T3, T4 | **T1, T1 ALT TERMINUS, T2, T4.** T3 and `psq_T3*.fasta` removed |
| Joining to the canonical | First occurrence of a 10-residue anchor at each end, even at an end that is already complete | A complete end (stop codon, gene start codon) is not joined; anchors of 7 residues must occur once, in order (PR 5) |

Cumulative effect on the bundled chr15 data:

| Fate | v0.3.5a | v0.4.0, `-r 2` (as v0.3.5a's `-r 1`) | v0.4.0, default `-r 1` |
|---|---|---|---|
| T1 | 185 | 220 | 301 |
| T1 ALT TERMINUS | – | 24 | 26 |
| T2 | 25 | 30 | 49 |
| T3 | 2 | – | – |
| T4 | 236 | 289 | 403 |
| FAIL | 360 | 153 | 207 |
| Skipped, low count | 297 | 297 | 0 |
| Skipped, outside the coding region | – | 92 | 119 |

Changes to the output files:
- **`psq_T1.fasta`:** 152 records in v0.3.5a; 179 in v0.4.0 at `-r 2`;
  244 at the new default `-r 1` (PR 6).
- **New files:** `psq_T1_alt_cterm.fasta` and `psq_T1_alt_nterm.fasta`
  hold 18 and 5 at `-r 2`.
- **`psq_canonical.fasta`:** 268 records up to PR 5; 269 from PR 8, with
  one record for every gene that has a canonical.

Re-running is recommended for any analysis that used RI output,
C-terminal events or lower tiers.

**Decisions recorded during PR 4b.** The candidate carried two switches for
the audit only, `three_frame_fallback` and `frame_for_start_less_slice`.
Neither was released.

- **The three-frame fallback is removed.** Once frames were read correctly,
  it changed only 15 junctions, and no annotated sequence depended on it.
- **Slices lacking the start codon keep their annotated frame.** Withholding
  it removed real isoforms. An annotated-start rule is planned instead
  (PR 4c).

See
[reports/2026-11-10-phase-audit.md](reports/2026-11-10-phase-audit.md) §7
and §9.

### Added
- **`-p/--processes N`: translate junctions in N processes (PR 15).**
  - **Same output:** the output files and the log are the same as with one
    process. Workers only translate; the main process writes everything in
    rMATS file order.
  - **Measured on a whole-genome sample:** 16 processes give 736 junctions/s,
    11.3× one process, so a 600,000-junction run takes about 15 minutes.
  - **Default 1.**
  - **Requires `--uniprot`:** without it, canonicals come from the UniProt
    API and its local cache, which cannot be shared between processes.
  - **macOS and Linux only:** Windows runs with one process. The commented-out
    2020 thread-pool code is removed.

  See
  [reports/2026-11-12-parallel-audit.md](reports/2026-11-12-parallel-audit.md).
- **`--uniprot FASTA`: canonicals from a supplied UniProt FASTA instead of the
  UniProt API (PR 14).** Recommended mode: offline, and the file records
  which UniProt release was used.
  - **Choosing the canonical:** each gene's canonical is the UniProt
    canonical whose sequence equals one of the gene's GTF-translated
    proteins, else the one whose gene name (`GN=`) is the gene symbol.
    Either must be unique, and reviewed entries come before unreviewed
    ones. A gene with neither uses the GTF canonical.
  - **New output files:**
    - `uniprot_canonicals.tsv`: the canonical chosen per gene, and how.
    - `psq_isoforms.tsv`: JCAST records identical to a UniProt canonical
      or isoform.
    - `jcast_uniprot_merged.fasta`: the supplied FASTA, plus the new
      `psq_T1` and `psq_T1_alt_*` records joined to a Swiss-Prot canonical,
      renamed `sp|<accession>-J<n>|<entry name>`. This replaces the post-run
      R merge; on chr15 it gives the same records.
  - **Unchanged:** without `--uniprot`, the API path is unchanged, and
    the `psq_*.fasta` files keep their format.
  - **On chr15:** with Swiss-Prot 2026_04, the canonical equals the cached
    API answer for 183 of 184 genes; no junction changes fate.

  See
  [reports/2026-11-11-uniprot-fasta-evaluation.md](reports/2026-11-11-uniprot-fasta-evaluation.md).
- **Alternative C- and N-termini (affects output).**
  - A slice whose first stop codon is an annotated protein-coding stop
    codon is complete at that end, not a premature stop. It is joined to
    the canonical at its N-terminus only.
  - A slice that lacks the gene's start codon (present in the other slice)
    but contains another annotated start codon is translated from that start
    and joined at its C-terminus only.
  - Such junctions get the new fate `SUCCESS TIER 1 ALT TERMINUS`, and their
    slices go to `psq_T1_alt_cterm.fasta`, `psq_T1_alt_nterm.fasta` or
    `psq_T1_alt_both.fasta` (with `_orphan` files) by which terminus differs
    from the canonical.
  - Annotated sites use the same TSL threshold as the coding bounds.

  On chr15, 26 junctions change fate: 16 come up from FAIL or T4, and 9 were
  T1 whose start-less slice recreated the canonical and now gives the
  annotated isoform. 22 alt records are written, 19 of them exact
  annotated proteins. No `psq_T1.fasta` record moves, and no annotated
  sequence is lost. See
  [reports/2026-11-10-alt-termini-audit.md](reports/2026-11-10-alt-termini-audit.md).
- Optional `jcast/gui/` PySide6 skeleton; launch via `python -m jcast --gui`
  or `jcast-gui`. Requires `pip install jcast[gui]`.
- [REVIEW.md](REVIEW.md): comprehensive critical review covering architecture,
  correctness bugs, scientific integrity, dependencies, and roadmap.
- [CHANGELOG.md](CHANGELOG.md) (this file), replacing the older `CHANGES.txt`.

### Changed
- **Read-count cutoffs mean what they say (PR 6).** A junction is kept if
  its read count is **at least** the cutoff; it used to need one read more.
  - **`-r N`** keeps junctions with N reads. **The default `-r 1` now keeps
    every junction:** rMATS reports no junction with 0 reads. Before, it
    dropped those with 1 read.
    - On chr15, 297 more junctions are translated and `psq_T1` grows from
      179 to 244 records.
    - For the previous default behaviour, pass `-r 2`. An explicit
      `-r N` now keeps one read lower than before.
  - **`-m` cuts at the model's own boundary.** It used to cut two reads
    above it: once because the model is fitted to ln(count + 1), once from
    the comparison.
    - The cutoff is logged, and drawn in `model.pdf`, in reads.
    - On the 600,000-junction whole-genome test run, 3,410 more
      junctions pass (cutoff 40 reads, was effectively 42).
  - **The count is unchanged but now named for what it is:** the smaller of
    the summed inclusion and skipping junction counts (`Junction.min_jc`,
    was `sum_sjc`).
    - Help text, README header field 11 and the plot said "skipped
      junction count".
    - The FASTA header's `r<count>` value is unchanged.
  - The effective cutoff is logged at the start of every run.

  See [reports/2026-11-12-thresholds-runs-canonical-audit.md](reports/2026-11-12-thresholds-runs-canonical-audit.md) §2 and §9.
- **Filters first (PR 6).** The read-count and FDR filters run before
  trimming, translation and the canonical lookup. A skipped junction only
  writes its gene's canonical with `-c`.
  - Before, every junction was translated and logged first. On a
    whole-genome `-m` run, 79% of junctions were skipped, after taking
    about 78% of the run time.
  - No output file changes. A junction that is both below the cutoff and
    masked now reports low coverage.
  - Without `-c`, `uniprot_canonicals.tsv` lists only the genes of
    translated junctions.

  See [reports/2026-11-12-log-audit.md](reports/2026-11-12-log-audit.md).
- **Each run writes a fresh output directory (PR 7).**
  - `-o` is created if it does not exist. Before, it had to exist, even the
    default `out/`.
  - **A directory that already holds JCAST output is refused** unless
    `--overwrite` is given (new option). `--overwrite` removes only JCAST's
    own files: `logfile.log`, `psq_*.fasta`, `model.pdf` and the
    `--uniprot` outputs.
  - Before, output was appended. Sequences already in a file were skipped,
    so a rerun into the same directory kept the first run's records: a
    `-r 30` run into a `-r 1` directory left the `-r 1` files byte for byte.
- **One failing junction no longer stops the run (PR 7).**
  - An exception while translating a junction gives it the fate
    `FAILURE (exception). <type>: <message>`, and the run continues.
  - The junction writes nothing. Its traceback goes to the log.
  - The first failure is shown on the console. At the end, a warning
    counts the failed junctions. The exit status stays 0.
  - The same holds with `-p`, where the failure now costs only that
    junction instead of its batch of 20.
- **Caches moved to the user cache directory (PR 7).**
  - The location is `~/Library/Caches/jcast` on macOS or `~/.cache/jcast`
    on Linux, or `--cache_dir` (new option).
  - New dependency: `platformdirs`.
  - **GTF cache:**
    - It used to be `<gtf>.cached` next to the GTF. That failed in
      read-only directories, was used even after the GTF changed, and
      could be read half-written by a second run.
    - It is now keyed by the GTF's path, size and modification time, and
      renamed into place once complete.
    - It is stored as a pickle of the same table. Reading the GENCODE v50
      cache takes 8.5 s instead of 49 s.
    - **The first run on each GTF builds it again;** old `.cached` files can
      be deleted.
  - **UniProt API cache:** it was `cache/uniprot-cache.db` under the
    working directory. To keep earlier answers, move that file into the
    cache directory.
  - If a cache cannot be written, JCAST warns and runs without it.

  See [reports/2026-11-12-thresholds-runs-canonical-audit.md](reports/2026-11-12-thresholds-runs-canonical-audit.md) §3 and §8.
- **About 25× faster on a whole genome, with identical output (PR 9).**
  - **GTF lookups:** the GTF is indexed by gene once, instead of being
    filtered whole (11 M rows for GENCODE) several times per junction.
  - **FASTA writes:** written sequences are kept in memory, instead of the
    output file being re-read before every write.
  - **Measured:** on 2,000 junctions of a whole-genome run, 1.9 → 46.6
    junctions/s, with every output file and log message identical. See
    [reports/2026-11-12-parallel-audit.md](reports/2026-11-12-parallel-audit.md).
- **Read-count cutoff model (`-m`) no longer uses `pomegranate`.** pomegranate
  1.x removed the API JCAST used (`ImportError: cannot import name
  'GeneralMixtureModel'`), and 0.x no longer builds on current Python. The same
  Gamma + Normal mixture on ln(count + 1) is now fitted by a deterministic EM in
  numpy/scipy. On the bundled datasets, cutoffs agree with pomegranate 1.x within
  1–2 reads wherever 1.x converges; see
  [REVIEW.md → §2](REVIEW.md#2-the-pomegranate-break-and-the-immediate-fix).
  When no stable two-component fit exists, JCAST logs a warning and falls back
  to `-r` instead of crashing. Fitted parameters are now written to the log.
- scikit-learn and matplotlib are imported lazily, so a missing optional
  package no longer breaks every CLI invocation.
- Dependency bounds aligned between `setup.py` and `requirements.txt`:
  `gtfparse>=1.3,<2` (3.x is polars-based and incompatible), `pandas>=1.5,<4`,
  and `numpy` / `scipy` declared explicitly.
- **Anchor length for joining to the canonical: 10 → 7 residues**
  (`params.aa_stitch_length`). 7 is where repeated anchors and chance
  matches level off across the MANE Select proteome. On chr15, joins
  missed for want of a long-enough shared end fall from 17 to 3.
  Alternative-termini file routing still compares 10 residues (new
  `params.aa_terminus_length`).
- **`--mask` is deprecated and ignored**, and logs a warning. Every base
  outside `ACGT` is now skipped, which covers its old default (`NnXx`). It
  will be removed in the next release.
- The fate for skipped slices reads `masked or ambiguous bases` (was
  `masked regions`).
- The README no longer recommends hard-masked genomes. On chr15 they
  silently skip 19% of junctions.
- [README.md](README.md) rewritten to consolidate the install, quick-start,
  output-format, and FAQ content previously split between `README.md`,
  `TODO.md`, and the Quarto sources under `jcast_website/`.

### Fixed
- **`-c` now writes every gene's canonical (PR 8).**
  - Before, the T1, T2 and T4 paths wrote it only when a translated slice
    equalled it. Masked junctions never wrote it: the masked-base check
    came before the canonical lookup.
  - Missing before: 1 chr15 gene at `-r 2` and 4 at `-r 1` (with
    `--uniprot`), and 5 of
    13,622 genes in a whole-genome run (REVIEW.md MED-14).
- **`write_canonical` no longer asks the UniProt API again** for an empty
  canonical (PR 8).
  - It sent two requests and discarded the answers, so it never wrote
    anything; it also did so in GTF-only mode (MED-9).
  - On chr15 this was 244 of 1,008 requests.
- **UniProt API failures no longer abort the run or silently change it
  (PR 8; MED-10).**
  - A connection error, timeout or server error after the retries now
    fails that junction (`FAILURE (exception)`, PR 7). Before, it aborted
    the run.
  - So does any refusal, such as HTTP 429 (too many requests). Before, a
    refusal was taken as "no canonical", which moved all the gene's
    records to orphan files: chr15 `psq_T1` 149 → 0 when every request
    was refused.
  - Requests now time out after 30 s (`params.uniprot_timeout`); they
    could hang forever. 429 is retried like server errors.
  - The end-of-run warning counts API failures and suggests `--uniprot`.

  See [reports/2026-11-12-thresholds-runs-canonical-audit.md](reports/2026-11-12-thresholds-runs-canonical-audit.md) §4 and §10.
- **`--uniprot` matches a gene name only if no other gene in the GTF has
  it (PR 8).**
  - GENCODE v50 has 484 names shared by several genes. For example, two
    genes are named PINX1, and the second took PINX1's Swiss-Prot entry by
    name.
  - In a whole-genome run, 77 genes lose a by-name canonical. 72 are
    duplicated lncRNA or pseudogene loci; 5 are coding loci named like an
    established gene, which keeps the entry by sequence.
  - 2 T4 records become orphans. 49 Swiss-Prot entries of duplicated
    pseudogene names leave `psq_canonical.fasta`; they stay in the merged
    database. See [reports/2026-11-12-uniprot-mapping-audit.md](reports/2026-11-12-uniprot-mapping-audit.md).
- **A second run in the same process logged every line twice,** once
  into the first run's log (REVIEW.md MED-13). Matters for the GUI,
  notebooks and tests. Each run now closes its log when it ends (PR 7).
- **Joining slices to the canonical (affects output; REVIEW.md MED-3).**
  On chr15, 12 of 998 joins were wrong, and 10 of them gave exactly the
  canonical, so the isoform was lost:
  - A Tier 4 peptide ending at its premature stop codon was joined at its
    C-terminus too. It is now joined at its N-terminus only.
  - A peptide starting at the gene's start codon was joined at its
    N-terminus. It is now joined at its C-terminus only.
  - An anchor occurring more than once took its first occurrence. An anchor
    is now used only if it occurs once, and the start anchor must come
    before the end anchor; otherwise the slice is written as an orphan.

  Verified against genomic coordinates on chr15: 0 wrong joins (was 12).
  One junction changes fate (GTF2A2 SE 11945, T1 ALT TERMINUS → T1). Its
  7-residue peptide can now be joined and gives the canonical. See
  [reports/2026-11-11-stitch-anchor-audit.md](reports/2026-11-11-stitch-anchor-audit.md).
- **Soft-masked genomes aborted the run** (REVIEW.md HIGH-2/HIGH-3).
  Lowercase bases raised `KeyError` in translation, and nothing caught it. On
  chr15 with Ensembl's soft-masked genome, 147 junctions raised in a default
  run, and the run stopped at the third junction. Genome sequence is now
  uppercased when read, and gives exactly the output of the unmasked genome.
- **Ambiguous bases no longer raise.** A slice containing any base outside
  `ACGT` (N or another IUPAC code) is skipped, as N already was. Elsewhere,
  e.g. the canonical translated from the GTF, a codon with an ambiguous base
  makes the translation empty rather than raising. `get_complementary` now
  uses Biopython's IUPAC-aware reverse complement; its fallback was
  unreachable. See
  [reports/2026-11-11-masked-bases-audit.md](reports/2026-11-11-masked-bases-audit.md).
- **Gzipped genomes (`.fa.gz`) crashed every run** with `UnicodeDecodeError`.
  The GTF and genome arguments are now plain paths (no `argparse.FileType`),
  and missing files are reported by argparse.
- Test suite repaired and moved to pytest fixtures. Adds a chr15
  end-to-end golden test pinning per-tier outcome counts and output FASTA
  contents, plus model cutoff tests. Tests no longer write caches into the
  repository.
- **Translation-end bound taken from start codons (affects output).** The
  bound used to trim splice events to the coding region was computed from
  start codons instead of stop codons (since 2020-07). As a result:
  - splice events spanning a protein's C-terminal exon were almost never
    translated;
  - some events containing a start-codon exon yielded 5′-UTR-derived or
    event-less sequences in lower tiers.

  The bound now uses annotated stop codons, and each side falls back to no
  trimming when it cannot be determined. On the bundled chr15 data, 48 of
  1,105 junctions change tier. `psq_T1.fasta` gains 12 sequences and loses
  none, and `psq_canonical.fasta` is unchanged. Re-running is recommended to
  recover C-terminal isoforms. See
  [reports/2026-11-03-tx1-translation-bounds-audit.md](reports/2026-11-03-tx1-translation-bounds-audit.md).
- **Splice events are now trimmed to the coding region exon by exon (affects
  output).**
  - Exons lying entirely in the UTR are removed. Previously 172 of 1,105
    chr15 junctions kept one, and some 5′-UTR translations were reported as
    successes.
  - A5SS/A3SS events on the `−` strand are trimmed at the correct bound.
  - The retained region of RI events is trimmed.
  - Events lying entirely outside the coding region are skipped with the new
    fate `SKIPPED. Splice event lies outside the annotated coding region.`;
    the canonical is still written with `-c`.

  On chr15, 103 junctions change fate, 90 of them to the new skip.
  `psq_T1.fasta` gains 5 sequences and loses none, and no sequence matching
  an annotated protein is lost. See
  [reports/2026-11-10-utr-exon-trimming-audit.md](reports/2026-11-10-utr-exon-trimming-audit.md).
- **Retained-intron (RI) slices contained both flanking exons twice (affects
  output).** rMATS's retained exon already spans the upstream and downstream
  exons, but JCAST added them again. So every intron-retaining RI sequence
  was translated from a duplicated sequence, in every release since at least
  v0.2.6.

  The retained slice is now the retained exon alone. On chr15, 5 of 88 RI
  junctions change fate, and 58 RI records, mostly in `psq_T4.fasta`, are
  replaced by correct translations. `psq_T1.fasta` and
  `psq_canonical.fasta` are unchanged. Re-running is recommended if RI
  output was used. See
  [reports/2026-11-10-ri-slice-audit.md](reports/2026-11-10-ri-slice-audit.md).
- **Each slice is translated in its own annotated reading frame (affects
  output).** The frame used to come from one fixed exon that had to match a
  CDS row exactly, and was shared by both slices. After trimming, that
  lookup failed for 184 of 716 non-skipped chr15 junctions.
  - The frame is now read at each slice's 5′ end from the protein-coding CDS
    rows covering it, so slices starting in different exons get their own
    frames.
  - Where the old lookup found a frame (532 junctions), the new one agrees
    every time.
  - FASTA headers report each slice's own frame.

  On chr15, 88 junctions change fate, and FAILURE falls from 238 to 157.
  No sequence matching an annotated protein is lost, and T1/T2 gain the
  annotated GTF2A2-208 and LINS1-207. `psq_canonical.fasta` is unchanged.
  See
  [reports/2026-11-10-phase-audit.md](reports/2026-11-10-phase-audit.md).
- **The T1/T2 frameshift test ignored per-slice frames (affects output;
  regression from the per-slice frame change).** It compared raw slice
  lengths, which is wrong when the two slices start in different frames.
  It now compares coding lengths.

  On chr15, 17 junctions are reclassified: 16 from T2 to T1, including the
  annotated LINS1-207, and 1 from T1 to T2. No sequence is lost or created;
  records only move between the T1 and T2 files. See
  [reports/2026-11-10-frameshift-audit.md](reports/2026-11-10-frameshift-audit.md).
- GTFs without `transcript_support_level`: the translation-bound lookup now
  checks for the column instead of relying on exception handling, which
  could pass `None` to `len()` (REVIEW.md MED-7).
- **rMATS rows without a gene symbol** (`NA`) are named
  `NoGeneName_<gene id>` in FASTA headers, as v0.3.5a intended. pandas reads
  `NA` as `NaN`, so the v0.3.5a check (`is None`) never fired. Headers carried
  `nan` and `str(junction)` raised `TypeError`.
- `-m` on pandas ≥ 2: `DataFrame.append` replaced with `pd.concat`.
- `-m` on matplotlib ≥ 3.6: removed the `seaborn-white` style call.
- GTFs without `transcript_support_level` on pandas ≥ 2: catch
  `pandas.errors.UndefinedVariableError` (the old
  `pandas.core.computation.ops` location no longer exists).

### Removed
- **Tier 3 and the three-frame fallback** (`psq_T3.fasta`,
  `psq_T3_orphan.fasta`, `fates.TIER3`, and `Sequence.translate(use_phase=...)`).
  They guessed a reading frame when the annotated one was missing or failed,
  keeping the last frame tried (REVIEW.md MED-2, now closed). A slice without
  an annotated frame is no longer translated.
- The retry with the gene's second start codon in `trim_cds`. Once UTR exons
  are removed, it deleted exons coding in well-supported transcripts.
- `pomegranate` dependency.
- `CHANGES.txt` — content migrated to this file.
- `TODO.md` — items (gff3 input, APPRIS-tag canonical, UniProt-failure
  fallback) absorbed into [REVIEW.md → §9](REVIEW.md#9-v040-hardening-roadmap-revised).

### Planned (not yet implemented in this milestone)

The following items are scheduled for v0.4.0 but live in subsequent focused
PRs. See [REVIEW.md → §9](REVIEW.md#9-v040-hardening-roadmap-revised) for the full
table.

- Catch an exception in one junction, record it as that junction's fate,
  and continue, instead of aborting the whole run (PR 7).
- Remove the deprecated `--mask` option (next release).
- An option to add T2 and T4 records to `jcast_uniprot_merged.fasta`.
- `--uniprot` (and, from a subprocess or with another start method, `-p`)
  in the GUI.
- A helper that downloads a UniProt proteome and names the file by its
  release (`X-UniProt-Release`).
- With `-c`, write the canonical for every gene, whatever its junctions'
  tiers (REVIEW.md MED-14; 2 chr15 genes affected).
- Output-directory creation guarantee; route the UniProt cache through
  `args.out` (or the platform user-cache dir) instead of CWD.
- Bump Python floor to 3.10.
- GitHub Actions CI matrix across 3.10–3.13.
- Drop the dev-only `--g_or_ln` flag.

## [0.3.5a] — 2022-10-28

### Added
- Skip a junction if its nucleotide sequence contains masked bases.
- `--mask` argument to configure masked-nucleotide letters (default: `NnXx`).
- `-s / --splice_type` argument to translate only a subset of splice types
  (MXE, RI, SE, A3SS, A5SS).

## [0.3.5]

### Fixed
- A small number of slices were not translated.
- JCAST quit if the `--model` flag was not set.

## [0.3.3]

### Fixed
- Bug when the model flag is not set.

## [0.3.2]

### Changed
- `-q / --qvalue` now works by **including** splice junctions within the
  specified range. To preserve the prior behavior, use `-q 0.01 1`. The
  default `-q 0 1` disables filtering.
- `-c / --canonical` now defaults to `False`. Set the flag explicitly to
  translate canonical proteins detected in the RNA-seq experiments even when
  isoform translation fails.
- `-c / --canonical` now writes a canonical sequence even if a splice event
  has low skipped-junction counts or falls outside the specified FDR range.
- `-m / --model` automatically models the likely translatable junction-count
  threshold. Previously the cutoff was supplied manually via `-r / --read`.
- The skipped-junction counts used for modeling now consider the **sum** of
  all biological and technical replicates rather than the minimum.

## [0.3.1]

### Fixed
- Reading Gencode GTF files that lack `transcript_biotype` annotations.
- Retrieval of UniProt canonical sequences for Ensembl transcript versions.
- Compatibility with biopython ≥ 1.78.

## [0.2.8]

### Added
- Option to use the longest coding sequence from the GTF as the canonical
  protein.

### Fixed
- `jcast` being called twice from the direct CLI entry point.

## [0.2.6]

### Fixed
- Some A5SS minus-strand transcripts retrieved incorrect phases from the GTF.
- Some transcripts using a secondary translation start site were not trimmed.
