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

### Tiers in v0.4.0 (summary of PRs 4–4c)

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

Cumulative effect on the bundled chr15 data:

| Fate | v0.3.5a | v0.4.0 |
|---|---|---|
| T1 | 185 | 219 |
| T1 ALT TERMINUS | – | 25 |
| T2 | 25 | 30 |
| T3 | 2 | – |
| T4 | 236 | 289 |
| FAIL | 360 | 153 |
| Skipped, low count | 297 | 297 |
| Skipped, outside the coding region | – | 92 |

`psq_T1.fasta` grows from 152 to 173 records. The new
`psq_T1_alt_cterm.fasta` and `psq_T1_alt_nterm.fasta` hold 18 and 4.
`psq_canonical.fasta` is unchanged (268). Re-running is recommended for any analysis that used RI
output, C-terminal events or lower tiers.

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
- With `-c`, write the canonical for every gene, whatever its junctions'
  tiers (REVIEW.md MED-14; 2 chr15 genes affected).
- Output-directory creation guarantee; route the UniProt cache through
  `args.out` (or the platform user-cache dir) instead of CWD.
- Stitch-anchor uniqueness check in `stitch_to_canonical_aa`.
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
