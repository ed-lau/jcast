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

### Added
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
- [README.md](README.md) rewritten to consolidate the install, quick-start,
  output-format, and FAQ content previously split between `README.md`,
  `TODO.md`, and the Quarto sources under `jcast_website/`.

### Fixed
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
- GTFs without `transcript_support_level`: the translation-bound lookup now
  checks for the column instead of relying on exception handling, which
  could pass `None` to `len()` (REVIEW.md MED-7).
- `-m` on pandas ≥ 2: `DataFrame.append` replaced with `pd.concat`.
- `-m` on matplotlib ≥ 3.6: removed the `seaborn-white` style call.
- GTFs without `transcript_support_level` on pandas ≥ 2: catch
  `pandas.errors.UndefinedVariableError` (the old
  `pandas.core.computation.ops` location no longer exists).

### Removed
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

- Treat soft-masked (lowercase) and ambiguous IUPAC bases safely in
  `make_pep`.
- Build the RI retained-intron slice from the retained exon alone; it
  currently contains both flanking exons twice (PR 4a2).
- Take the translation phase from the first retained coding exon (PR 4b).
- Treat annotated alternative stop codons as natural C-termini, written to
  `psq_T1_alt_cterm.fasta` (PR 4c).
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
