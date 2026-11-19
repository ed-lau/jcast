# JCAST audit reports

Scientific-integrity audits behind the v0.4.0 hardening series. Each
output-changing PR has an audit, written before or alongside the code
change, with:

- the evidence files and the script that regenerates them;
- regression tests that fail before the change;
- the commit at which the evidence reproduces exactly.

The roadmap is [REVIEW.md §9](../REVIEW.md#9-v040-hardening-roadmap-revised).
User-facing changes are in [CHANGELOG.md](../CHANGELOG.md).

Most audits use the bundled chr15 test set (1,105 rMATS junctions, Ensembl
GRCh38.89) in GTF-only mode with UniProt stubbed out, so they need no
network. There are two exceptions:
- **The UniProt FASTA evaluation** uses a Swiss-Prot release and the local
  UniProt API cache.
- **The parallel audit** uses a whole-genome run (GENCODE v50, GRCh38,
  about 600,000 junctions). Those inputs are not bundled; its scripts take
  the paths as arguments.

The packaging check (PR 10) changes no output. It runs the test suite on
each Python version instead.

The scope of v0.5.0 (the GUI and the command line on Typer, PRs 18, 19 and
13) is here too. It changes no output; its evidence is Typer's behaviour
and the cost of a model preview.

## Audits

| Report | PR | Status | Evidence reproduces at |
|---|---|---|---|
| [tx0/tx1 translation bounds](2026-11-03-tx1-translation-bounds-audit.md) | 4 | Implemented in `da72a0a` | `a555446` (`python reports/tx1_audit/audit_tx1.py`) |
| [Trimming to the coding region (UTR exons)](2026-11-10-utr-exon-trimming-audit.md) | 4a | Implemented in `89e1a82` | `89e1a82` (`python reports/utr_audit/audit_utr.py`) |
| [RI retained-intron slices](2026-11-10-ri-slice-audit.md) | 4a2 | Implemented in `bc7d894` | `bc7d894` (`python reports/ri_audit/audit_ri.py`) |
| [Reading frame and three-frame fallback](2026-11-10-phase-audit.md) | 4b | Candidate `9e27815`; final `2679edb` | four variants at `9e27815`; final check with `--final` |
| [Frameshift test from coding lengths](2026-11-10-frameshift-audit.md) | 4b2 | Implemented in `a0f3338` | `a0f3338` (`python reports/frameshift_audit/audit_frameshift.py`) |
| [Alternative C-termini (first evaluation)](2026-11-10-alt-cterm-tiering-audit.md) | 4c | Superseded by the next row. Measured on PR 4 alone. | `a555446` (`python reports/cterm_audit/audit_cterm.py`) |
| [Alternative C- and N-termini on the current pipeline](2026-11-10-alt-termini-audit.md) | 4c (and 4b2) | **Implemented in PR 4c** (§10), after PR 4b2, which fixes the regression it found (§6). | evaluation: `919f6c7` (`audit_termini.py`); expected after 4b2: `a0f3338` (`audit_termini.py --out …/expected_after_4b2`); verification: `e11cc86` (`verify_4c.py`) |
| [Soft-masked and ambiguous genome bases](2026-11-11-masked-bases-audit.md) | 3 | **Implemented in `fcdc0ee`** (§8); decisions approved 2026-11-11 (§7) | audit: `819cc53`; verification: `fcdc0ee` (`python reports/mask_audit/audit_mask.py --softmasked …dna_sm.chromosome.15.fa.gz [--verify]`); genome survey with `survey_genome.py` (GENCODE GRCh38 + v50 GTF, not bundled) |
| [Joining slices to the canonical; anchor length](2026-11-11-stitch-anchor-audit.md) | 5 | **Implemented in `2142711`** (§9); decisions approved 2026-11-11 (§8) | audit: `887a98a` (`python reports/stitch_audit/audit_stitch.py`); verification: `2142711` (`verify_5.py`); proteome survey with `survey_anchors.py` (GENCODE v50, not bundled) |
| [Run time on a whole genome; parallel junctions](2026-11-12-parallel-audit.md) | 9, 15 | **PR 9 implemented in `cde553d`** (§8); **PR 15 implemented in `9d2865e`** (§9); decisions approved 2026-11-12 (§7) | audit: `000deb1` (`python reports/parallel_audit/audit_parallel.py --rmats … --gtf … --genome … --uniprot …`); PR 9 verification: `cde553d` (`verify_runs.py`, §8); PR 15 verification: `9d2865e` (`verify_runs.py … --run p16=.:16`, §9); whole-genome inputs not bundled |
| [Read-count thresholds; robustness of real runs; canonicals under `-c`](2026-11-12-thresholds-runs-canonical-audit.md) | 6, 7, 8 | Decisions taken 2026-11-12 (§7). **PR 7 implemented in `460005d`** (§8); **PR 6 in `05a6931`** (§9); **PR 8 in `890d631`** (§10) | audit: `8b845aa` (`python reports/runs_audit/audit_runs.py [--rmats … --jcast_out … --uniprot …]`; `bench_gtf_cache.py --cached … --scratch …`); PR 7 verification: `460005d` (`verify_7.py --before … --after …` and `reports/parallel_audit/verify_runs.py`, §8); PR 6 verification: `05a6931` (`verify_6.py --pair …`, §9); PR 8 verification: `890d631` (`verify_8.py --before … --after … [--genome_before … --genome_after … --uniprot_full …]`, §10); whole-genome inputs not bundled |
| [The run log's size and content](2026-11-12-log-audit.md) | 6, 16 | Decisions approved 2026-11-12 (§6): filter first, **done in PR 6 (`05a6931`)**; `junctions.tsv`, a small log and `--verbose`, **done in PR 16 (`a85f8fc`)** (§7) | audit: `8b845aa` (`python reports/log_audit/audit_log.py --log …`; `time_filter_first.py`); PR 16 verification: `a85f8fc` (`verify_16.py --before … --after … [--rmats … --gtf … --genome … --uniprot …]`, §7); whole-genome inputs not bundled |
| [The GTF canonical from GENCODE's `Ensembl_canonical` tag](2026-11-12-ensembl-canonical-audit.md) | 17 | Decisions approved 2026-11-12 (§9): the fallback without readthrough transcripts. **Implemented in PR 17 (`8ddf8c3`)** (§10). Prototype measured on the whole genome (§6) | audit: `ae3d937` (`python reports/canonical_audit/audit_canonical.py genes --gtf … --run … --rmats …`; `… compare --before … --after … --label …`; prototype in `prototype.diff`); PR 17 verification: `8ddf8c3` (`verify_17.py --before … --after … [--genome_before … --genome_prototype … --genome_after …]`, §10); whole-genome inputs not bundled |
| [How `--uniprot` maps genes to UniProt entries](2026-11-12-uniprot-mapping-audit.md) | 8, 17 | Decisions approved 2026-11-12 (§7): mapping kept; unique-name rule and TrEMBL note **done in PR 8 (`890d631`)** (its measured effect: §8); `Ensembl_canonical` fallback **done in PR 17 (`8ddf8c3`)**, without readthrough transcripts ([its audit](2026-11-12-ensembl-canonical-audit.md)) | audit: `8b845aa` (`python reports/mapping_audit/audit_mapping.py --jcast_out … --uniprot … --gtf_cached … --xrefs … --gtf … --genome …`); whole-genome inputs and UniProt cross-references (release 2026_04) not bundled |
| [A supplied UniProt FASTA; joining by coordinates](2026-11-11-uniprot-fasta-evaluation.md) | 14 | **Implemented in `0f3b409`** (§8); decisions approved 2026-11-11 (§7) | evaluation: `8b0c024` (`python reports/uniprot_eval/evaluate_uniprot_fasta.py --uniprot <Swiss-Prot FASTA>`); verification: `0f3b409` (`verify_14.py --uniprot <Swiss-Prot FASTA>`); both need the local UniProt cache |
| [Packaging: Python versions and dependency bounds](2026-11-13-packaging-check.md) | 10, 11 | Decisions taken 2026-11-13 (§5). **Implemented in PR 10 (`4a039a9`)** (§2); **CI in PR 11 (`dee41b6`)** (§6); no output change | PR 10: `4a039a9` (`python reports/packaging/check_install.py --uv … --python … [--resolution lowest]`); needs uv and an interpreter per Python version |
| [Command-line checks: the rMATS files read, and `-m`](2026-11-13-command-line-checks.md) | 12 | Decisions taken 2026-11-13 (§3). **Implemented in PR 12 (`d84f49d`)** (§4); no output change for any run that worked before (§5) | PR 12: `d84f49d` (`python reports/cli_checks/verify_12.py --before … --after …`) |
| [GUI and command-line scope for v0.5.0](2026-11-13-gui-cli-scope.md) | 18, 19, 13 | Decisions taken 2026-11-13 (§7). **PR 18 implemented in `5479ee3`** (§9); **PR 19 implemented in `ec1068a`** (§10); no output change. **PR 20 in `2f0e54c`** (§11): `-m` without a fit falls back to `-r`. **PR 13 implemented in `12e0576`** (§12): the GUI; no output change. **PR 21 in `4cc72c2`** (§13): the GUI after its first run on a Mac | scope: `f37faf9` (`python reports/gui_scope/scope_evidence.py --work … [--python …] [--rmats …]`), needs PyPI, whole-genome rMATS folder not bundled; PR 18 verification: `5479ee3` (`python reports/cli_checks/verify_18.py --before … --after …`); PR 19 verification: `ec1068a` (`python reports/gui_scope/verify_19.py --before … --after … [--rmats …]`); PR 20 verification: `2f0e54c` (`python reports/gui_scope/verify_20.py --before … --after … [--earlier …]`); PR 13 verification: `12e0576` (`python reports/gui_scope/verify_13.py --before … --after … [--rmats …] [--output …]`), needs PySide6; PR 21 verification: `4cc72c2` (`python reports/gui_scope/verify_21.py --before … --after … [--build-python …]`), needs PySide6 and build |
| [The `-m` read-count model](2026-11-15-read-count-model-audit.md) | 23 | Measured on three data sets; §13 and §15.7 answered 2026-11-18; §15.7 done in PR 24 (§17). The known isoforms a cutoff loses on the human cells with rMATS 4.1.0 (§3) and 4.4.0 (§4), on mouse ventricles (§5), and genes against PaxDb's heart proteins (§6); a cutoff on novel events only (§§3.6, 4.5, 5.6); the model's fit (§8: a depth-scaled share, no boundary); two cutoffs, annotated and novel (§9); NMDI14 (§10, inconclusive); a third data set, iPSC differentiation (§14: no conclusion changes); the events rMATS's statistics drop before JCAST reads them, and how far `--statoff` moves `-m` (§15); the search against spectra and the Model page (§16); PR 24, folders made with `--statoff` or from one sample (§17) | §§3–7: `f368783`. A run at `-r 1` with [`audit_records.diff`](model_audit/audit_records.diff) and `JCAST_AUDIT_RECORDS` set, then `python reports/model_audit/cutoff_recovery.py --records … --run … --rmats … --uniprot … --gencode … --gtf … [--salmon …] [--paxdb …] [--check <dir>:<cutoff> …] --out …`; §4.1: `python reports/model_audit/rmats_versions.py --old … --new … --out …`; §8: `python reports/model_audit/model_fit.py --rmats … --out …` (or `--bundled tests/data/model`, offline); §9: `python reports/model_audit/two_cutoffs.py` with `cutoff_recovery.py`'s inputs; §10: `python reports/model_audit/nmd_response.py --rmats … --run … --targets … --dte … --gtf … --out …`; §15: rMATS's post step with `--statoff`, then `python reports/model_audit/dropped_events.py --rmats … --statoff … --dropped … --out …`, a `-r 1` run on the dropped folder, and `dropped_events.py` with `--records … --dropped_records … --dropped_run …`; §16.1: [`novel_only.diff`](model_audit/novel_only.diff); §17: `python reports/model_audit/verify_24.py --before … --after … [--folder <label>=<dir> …] [--union …]`; whole-genome inputs, GENCODE translations and PaxDb not bundled |

Each script compares a baseline commit with the working tree, or patches
the candidate in at runtime, so later code changes its output. To reproduce
a report's evidence files exactly, check out the commit in the last column
first.

## Decisions log

| Date | Decision | Evidence |
|---|---|---|
| 2026-11-10 | Drop exons lying entirely in the UTR; skip events entirely outside the coding region (new fate) | [4a audit](2026-11-10-utr-exon-trimming-audit.md) |
| 2026-11-10 | Defer one coding region per gene vs per transcript | [tx1 audit §8](2026-11-03-tx1-translation-bounds-audit.md#8-residual-issues-out-of-scope-for-pr-4-present-under-both-versions) item 2 |
| 2026-11-10 | Second-start-site retry: deferred, then **removed in PR 4a**, because with UTR-exon removal it deleted exons coding in TSL 1–2 transcripts (20 junctions) | [4a audit §4.3](2026-11-10-utr-exon-trimming-audit.md#43-why-the-retry-had-to-be-removed) |
| 2026-11-10 | Read each slice's frame at its 5′ end from the first retained coding exon; translate in one frame only | [4b audit](2026-11-10-phase-audit.md) |
| 2026-11-10 | **Remove the three-frame fallback and Tier 3** (switch `three_frame_fallback`, audit only, never released) | [4b audit §5.3, §7](2026-11-10-phase-audit.md#53-effect-of-the-three-frame-fallback-fallback--one_frame) |
| 2026-11-10 | **Keep a frame for slices lacking the start codon** (switch `frame_for_start_less_slice`, audit only, never released); add an annotated-start rule in PR 4c | [4b audit §5.4, §7](2026-11-10-phase-audit.md#54-slices-that-lack-the-start-codon) |
| 2026-11-10 | Alternative C-terminal isoforms go to a separate `psq_T1_alt_cterm.fasta` | [C-terminal audit §7](2026-11-10-alt-cterm-tiering-audit.md) |
| 2026-11-10 | PR 4c: annotated sites at TSL ≤ 2; one file per terminus (`psq_T1_alt_cterm`, `psq_T1_alt_nterm`, `psq_T1_alt_both`, each with `_orphan`); one fate `SUCCESS TIER 1 ALT TERMINUS`; frameshift test waived only for a slice ending at an internal annotated stop | [termini audit §7](2026-11-10-alt-termini-audit.md#7-decisions-needed) |
| 2026-11-10 | PR 4b2: fix the frameshift test (coding lengths) before PR 4c | [termini audit §6](2026-11-10-alt-termini-audit.md#6-regression-found-the-frameshift-test-ignores-per-slice-frames-proposed-pr-4b2) |
| 2026-11-11 | PR 3: uppercase genome sequence at retrieval; skip a slice with any base outside `ACGT` (as for N); keep `--mask` but deprecate it (no effect; remove in the next release); catch an exception per junction, record it and continue (in PR 7) | [masked-bases audit §7](2026-11-11-masked-bases-audit.md#7-decisions-requested) |
| 2026-11-11 | PR 5: join a peptide only at an end that is not already complete (stop codon, gene start codon); use an anchor only if it occurs once, start before end; anchor length 7; terminus comparison stays at 10 | [stitch audit §8](2026-11-11-stitch-anchor-audit.md#8-decisions-requested) |
| 2026-11-11 | PR 14: add `--uniprot FASTA` (canonical by sequence, then gene name; GTF canonical if neither) and a built-in Swiss-Prot merge of `psq_T1` + `psq_T1_alt_*` without `tr|` records; keep the API path without `--uniprot`; defer joining by coordinates | [UniProt evaluation §7](2026-11-11-uniprot-fasta-evaluation.md#7-decisions-2026-11-11) |
| 2026-11-12 | PR 9: per-gene GTF index and in-memory FASTA de-duplication (no output change). PR 15: parallel junctions with `-p/--processes` (default 1; forked workers return records, the main process writes them in junction order); more than one process only with `--uniprot`; serial on Windows | [parallel audit §7](2026-11-12-parallel-audit.md#7-decisions-2026-11-12) |
| 2026-11-12 | PR 6: keep a junction if its count ≥ the threshold; default `-r 1` (no filter); `-m` cut at the model's boundary and reported in reads; the count stays min(summed IJC, summed SJC), renamed `min_jc`. PR 7: create the output directory; refuse one holding JCAST output unless `--overwrite`; a failing junction writes nothing, gets `FAILURE (exception)`, and the run continues (exit status 0); caches in the user cache directory via `platformdirs` with `--cache_dir`, GTF cache keyed by path, size and mtime, written atomically, as a pickle. PR 8: `-c` writes the canonical on every path; a failed API lookup fails the junction (30 s timeout, 429 retried); `write_canonical`'s API fallback removed | [PR 6–8 audit §7](2026-11-12-thresholds-runs-canonical-audit.md#7-decisions-2026-11-12) |
| 2026-11-12 | Log: apply the read-count and FDR filters first (with PR 6); per-junction `junctions.tsv`, a small log with an end-of-run summary, sequences only with `--verbose`, no `~module~function` suffix (PR 16) | [log audit §6](2026-11-12-log-audit.md#6-decisions-2026-11-12) |
| 2026-11-12 | `--uniprot` mapping kept (sequence, then gene name, then GTF); gene names used only when unique in the GTF (PR 8); README documents a Swiss-Prot + TrEMBL FASTA (not expected for human); GTF canonical from `Ensembl_canonical` when no transcript passes the filter, as PR 17 after PR 16, audit first | [mapping audit §7](2026-11-12-uniprot-mapping-audit.md#7-decisions-2026-11-12) |
| 2026-11-12 | PR 17: when no transcript passes the GTF canonical rule, take the gene's `Ensembl_canonical` transcript if it is protein_coding, has start and stop codons and is not tagged `readthrough_transcript` (any TSL); the README says readthrough genes get no canonical, and why | [`Ensembl_canonical` audit §9](2026-11-12-ensembl-canonical-audit.md#9-decisions-2026-11-12) |
| 2026-11-13 | PRs 10–12: Python ≥ 3.11, tested on 3.11–3.15, with the dependency bounds checked; CI on Linux and macOS; version 0.4.0, to be released on GitHub after PRs 10–12; in the order 10, 11, 12; PR 12 also removes `--g_or_ln` and the unused `get_rest_nuc` | [packaging check §5](2026-11-13-packaging-check.md#5-decisions-2026-11-13) |
| 2026-11-13 | PR 12: only the rMATS files of the `-s` splice types must exist; with `-m`, the others present are read too, so the model fits on every junction in the folder (fitting on the selected types alone would move the chr15 `-s SE -m` cutoff from 26 to 6); the log lists every rMATS file read; remove the unused scikit-learn Gaussian mixture and scikit-learn | [command-line checks §3](2026-11-13-command-line-checks.md#3-decisions-2026-11-13) |
| 2026-11-13 | v0.5.0: `jcast` a group of commands (`run`, `gui`), `jcast RMATS GTF GENOME …` still running; `-s SE MXE` still accepted; `typer>=0.26,<1`; option spellings kept, no hyphenated aliases; the GUI runs `python -m jcast run …` in a child process; it shows the model before the run and the run summary after; junction browser in v0.5.x; no standalone bundles in v0.5.0 | [GUI and command-line scope §7](2026-11-13-gui-cli-scope.md#7-decisions-2026-11-13) |
| 2026-11-13 | PR 18: argparse's abbreviated long options (`--overw` for `--overwrite`) are not restored; long options are written in full | [GUI and command-line scope §9.3](2026-11-13-gui-cli-scope.md#93-abbreviated-long-options) |
| 2026-11-13 | PR 20: when the `-m` model has no fit, for any reason (also an exception, such as rMATS files without junctions), the run warns, falls back to `-r` (1 by default) and continues | [GUI and command-line scope §11.1](2026-11-13-gui-cli-scope.md#111-decision-2026-11-13) |
| 2026-11-14 | v0.5.0 widened before its release: the junction browser (was v0.5.x) becomes a view of each junction and its peptides, designed after a feasibility analysis; the `-m` read-count model is assessed (its fit to junction counts, a cutoff per splice type, the protein-coding isoforms a cutoff loses); each starts with an audit, and the release follows them | [REVIEW.md §9](../REVIEW.md#v050-the-gui-and-the-command-line), rows 22 and 23; HANDOFF.md |
| 2026-11-15 | PR 13: the GUI's choices that the scope did not settle are kept (`-r 1` until the model's cutoff is chosen; another folder resets the cutoff; the fit starts when the Model page opens; `-p` only with `--uniprot`; "required" shown once a run starts; Stop and quitting ask first; Delete removes only a run's files, if it wrote its log; settings not remembered) | [GUI and command-line scope §12.1](2026-11-13-gui-cli-scope.md#121-what-pr-13-changes) |
| 2026-11-15 | Row 23: testing the cutoff against mass spectrometry is deferred, as a future research item: no RNA-seq here has spectra from the same samples (not even ENCODE's), and the result would depend on the spectra's depth. Instead, measure at each cutoff the known protein-coding isoforms (Swiss-Prot, GENCODE) that a JCAST database reproduces. The whole genome's samples are iPSC-derived cardiomyocytes, DMSO against the NMD inhibitor NMDI14 | [REVIEW.md §5.6](../REVIEW.md#56--the-modeled-read-cutoff-is-a-weakly-identified-heuristic--new); HANDOFF.md |
| 2026-11-16 | Row 23: a cutoff for novel events, with one for annotated events, is favoured over a cutoff per splice type (not pursued); the GUI's Model page is to show what a cutoff keeps, with charts like the audit's, so users can choose; a third data set (the 2019 paper's deep iPSC differentiation series, rMATS 4.4.0) is to join the comparisons | [read-count model audit §§9, 11, 12](2026-11-15-read-count-model-audit.md#9-two-cutoffs-annotated-and-novel-events) |
| 2026-11-18 | Row 23: `-m` is kept. The pair of cutoffs (`-m` on novel events only) is to be decided after a quantms search of two deep human heart proteomes, with a database each way. A cutoff per splice type is dropped for now. The Model page is to show the database-size charts, the alternative region by read count, and the records by kind if quick. NMDI14 is inconclusive, and a science report on it is to come. The GUI remembers its settings in v0.5.0. | [read-count model audit §§13, 16](2026-11-15-read-count-model-audit.md#13-decisions-requested) |
| 2026-11-18 | Row 23: JCAST is to read rMATS's `--statoff` and one-group folders (an `NA` FDR passes the default `-q`; an empty count field counts 0 for its group; events with no inclusion or no skipping read are left out as the folder is read, so only the one-group events come back, and `-m` is fitted on what remains), and the documentation is to recommend `--statoff`. The maintainer will make the three folders again with it. The Model page shows only the charts that need no run for now; those from a run's records wait. | [read-count model audit §§15.7, 15.8, 16.2](2026-11-15-read-count-model-audit.md#157-decisions-requested) |

## Current tiers (after PR 4c)

Before tiering, each junction is trimmed to the gene's annotated coding
region, and exons entirely in the UTR are removed. Junctions whose event
lies outside that region are skipped. Each slice is then translated once,
in its own annotated frame; a slice with no annotated frame is not
translated.

| Tier | Meaning |
|---|---|
| **T1** | Both slices translate with no premature stop codon, and their coding lengths (from each slice's own frame) differ by a multiple of 3. |
| **T1 ALT TERMINUS** | At least one slice ends at an annotated alternative stop codon (C-anchored, joined at its N-terminus only) or starts at an annotated alternative start codon (N-anchored, joined at its C-terminus only). Files: `psq_T1_alt_cterm`, `psq_T1_alt_nterm`, `psq_T1_alt_both` (each with `_orphan`). |
| **T2** | As T1, but with a frameshift. |
| **T4** | One slice translates. The other is translated up to its first stop codon; that product is joined N-terminally if it is ≥ `ptc_threshold` × the translated slice's length, otherwise it is written as an orphan. |

T3 no longer exists. The full definitions are in
[REVIEW.md §3](../REVIEW.md#tier-classification-current-after-pr-4c).

## Golden-test history

`tests/test_pipeline.py::test_chr15_golden` pins every junction's fate and
every output file's contents. Its values at each stage:

| Stage | Commit | FAIL | T4 | Skip: low count | Skip: outside CDS | T1 | T2 | T3 | `psq_T1` records |
|---|---|---|---|---|---|---|---|---|---|
| v0.3.5a code (PR 1 golden) | `3ed17b0` | 360 | 236 | 297 | – | 185 | 25 | 2 | 152 |
| PR 4: stop-codon bounds | `da72a0a` | 330 | 261 | 297 | – | 189 | 26 | 2 | 164 |
| PR 4a: trim to the coding region | `89e1a82` | 240 | 259 | 297 | 90 | 190 | 27 | 2 | 169 |
| PR 4a2: RI slices | `bc7d894` | 238 | 258 | 297 | 92 | 189 | 28 | 3 | 169 |
| PR 4b candidate: per-slice frames, fallback kept | `9e27815` | 165 | 289 | 297 | 92 | 214 | 45 | 3 | 170 |
| PR 4b final: fallback and T3 removed | `2679edb` | 157 | 302 | 297 | 92 | 212 | 45 | – | 170 |
| PR 4b2: frameshift test from coding lengths | `a0f3338` | 157 | 302 | 297 | 92 | 227 | 30 | – | 171 |
| PR 4c: alternative termini (+ T1 ALT 25) | `e11cc86` | 153 | 289 | 297 | 92 | 219 | 30 | – | 173 |
| PR 5: joining to the canonical (+ T1 ALT 24) | `2142711` | 153 | 289 | 297 | 92 | 220 | 30 | – | 179 |
| PR 6: `-r 1` keeps every junction (+ T1 ALT 26) | `05a6931` | 207 | 403 | 0 | 119 | 301 | 49 | – | 244 |
| PR 8: `-c` writes the canonical on every path (fates and `psq_T*` unchanged) | `890d631` | 207 | 403 | 0 | 119 | 301 | 49 | – | 244 |

`psq_canonical.fasta` has 268 records at every stage up to PR 5, and 266
from PR 6: two genes whose junctions are now translated lost the canonical
their skipped junctions wrote (MED-14). From PR 8 it has 269, one for every
gene with a canonical. Up to PR 4c, no
`psq_T1.fasta` sequence matching an annotated protein is lost at any stage.
PR 5 removes one: CD276-201. It came only from a join at the wrong copy of
a repeated anchor
([stitch audit §9](2026-11-11-stitch-anchor-audit.md#9-implementation-and-verification-pr-5)).

PR 6 changes the test because its `read=1` now keeps every junction; with
`-r 2` the output is byte for byte that of PR 5
([PR 6–8 audit §9](2026-11-12-thresholds-runs-canonical-audit.md#9-implementation-and-verification-pr-6-05a6931)).

PR 14 leaves this test unchanged: it does not use `--uniprot`. The new
`tests/test_uniprot.py::test_chr15_uniprot` pins the same run with
`--uniprot` on a Swiss-Prot 2026_04 subset: `psq_T1` 186
records, `psq_canonical` 275, and `jcast_uniprot_merged.fasta` 837 (711
UniProt entries + 126 JCAST records)
([UniProt evaluation §8](2026-11-11-uniprot-fasta-evaluation.md#8-implementation-and-verification-pr-14)).
- **From PR 6 (`-r 1` keeps every junction):** `psq_T1` 253, `psq_canonical`
  272, and `jcast_uniprot_merged.fasta` 893.
- **From PR 8:** `psq_canonical` 276, one for every gene with a canonical.

PR 17 leaves both tests unchanged: the bundled Ensembl 89 GTF has no
`Ensembl_canonical` tag
([`Ensembl_canonical` audit §10](2026-11-12-ensembl-canonical-audit.md#10-implementation-and-verification-pr-17-8ddf8c3)).

## New defects found during the audits

- **HIGH-11**, found in the PR 4a audit: RI slices contained both flanking
  exons twice. Fixed in PR 4a2.
- **A5SS/A3SS on the `−` strand were trimmed at the wrong bound**, found in
  the PR 4a audit. Fixed in PR 4a.
- **The frameshift test ignored per-slice frames after PR 4b,** found in
  the PR 4c audit (§6). It misclassifies T1/T2 for 17 chr15 junctions.
  Fixed in PR 4b2.
- **`-c` does not write the canonical for genes whose junctions all
  succeed** (REVIEW.md MED-14; 2 chr15 genes), found while verifying
  PR 4c. The PR 6–8 audit found that masked junctions never wrote it
  either. Fixed in PR 8.
- **Hard-masked genomes silently drop junctions**, found in the PR 3 audit
  (§4.2). The README recommended them as a workaround for soft-masking;
  on chr15, 206 of 1,105 junctions (20 successes) are skipped as masked.
- **One exception aborts the whole run**, found in the PR 3 audit: the
  per-junction loop catches nothing. Fixed in PR 7.
- **Joins at an already-complete end returned the canonical**, found in
  the PR 5 audit. A T4 peptide ending at its stop codon was joined
  C-terminally; a peptide starting at the gene's start codon was joined
  N-terminally. 10 of 12 wrong joins on chr15 gave exactly the canonical.
  Fixed in PR 5.
- **A half-written GTF cache was read without error**, found in the
  [PR 6–8 audit](2026-11-12-thresholds-runs-canonical-audit.md) (§3.5). A run started while another wrote the cache
  read half the genes. Fixed in PR 7.
- **A refused UniProt API request silently orphaned every record**, found
  in the [PR 6–8 audit](2026-11-12-thresholds-runs-canonical-audit.md) (§4.3). HTTP 429 and other refusals were
  taken as "no canonical". Fixed in PR 8.
- **`-m` cut two reads above the model's boundary**, found in the
  [PR 6–8 audit](2026-11-12-thresholds-runs-canonical-audit.md) (§2.3): 3,410 of 131,233 whole-genome junctions
  dropped. Fixed in PR 6.
- **The T4 length threshold was mis-documented.** It decides whether the
  truncated product is joined to the canonical, not whether it is written.
  Corrected in the README and REVIEW.md with PR 4b.
