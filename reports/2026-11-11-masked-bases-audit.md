# Audit: soft-masked and ambiguous genome bases (PR 3; REVIEW.md HIGH-2/HIGH-3)

| | |
|---|---|
| **Date** | 2026-11-11 |
| **Status** | **Implemented in PR 3** (§8). Decisions in §7 approved 2026-11-11 (option A; deprecate `--mask`; per-junction exceptions in PR 7). Baseline is commit `a784030`. |
| **Code** | `helpers.get_local_nuc`, `helpers.get_complementary`, `helpers.make_pep`, the `--mask` gate in `main._translate_one` |
| **Environment** | Python 3.11.13, pandas 2.3.1, numpy 1.26.4, biopython 1.78, gtfparse 1.3.0 |
| **Reproduce** | chr15: `python reports/mask_audit/audit_mask.py --softmasked Homo_sapiens.GRCh38.dna_sm.chromosome.15.fa.gz` (add `--verify` at the PR 3 commit for §8) (≈ 10 min, no network once the file is downloaded). Genome survey: `python reports/mask_audit/survey_genome.py --genome GRCh38.primary_assembly.genome.fa --gtf gencode.v50.primary_assembly.annotation.gtf --out reports/mask_audit/gencode_v50` (≈ 5 min) |
| **Evidence** | [mask_audit/baseline/summary.txt](mask_audit/baseline/summary.txt) · [mask_audit/baseline/junctions.tsv](mask_audit/baseline/junctions.tsv) · [mask_audit/gencode_v50/](mask_audit/gencode_v50/) |

## 1. Summary

- **A soft-masked genome aborts a run.**
  - On chr15, 147 of 1,105 junctions (82 genes) raise `KeyError` in a
    default run. With `use_gtf_only` the count rises to 393 (130 genes).
  - The per-junction loop catches nothing, so the run stops at the first
    one: MXE 430, the third junction processed.
  - Soft-masked genomes are common: Ensembl `dna_sm` and UCSC `hg38.fa`
    are both soft-masked. On chr15, 712 of 2,678 protein-coding transcripts
    (290 genes) have at least one soft-masked base in their CDS.
- **Soft-masking carries no sequence information that JCAST needs.** Ensembl's
  `dna_sm` chr15 equals its unmasked `dna` chr15 once uppercased. The
  junctions that do not raise get exactly their unmasked fates.
- **Hard-masked genomes are not a safe workaround.** The README currently
  recommends them. In a default run they don't crash, but 206 of 1,105
  chr15 junctions (19%) are silently skipped as "masked". Those include
  20 that were successes (13 T1, 6 T1 ALT TERMINUS, 1 T2). With
  `use_gtf_only`, the same junctions are skipped and a further 239 raise.
- **Ambiguity codes don't occur in the human reference.**
  - GENCODE's GRCh38 primary assembly is unmasked: no lowercase bases and
    no IUPAC codes other than N.
  - None of its 151 Mb of N lies in any of the 370,915 annotated CDS
    (GENCODE v50).
  - For human, the IUPAC policy therefore matters only for robustness:
    other organisms, custom assemblies and personal genomes.

## 2. The defects

`make_pep` looks codons up in a dictionary of the 64 uppercase codons, and
`get_complementary` maps only `A/C/G/T`. Any other character raises
`KeyError`. The fallback in `get_complementary` (`comp += 'X'`) follows a
`raise` and is unreachable (HIGH-3). Reproduced at `a784030`:

| Input | Strand | Result |
|---|---|---|
| `ATGAAATAG` | + | `'MK'` |
| `atgaaatag` | + | `KeyError: Seq('atg')` |
| `ATGNNNTAG` | + | `KeyError: Seq('NNN')` |
| `ATGRAATAG` | + | `KeyError: Seq('RAA')` |
| `ctatttcat` | − | `KeyError: 'Invalid nucleotide: t'` |
| `CTARTTCAT` | − | `KeyError: 'Invalid nucleotide: R'` |

Translation goes through `make_pep` in three places:

1. The two slices (`Sequence.translate`).
2. The alternative-termini path (`route_alternative_termini`), which
   uppercases for its own codon checks but passes the original sequence to
   `make_pep`.
3. The canonical from the GTF (`translate_annotated_transcriptfrom_gtf`,
   only with `use_gtf_only`).

The `--mask` gate (default `NnXx`) protects only the slices, only against
the listed letters, and only by skipping the junction. It does not protect
the GTF canonical, which is why a hard-masked genome still crashes with
`use_gtf_only` (§4.2).

## 3. Survey: GENCODE GRCh38 primary assembly

Genome `GRCh38.primary_assembly.genome.fa` (GENCODE; md5
`49bdb80d21a64dcb16acfc941843356e`). Annotation
`gencode.v50.primary_assembly.annotation.gtf` (Ensembl 116; md5
`d4d28a69afdfc3c13cc08aba4217eab5`). Every base was counted, and every CDS
row intersected with every non-ACGT and lowercase base.
[survey_genome.py](mask_audit/survey_genome.py) was checked on a synthetic
genome first.

| | |
|---|---|
| Sequences | 194 (25 chromosomes, 169 scaffolds); 3,099,750,718 nt |
| Uppercase A/C/G/T | 2,948,627,755 |
| N | 151,122,963 (4.88%), in 45 sequences |
| Other IUPAC codes (R, Y, K, M, S, W, B, D, H, V) | 0 |
| Lowercase (soft-masked) | 0 |
| CDS transcripts with any N, IUPAC or lowercase base | 0 of 370,915 (all CDS); 0 of 278,484 protein-coding; 0 of 19,229 MANE Select |

The genome survey files (`composition.tsv`, `iupac_sites.tsv`,
`cds_affected.tsv`) are in [mask_audit/gencode_v50/](mask_audit/gencode_v50/).
The last two are empty apart from their header. Only this distribution was
surveyed. Whether UCSC or NCBI builds of GRCh38 contain IUPAC codes was not
measured.

## 4. chr15: what a masked genome does to a run

[audit_mask.py](mask_audit/audit_mask.py) runs every chr15 junction
through `main._translate_one` with three genomes:

- the unmasked test genome (Ensembl 89 `dna`);
- the soft-masked genome of the same release (Ensembl 89 `dna_sm`; md5
  `0c7f6f1ce8ce5138a8cbe999f8a83276`), which is identical to the unmasked
  one apart from case, with 43.8% of chr15 lowercase;
- a hard-masked genome derived from it (every soft-masked base → N).

Each genome is run twice: with the canonical translated from the GTF
(`use_gtf_only`), and as a default run, where the canonical comes from
UniProt. UniProt is stubbed out in the audit, so in the second mode only
the slices are translated. A real run stops at the first exception; the
audit catches each one to count them all.

### 4.1 Soft-masked

| Mode | Junctions raising `KeyError` | Genes | First (run order) |
|---|---|---|---|
| Default (slices only) | 147 | 82 | MXE 430 |
| `use_gtf_only` (slices and canonical) | 393 | 130 | MXE 430 |

- **Where the 393 raise:** 73 in a slice's translation and 320 in the GTF
  canonical's. 216 of them fail in `get_complementary` (`−` strand) and 177
  directly in `make_pep`.
- **What those 393 junctions are when unmasked:** 74 T1, 10 T1 ALT
  TERMINUS, 7 T2, 119 T4, 29 FAILURE, 120 low coverage, 34 outside the
  coding region. Translation runs before the coverage filter, so
  low-coverage junctions crash too.
- **The other 712 junctions:** each gets exactly its unmasked fate.

Soft-masked bases in chr15 coding sequences (Ensembl 89 GTF, CDS plus stop
codon):

| Transcripts | Total | With ≥ 1 soft-masked CDS base | Soft-masked CDS bases |
|---|---|---|---|
| protein-coding | 2,678 (583 genes) | 712 (290 genes) | 76,172 of 3,487,082 (2.18%) |
| protein-coding, TSL 1–2 | 1,389 (520 genes) | 420 (214 genes) | 32,687 of 2,272,536 (1.44%) |

### 4.2 Hard-masked

| Mode | Raising | Skipped as masked | Of those skipped, successes when unmasked |
|---|---|---|---|
| Default | 0 | 206 (19%) | 20 (13 T1, 6 T1 ALT TERMINUS, 1 T2) |
| `use_gtf_only` | 239 (76 genes), first MXE 431 | 206 | 20 (14 T1, 5 T1 ALT TERMINUS, 1 T2) |

- **Default mode, the 206 skipped junctions** were, unmasked: 83 T4,
  43 FAILURE and 60 low coverage, plus the 20 successes.
- **`use_gtf_only`:** the same 206 junctions are skipped, because the gate
  runs first. Of the rest, 239 raise: the GTF canonical is translated
  without passing the gate, so it fails on `N` (e.g. `KeyError: Seq('NNN')`).
- **This loss is inherent.** Sequence under N is unknown, so no fix can
  recover these junctions. The README's advice to use a hard-masked genome
  should be withdrawn (§6.4).

## 5. Policy options for non-ACGT bases

Lowercase needs no policy. It is annotation, and uppercasing is lossless
(§4.1). Ambiguous bases (N, and IUPAC codes where present) do need one:

| Option | Slice with an ambiguous base | Pros | Cons |
|---|---|---|---|
| **A. Skip (recommended)** | Junction skipped (`SKIPPED_MASKED`), as N is today | Same as the current N behaviour; no invented residues; output unchanged on every genome surveyed | Loses the junction (on GRCh38, no CDS is affected) |
| B. Translate with `X` | Ambiguous codon → `X`; unambiguous ones resolved, as Biopython does (`CTN` → L, `RAA` → X) | Keeps the junction | `X` already means "stop" inside `make_pep`; `X` residues in a search database behave differently between search engines; the `--mask` gate already skips N, so A is needed anyway for consistency |
| C. Raise a clear error | Run stops with the base and position | Explicit | Same as today in effect: one base aborts a whole run |

With any option, `make_pep` and `get_complementary` should never raise.
`get_complementary` would use Biopython's IUPAC-aware
`reverse_complement`. A codon that `make_pep` cannot translate would make
the translation untranslatable (`''`), the same outcome as a premature
stop. This only arises outside the gate, i.e. for the GTF canonical. An
empty canonical is handled exactly like a missing UniProt entry, so
`use_gtf_only` can no longer crash on hard-masked genomes.

## 6. Proposed change (PR 3)

1. **Uppercase at retrieval.** `get_local_nuc` is the only function that
   reads the genome. It returns the slice uppercased (`SeqRecord.upper()`).
   This covers the slices, the GTF canonical and the alternative-termini
   path at once. Uppercasing the whole genome at load would also work, but
   would briefly double the memory used by a 3 GB genome.
2. **Option A for ambiguous bases.** The gate skips a slice containing any
   base outside `ACGT`, which subsumes the default `NnXx`.
3. **`make_pep` and `get_complementary` never raise** (§5).
   `get_complementary` becomes `Bio.Seq.reverse_complement` (fixes HIGH-3).
   An untranslatable codon returns `''`.
4. **README:** soft-masked genomes are supported and give exactly the
   unmasked output. Remove the known issue and the "hard-masked" advice.
   Recommend unmasked or soft-masked genomes, and warn that hard-masked
   genomes skip every junction touching a repeat (19% on chr15).
5. **Tests:**
   - unit tests for the six inputs in §2 on both strands;
   - an integration test translating chr15 junctions from an in-memory
     lowercased copy of the test genome (MXE 430, the first to crash, plus
     one junction on each strand), asserting the unmasked fates and
     sequences;
   - a gate test with an IUPAC base.
6. **Verification:** `audit_mask.py --verify`. For both the unmasked and
   the soft-masked genome, the working tree must reproduce the baseline
   unmasked run exactly, in every fate and every FASTA record. The golden
   test must be unchanged.

**Expected output change:** none on unmasked genomes. Soft-masked genomes
change from aborting to producing output identical to the unmasked genome.

## 7. Decisions requested

**Approved 2026-11-11:** 1 → option A; 2 → (a), deprecate and remove in the next
release; 3 → yes, in PR 7.

1. **Ambiguous-base policy:** option A (skip, as for N), B or C?
2. **`--mask` option:** with every non-ACGT base skipped, the option has no
   remaining effect.
   - (a) **Recommended:** keep accepting it, warn that it is deprecated,
     and remove it in a later release.
   - (b) Keep it as a way to *add* letters, which has no use once
     everything outside ACGT is skipped.
   - (c) Remove it now. That breaks scripts that pass it.
3. **Separate from this PR:** should an exception in one junction stop the
   whole run? Any future unexpected error would also abort everything,
   with no partial output. Proposal: catch it per junction, log it, record
   a new `FAILURE (exception)` fate, and continue. That is MED-level
   robustness, better placed in PR 7 than mixed into this fix.

## 8. Implementation and verification (PR 3)

Implemented as proposed in §6, with the decisions in §7:

- **Retrieval:** `helpers.get_local_nuc` returns the genome slice uppercased.
- **The gate** in `main._translate_one` skips a junction if either slice
  has a base outside `ACGT`. The fate now reads `SKIPPED. Sequence
  discarded due to masked or ambiguous bases.`
- **`helpers.make_pep`** uppercases its input and returns `''` for a codon
  it cannot translate.
- **`helpers.get_complementary`** is `Bio.Seq.reverse_complement`.
- **`--mask`** defaults to `None`. If given, it is ignored and a
  deprecation warning is logged.

**Verification** (`audit_mask.py --verify`; evidence in
[mask_audit/final/](mask_audit/final/)):

| Check | Result |
|---|---|
| Unmasked genome, working tree vs baseline: fates and FASTA files | identical (1,105 fates; all 10 FASTA files) |
| Soft-masked genome, working tree vs baseline **unmasked**: fates and FASTA files | identical |
| Exceptions with the working tree (unmasked, soft-masked, hard-masked; with and without the GTF canonical) | 0 |
| Hard-masked, default run: junctions skipped as masked | the same 206 as before; only the fate wording changed |
| Hard-masked, `use_gtf_only`: the 239 junctions that raised | now get ordinary fates (59 T1, 6 T1 ALT TERMINUS, 6 T2, 53 T4, 8 FAILURE, 73 low coverage, 34 outside the coding region) |
| Golden end-to-end test | unchanged |

**Regression tests** ([tests/test_masking.py](../tests/test_masking.py)).
14 of the 15 fail at `819cc53`; the plain uppercase `make_pep` case passes
on both.

| Test | Asserts |
|---|---|
| `test_make_pep_never_raises` | the six inputs of §2: lowercase translates as uppercase, and ambiguous codons give `''` |
| `test_get_complementary_iupac` | IUPAC bases are complemented and case is kept |
| `test_softmasked_genome_gives_unmasked_output` | MXE 430, SE 5739, SE 28030, SE 6592 on an all-lowercase chr15: same fate and FASTA records as unmasked |
| `test_slice_with_ambiguous_base_is_skipped` | `N`, `n`, `R` or `y` inside a slice: `SKIPPED_MASKED` |
