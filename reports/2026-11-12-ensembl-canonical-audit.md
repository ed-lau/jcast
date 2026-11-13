# Audit: the GTF canonical from GENCODE's `Ensembl_canonical` tag (PR 17)

| | |
|---|---|
| **Date** | 2026-11-12 |
| **Status** | Decisions approved 2026-11-12 (§9): option 1, the fallback without readthrough transcripts. **Implemented in PR 17** (`8ddf8c3`, §10). Audit measured at `ae3d937`, with a prototype (two audit-only variants) run on the whole genome (§6). |
| **Question** | When no transcript passes the GTF canonical rule (protein_coding, TSL ≤ 2, annotated start and stop codons), should JCAST take the transcript GENCODE tags `Ensembl_canonical`? Scheduled as PR 17 by the [mapping audit](2026-11-12-uniprot-mapping-audit.md#7-decisions-2026-11-12) decision 3, which counted 173 of 384 coding genes that it would cover. |
| **Reproduce** | `python reports/canonical_audit/audit_canonical.py genes --gtf <gtf> --run <output dir> --rmats <dir>` (≈ 1 min with the GTF cache) · `… compare --before <dir> --after <dir> --label <name>` · PR 17: `python reports/canonical_audit/verify_17.py --before <checkout> --after <checkout> [--genome_before <dir> --genome_prototype <dir> --genome_after <dir>]` (≈ 4 min) |
| **Inputs** | The user's whole-genome data (599,644 junctions, GENCODE v50, GRCh38, Swiss-Prot human with isoforms 2026_04), run with `-m -c --uniprot -p 16` at PR 16 (`a85f8fc`; output identical to `f181966`, [log audit §7](2026-11-12-log-audit.md#7-implementation-and-verification-pr-16-a85f8fc)). Not bundled. |
| **Evidence** | [canonical_audit/genome/summary.txt](canonical_audit/genome/summary.txt) · [genes.tsv](canonical_audit/genome/genes.tsv) (one row per gene of the run without a canonical) · [compare_all.txt](canonical_audit/genome/compare_all.txt) · [compare_no_readthrough.txt](canonical_audit/genome/compare_no_readthrough.txt) · [prototype.diff](canonical_audit/prototype.diff) · PR 17: [canonical_audit/verify_17/verify_17.txt](canonical_audit/verify_17/verify_17.txt) · [compare_pr17.txt](canonical_audit/verify_17/compare_pr17.txt) |

## 1. Summary

- **Where the GTF canonical matters.** With `--uniprot`, it is used only
  for a gene that matches no Swiss-Prot entry by sequence or by name. The
  API path never uses it, and GTF-only runs (`params.use_gtf_only`) are
  not reachable from the command line. A record joined to a GTF canonical
  is never added to `jcast_uniprot_merged.fasta`, which takes only records
  joined to a Swiss-Prot canonical.
- **388 protein-coding genes of the run have no canonical** (384 in the
  mapping audit, before PR 8's unique-name rule; the 4 more are its
  same-name coding loci). They hold 384 orphan records.
- **370 of the 388 are readthrough genes.** Their `Ensembl_canonical`
  transcript is a readthrough transcript, spanning two or more adjacent
  genes (GENCODE tag `readthrough_transcript`):
  - 211 have a CDS only in nonsense-mediated-decay transcripts, and their
    `Ensembl_canonical` is one of them;
  - 159 have protein_coding transcripts, of TSL 3–5 or NA, or without both
    codons.

  They hold 376 of the 384 orphan records.
- **The other 18 are ordinary genes,** mostly recent loci without a gene
  name or a TSL:
  - 17 have a complete protein_coding `Ensembl_canonical` transcript;
    8 of these are MANE Select (e.g. HMSD, TSL 3);
  - they hold 8 orphan records, in 4 genes.
- **A readthrough canonical would make fusion records.**
  - rMATS reports each event under one gene only (all 260,199 SE events
    are distinct).
  - 45% of the readthrough genes' SE events skip an exon annotated in
    another protein-coding gene, typically one of the genes the readthrough
    spans.
  - Their slices, joined to a readthrough canonical, would read as fusion
    proteins, not as isoforms of the gene the event belongs to.
- **The mapping audit's 173** (177 now, with PR 8's 4) are the genes
  whose `Ensembl_canonical` transcript is protein_coding, complete or not:
  159 readthrough and 18 ordinary. Its CDS is complete (start and stop
  codons annotated) for 99 of them: 82 readthrough and 17 ordinary.
- **Recommended:** the fallback for complete protein_coding
  `Ensembl_canonical` transcripts that are **not** readthrough. It gives
  17 genes a canonical and changes almost nothing else (§6).

## 2. How the GTF canonical is chosen, and where it is used

`Sequence.get_annotated_transcripts_gtf`
([sequences.py](../jcast/sequences.py)) takes, among the gene's
transcripts:
- of type `protein_coding`, with CDS rows;
- with transcript support level ≤ `params.tsl_threshold` (2);
- with an annotated start codon and stop codon.

The longest coding sequence is the canonical. If no transcript qualifies,
the gene has no GTF canonical.

Where it is used:
- **`--uniprot`:** the gene's canonical is a Swiss-Prot entry matched by
  sequence (against every CDS-bearing transcript of the gene, whatever its
  type or TSL), else by gene name if the name is unique in the GTF, else
  the GTF canonical, else none. In the whole-genome run: 13,072 by
  sequence, 372 by name, **102 from the GTF**, 3,522 none.
- **Without `--uniprot`:** the UniProt API; the GTF canonical is not used.
- **`params.use_gtf_only`:** GTF only. Used by the chr15 golden test;
  there is no command-line option for it.
- **Records:** a slice joined to a GTF canonical is named
  `<gene id>|<symbol>|…` and goes to its tier's file. It is not added to
  `jcast_uniprot_merged.fasta`. With `-c`, the GTF canonical itself goes
  to `psq_canonical.fasta`.

Nothing else uses the chosen transcript. The coding region, frames and
annotated start and stop sites have their own TSL ≤ 2 lookups in
`junctions.py`.

## 3. The genes without a canonical

[summary.txt](canonical_audit/genome/summary.txt). Genes of the run
(17,068) whose canonical source is `none` in `uniprot_canonicals.tsv`:

| Genes | Why the GTF rule gives none | `Ensembl_canonical` transcript | Readthrough | Orphan records |
|---|---|---|---|---|
| 211 | CDS only in transcripts not of type protein_coding | nonsense_mediated_decay | 211 | 0 |
| 153 | protein_coding transcripts of TSL 3–5 or NA only | protein_coding | 135 | 320 (312 readthrough) |
| 24 | TSL ≤ 2 protein_coding transcripts lack a start or stop codon | protein_coding | 24 | 64 |
| **388** | protein-coding genes | | **370** | **384** (376 readthrough) |
| 3,134 | other genes: no transcript with a CDS (3,133), or a CDS only in a non-coding transcript (1) | lncRNA, pseudogenes, … | – | 0 |

The `Ensembl_canonical` transcripts of the 388:
- **type:** protein_coding 177, nonsense_mediated_decay 211;
- **TSL:** 1: 7; 2: 92; 3: 51; 4: 29; 5: 103; NA: 106;
- **codons:** start and stop 242; stop only 94; start only 34; neither 18;
- **MANE Select:** 8, all ordinary genes.

Every one of the 388 has an `Ensembl_canonical` transcript. Ensembl tags
one per gene, coding or not.

**The 18 ordinary genes**
([genes.tsv](canonical_audit/genome/genes.tsv), `readthrough` = no):
- **Names:** HMSD, ZKSCAN8P1, RHOXF1P3, CNTNAP3C, and 14 loci without a
  gene name (their symbol is the gene id).
- **Why none:** all 18 have only protein_coding transcripts of TSL 3–5 or
  NA. TSL NA means the transcript was not analysed for support, as for
  most recent models; it does not mean unsupported.
- **Complete CDS:** 17. The exception is CNTNAP3C (`cds_end_NF`).
- **Size:** 3 of them encode 12–25 residues.
- **Orphan records:** 8, in ENSG00000288678, ENSG00000288715,
  ENSG00000288841 and ENSG00000289766 (2 each).

## 4. Readthrough genes and rMATS events

**rMATS reports an event once.** In the run's `SE.MATS.JC.txt`, the
260,199 events are all distinct by their eight coordinates. An event in
a region shared by a readthrough gene and the genes it spans is reported
under one of them only.

**The readthrough genes' events belong to other genes as much as to the
readthrough.**
- The 370 readthrough genes without a canonical have 2,793 SE events.
- For 1,247 (45%), the skipped exon is an annotated exon of another
  protein-coding gene.
- Their translated slices are orphans today: 362 of their 376 orphan
  records have no record of the same event under another gene.

**What a readthrough canonical would do.** Joined to it, such a slice
becomes the readthrough protein with the alternative exon: a fusion of
two genes' proteins, around an event that may belong to only one of
them. Readthrough transcripts are mostly NMD targets or low-support
models (§3), and no Swiss-Prot entry matched these genes, by sequence
or by name.

The better join for these events would be the canonical of the gene the
event lies in, i.e. joining by genomic coordinates. That was evaluated
and deferred in the
[UniProt evaluation §4](2026-11-11-uniprot-fasta-evaluation.md#4-joining-by-coordinates-question-1).

## 5. `Ensembl_canonical` against the current rule, where the rule gives a canonical

For context only; no option below changes these genes.

- **GENCODE v50 overall:** 20,131 protein-coding genes; the rule takes a
  transcript for 17,162. Its choice has the same coding sequence as the
  `Ensembl_canonical` transcript for 14,002. For 1,937, the
  `Ensembl_canonical` transcript passes the rule but is shorter. For
  1,223, it does not pass the rule.
- **In the run:** of the 102 genes with a GTF canonical, 86 have the same
  coding sequence as their `Ensembl_canonical` transcript; 16 differ.

Making `Ensembl_canonical` the first choice would change those 16 genes
under `--uniprot`, and 3,160 genes in GTF-only runs. It was not asked
for, and is not proposed.

## 6. Effect on the output (prototype)

A prototype ([prototype.diff](canonical_audit/prototype.diff), on
`ae3d937`) adds the fallback to `get_annotated_transcripts_gtf`, behind
an audit-only switch:
- **`all`:** if no transcript passes the rule, take the transcript tagged
  `Ensembl_canonical`, of type protein_coding, with start and stop codons,
  at any TSL;
- **`no_readthrough`:** the same, unless it is tagged
  `readthrough_transcript`.

Each variant ran on the whole genome with `-m -c --uniprot -p 16`, and
was compared with the PR 16 run
([compare_all.txt](canonical_audit/genome/compare_all.txt),
[compare_no_readthrough.txt](canonical_audit/genome/compare_no_readthrough.txt)).
Under the `all` switch, both golden tests pass unchanged: the Ensembl 89
chr15 GTF has no such tag.

| | PR 16 | `no_readthrough` (option 1) | `all` (option 2) |
|---|---|---|---|
| Genes whose canonical source changes (`none` → `gtf`) | – | **17** | 99 (82 readthrough) |
| Their records, joined / orphan: before | | 0 / 8 | 0 / 200 |
| … after | | 2 / 3 | 107 / 34 |
| `psq_canonical.fasta` | 13,474 | 13,491 | 13,573 |
| `psq_T1.fasta` | 22,521 | 22,522 | 22,542 |
| `psq_T2.fasta` | 6,458 | 6,458 | 6,464 |
| `psq_T4.fasta` | 42,929 | 42,930 | 43,009 |
| Orphan files (T1, T2, T4) | 3,675 / 399 / 9,627 | 3,672 / 399 / 9,625 | 3,639 / 393 / 9,503 |
| `jcast_uniprot_merged.fasta`, `psq_isoforms.tsv` | | unchanged | unchanged |
| Fates in `junctions.tsv` | | unchanged (202 rows differ, in canonical and files only) | unchanged (2,097 rows) |

**Notes:**
- **Option 1, the 17 genes:**
  - Each gets its canonical in `psq_canonical.fasta` (`-c`).
  - Of their 8 orphan records, 2 now join (ENSG00000288841,
    ENSG00000289766). Three are identical to the new canonical once
    joined, so they go to `psq_canonical.fasta`, as for any gene. Three
    stay orphans.
- **Option 2, the 82 readthrough genes:** about 105 records now join to
  fusion canonicals (§4), most of them in T4.
- **Neither option changes** the merged database or any fate. They change
  only which canonical a slice can join to, and whether it is written.
- **The run time** is unchanged (3 min 03 s translating, against 3 min 02 s).

## 7. Options

1. **The fallback, without readthrough transcripts** (recommended). If no
   transcript passes the rule, take the gene's `Ensembl_canonical`
   transcript when it:
   - is of type protein_coding;
   - has annotated start and stop codons;
   - is not tagged `readthrough_transcript`.

   TSL is not considered. 17 genes get a canonical, 8 of them MANE Select.
   Readthrough genes keep no canonical, by design, and the README says so.
2. **The fallback for every complete protein_coding `Ensembl_canonical`
   transcript,** readthrough included, as scheduled. 99 genes, 82 of them
   readthrough. Their records would join to fusion proteins (§4).
3. **Accept incomplete or NMD `Ensembl_canonical` transcripts** (variant
   "any CDS", 388 genes). Not recommended. Their translations would be
   fragments, or the products of NMD transcripts. They would need a frame
   taken from the first CDS row instead of the start codon.
4. **No change.** Document that genes whose transcripts are all TSL 3–5 or
   NA have no GTF canonical. The 17 genes stay without one.

Whichever option is taken:
- **Tags are matched whole,** after splitting on commas: GENCODE also has
  `Ensembl_canonical_extended`, and gtfparse joins repeated tags with
  commas.
- **Other GTFs:** recent Ensembl GTFs carry the same tags. A GTF without
  them, such as the bundled Ensembl 89 chr15, is unaffected, so neither
  golden test changes.
- **Tests:** a small synthetic GTF, in the style of `tests/test_runs.py`.

## 8. Decisions requested

1. **Scope:** option 1, the fallback without readthrough transcripts
   (recommended), option 2 as scheduled, or option 4, no change?
2. **The README** to say that readthrough genes get no canonical, and why
   (recommended with option 1)?

## 9. Decisions (2026-11-12)

1. **Option 1, the fallback without readthrough transcripts.** Approved,
   as PR 17. Readthrough genes are not JCAST's main purpose; their events
   are better left as orphans than joined to fusion proteins.
2. **The README says that readthrough genes get no canonical, and why.**
   Approved, with PR 17.

## 10. Implementation and verification: PR 17 (`8ddf8c3`)

**What changed** ([sequences.py](../jcast/sequences.py)):
- `Sequence.get_annotated_transcripts_gtf` applies the rule as before.
  If no transcript passes, it takes the CDS rows of the gene's transcript
  tagged `Ensembl_canonical` (`_ensembl_canonical_cds`):
  - of type protein_coding;
  - not tagged `readthrough_transcript`;
  - at any TSL.

  The same start- and stop-codon test then applies
  (`_annotated_transcripts`, the rule's loop, moved out unchanged).
- Tags are split on commas and matched whole. A GTF without a `tag`
  column gets no fallback.
- This is the prototype's `no_readthrough` variant without its switch.
  The fallback rows are computed only for a gene the rule gives nothing;
  the prototype computed them for every gene.

**Tests** ([tests/test_gtf_canonical.py](../tests/test_gtf_canonical.py),
10 tests, on a synthetic GENCODE-style GTF):
- a gene with TSL NA transcripts only gets the tagged transcript, not a
  longer untagged one, translated from its start codon;
- no canonical when the tagged transcript is a readthrough or an NMD
  transcript, lacks a start or a stop codon, or carries only
  `Ensembl_canonical_extended`;
- a gene with a transcript that passes the rule keeps the rule's choice;
- no fallback for a GTF without tags. The bundled human and rat GTFs
  (Ensembl 89 and 100) carry no `Ensembl_canonical` tag.

With the fallback switched off, the 3 tests of the tagged gene fail. All
148 tests pass. Both golden tests are unchanged.

**Verification**
([verify_17.py](canonical_audit/verify_17.py),
[verify_17.txt](canonical_audit/verify_17/verify_17.txt)). The code
before is `21e178f` and the code after is its working tree, each frozen
in a git worktree.
- **chr15** (`-c --uniprot`, with `-r 1` and with `-m`): every output
  file is identical, and so is the log's summary.
- **Whole genome** (`-m -c --uniprot -p 16`):
  - The run after is byte-identical to the prototype's `no_readthrough`
    run (§6): all 16 output files but the log and `model.pdf`, and the
    log's summary.
  - Against the run before
    ([compare_pr17.txt](canonical_audit/verify_17/compare_pr17.txt)), as
    in §6:
    - 17 genes go from no canonical to a GTF canonical;
    - `psq_canonical` 13,474 → 13,491; `psq_T1` +1 and `psq_T4` +1;
      their orphan files −3 and −2;
    - 202 rows of `junctions.tsv` differ, in the canonical, the files
      and the records written only. Fates, frames and coding regions are
      unchanged;
    - `jcast_uniprot_merged.fasta` and `psq_isoforms.tsv` are identical.
  - **Still without a canonical:** 371 protein-coding genes, holding
    376 orphan records. 370 are readthrough genes, and hold all 376. The
    other is CNTNAP3C, whose `Ensembl_canonical` CDS is incomplete
    (`cds_end_NF`).
  - **Time:** no measurable change. Run back to back on a busy machine
    (load average above 30), the junctions took 4 min 12 s to translate
    before and 4 min 20 s after. A first run after took 3 min 47 s, also
    under load; the prototype took 3 min 03 s on a quieter machine. Both
    the repeated runs reproduce the first runs' output exactly.
