# Evaluation: a user-supplied UniProt FASTA, and joining by coordinates

| | |
|---|---|
| **Date** | 2026-11-11 |
| **Status** | Evaluation. Recommendations in §5; decisions requested in §6 and approved 2026-11-11 (§7), to be implemented as PR 14. No code changed. Measured at `8b0c024`. |
| **Question** | (1) Should slices be joined to the canonical by genomic coordinates instead of anchor sequences? (2) Should the canonical come from a UniProt FASTA supplied by the user instead of per-gene API calls, with the post-run renaming against Swiss-Prot isoforms (currently an R script) built in? |
| **Reproduce** | `python reports/uniprot_eval/evaluate_uniprot_fasta.py --uniprot <Swiss-Prot human FASTA with isoforms>` (≈ 2 min; needs the local `cache/uniprot-cache.db`) |
| **Inputs** | UniProt release 2026_04 (4 Nov 2026), Swiss-Prot human with isoforms, downloaded from `https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=(proteome:UP000005640)+AND+(reviewed:true)`: 42,545 entries, 20,414 canonical. Ensembl 89 chr15 test set. API answers from the local cache, which an earlier run filled with 6,519 genes, 184 of them chr15 test genes. |
| **Evidence** | [uniprot_eval/summary.txt](uniprot_eval/summary.txt) · [uniprot_eval/genes.tsv](uniprot_eval/genes.tsv) · [uniprot_eval/records.tsv](uniprot_eval/records.tsv) |

## 1. Summary

- **Joining by coordinates is feasible but not worth it now.**
  - Coordinates exist for every GTF canonical, and for 179 of 184 UniProt
    canonicals (97%), whose sequence equals a GTF-translated transcript of
    the gene.
  - After PR 5, the sequence join makes no wrong joins on chr15 and misses
    3 of 994. Joining by coordinates would recover at most those 3, at the
    cost of a second join path.
  - It becomes cheap if decision 2 is adopted, because mapping a FASTA
    entry by sequence identifies the matching transcript as a side effect.
- **A supplied FASTA reproduces the API's canonicals, needs no network,
  and records which UniProt release was used.**
  - Mapped by sequence first and gene name second, Swiss-Prot gives a
    canonical for 272 of 274 protein-coding chr15 genes.
  - It agrees with the API's answer for 183 of 184 genes. The exception
    has only an unreviewed (TrEMBL) entry.
- **Building the renaming step into JCAST is straightforward.** With the
  API canonicals, 83 chr15 output records are identical to a Swiss-Prot
  isoform and 184 to a Swiss-Prot canonical. These are the records the R
  merge renames or drops.

## 2. How the canonical is obtained today

`Sequence.get_canonical_aa_uniprot` queries
`https://www.ebi.ac.uk/proteins/api/proteins/Ensembl:<gene id>` for the
first reviewed entry, and falls back to the first unreviewed one. It keeps
the answer in `cache/uniprot-cache.db` under the working directory
(REVIEW.md MED-5). Consequences:

- **A network call per gene.** A failed call (after 10 retries) leaves the
  canonical empty. Every slice of that gene is then an orphan, and the
  only sign is a warning in the log.
- **The UniProt release is not recorded.** Cached answers are reused
  indefinitely, whatever the release.
  - On these data the drift is small: 183 of the 184 cached chr15
    sequences equal release 2026_04.
  - The one that doesn't is the TrEMBL entry, absent from the Swiss-Prot
    file.
  - Two runs at different times can still use different canonicals
    without any record of it.
- **Isoforms are not used.** JCAST records identical to a Swiss-Prot
  isoform keep JCAST headers. The post-run R script
  (`combine_JCAST_comprehensive_canonical.R`) handles this:
  - reads `psq_T1.fasta` and drops `tr|` records;
  - renames each record to `sp|<accession>-J<n>|<entry name>`;
  - writes Swiss-Prot (canonical + isoforms) plus the JCAST records not
    identical to a Swiss-Prot entry, as one search database.

## 3. Results (chr15)

### 3.1 Mapping a Swiss-Prot entry to each gene

Of the 294 chr15 rMATS genes, 274 have a protein-coding transcript in the
GTF.

| Mapping | Unique | Several | None |
|---|---|---|---|
| By gene name (`GN=` = rMATS gene symbol) | 259 | 0 | 15 |
| By sequence (a Swiss-Prot canonical equals a protein translated from one of the gene's GTF transcripts) | 267 | 0 | 7 |
| **Sequence, then name** | **272** | 0 | 2 (JMJD7-PLA2G4B, a readthrough gene; AC012651.1, TrEMBL only) |

Against the 184 cached API answers:

| | Agrees with the API | Differs |
|---|---|---|
| Name mapping | 175 | 9: 8 genes renamed by HGNC since Ensembl 89 (e.g. FAM214A → ATOSA, FAM96A → CIAO2A, ADAL → MAPDA, TARSL2 → TARS3), 1 TrEMBL-only |
| Sequence mapping | 178 | 6: 5 canonical sequences differ from every Ensembl 89 translation, 1 TrEMBL-only |
| **Sequence, then name** | **183** | 1 (AC012651.1, TrEMBL-only) |

The GTF (Ensembl 89, 2017) and UniProt (2026) are nine years apart. With
matching releases, sequence mapping should do better still.

### 3.2 Coordinates for the canonical

179 of 184 API canonicals (97%) equal a protein translated from one of the
gene's GTF transcripts, so each of their residues has genomic coordinates.
Every GTF canonical has coordinates by construction.

### 3.3 Output records identical to a Swiss-Prot entry

The run uses the cached API canonicals, as a default run would; genes
without a cached answer have no canonical.

| File | Identical to a Swiss-Prot isoform | Identical to a Swiss-Prot canonical | Neither |
|---|---|---|---|
| `psq_T1` | 54 | 0 | 95 |
| `psq_T1_alt_*` (incl. orphans) | 6 | 1 | 23 |
| `psq_T2` | 1 | 0 | 23 |
| `psq_T4` (incl. orphans) | 22 | 1 | 343 |
| `psq_canonical` | 0 | 182 | 1 |
| other orphans | 0 | 0 | 114 |
| **All** | **83** | **184** | **599** |

In `psq_T1`, 54 of 149 records (36%) are already Swiss-Prot isoforms.
Built-in renaming would label them with the isoform accession (e.g.
`Q9XXXX-2`), instead of leaving the R step to drop them.

## 4. Joining by coordinates (question 1)

**What it needs:**
- **The canonical's coordinates.**
  - For a GTF canonical they are known.
  - For a UniProt canonical, the transcript whose translation equals it:
    97% of genes here.
  - Otherwise the sequence join is the fallback.
- **The slice's codon positions.** These are known (PR 4b/4c).
- **The join point:** the canonical residue whose codon starts, in frame,
  where the slice's first (last) codon starts. A codon is placed by its
  first base, which tolerates NAGNAG-type splice-site shifts.
- **The PR 5 rule that complete ends are not joined** still applies.

The audit's ground truth
([stitch audit §3](2026-11-11-stitch-anchor-audit.md#3-method)) is exactly
this.

**What it would gain on chr15:** after PR 5, the sequence join is correct
for 927 of 930 joinable attempts, with no wrong joins. Joining by
coordinates would add the 3 misses:
- CD276 A3SS 1316, both slices, unjoined by decision;
- NPTN A3SS 1317, whose shared end is shorter than 7 residues.

**Cost:**
- **A second join path,** chosen per gene by whether coordinates are
  available, plus a transcript lookup for UniProt canonicals.
- **Tests for both paths.**
- **About 100 lines.** The audit already contains the core functions
  (`coding_positions`, `true_anchor`).

**Assessment:** not worthwhile now; the gain is about 0.3% of joins. Two
things change that:
- A FASTA input mapped by sequence (§5) gives each canonical's transcript
  for free.
- The audit's coordinate truth can serve as a standing regression check.
  It already runs offline in about 4 minutes.

## 5. Recommendations

1. **Add `--uniprot FASTA`** (Swiss-Prot, with isoforms, as downloaded from
   UniProt).
   - **Choosing the canonical:** map each gene to a Swiss-Prot canonical by
     sequence first (an entry equal to one of the gene's GTF-translated
     proteins), then by a unique gene name.
   - **Gene left unmapped:** use the GTF canonical. Alternatively, leave
     it without one, which writes its slices as orphans as today; that is
     decision 3.
   - **The log:** record the file's name, entry count and checksum. The
     UniProt release is not in the FASTA itself, so the README should ask
     users to keep the download date or release.
   - **Header:** JCAST records keep the UniProt accession and entry name,
     as today.
2. **Build in the Swiss-Prot merge (with `--uniprot`).**
   - **New file `psq_isoforms.tsv`:** every JCAST record identical to a
     Swiss-Prot canonical or isoform, with that accession.
   - **New file `jcast_uniprot_merged.fasta`:**
     - all Swiss-Prot entries;
     - plus the JCAST records from the chosen tiers that are not identical
       to any Swiss-Prot entry;
     - each renamed `sp|<accession>-J<n>|<entry name> <JCAST header>`, as
       in the R script but without its `_HUMAN` assumption.
   - **Unchanged:** the individual `psq_*.fasta` files.
3. **Keep the API path when `--uniprot` is not given,** for compatibility.
   Document `--uniprot` as the recommended mode: offline, reproducible,
   faster.
4. **Defer joining by coordinates** (§4). Revisit once 1 is in, since the
   sequence mapping provides the transcripts it needs.

## 6. Decisions requested

1. Proceed with `--uniprot FASTA` and the built-in merge (§5.1–5.2) as the
   next PR?
2. Without `--uniprot`, keep the API path (recommended), or require the
   FASTA in v0.4.0?
3. A gene with no Swiss-Prot entry (2 of 274 here): fall back to the GTF
   canonical (recommended), or leave it without a canonical as today?
4. The merged database: which tiers?
   - The R script uses `psq_T1` only.
   - Recommended default: `psq_T1` and `psq_T1_alt_*`, both complete
     isoforms of annotated structure.
   - A `--merge_tiers` option could widen it.
   - Include TrEMBL (`tr|`) records? The R script drops them; recommended
     to drop them too.
5. Joining by coordinates: defer (recommended)?

## 7. Decisions (2026-11-11)

1. **Add `--uniprot FASTA` and the built-in merge (§5.1–5.2), as PR 14.**
   - Considered instead: naming a proteome (e.g. `UP000005640`) and having
     JCAST download it. Rejected for now because it ties a run to the
     download date. The FASTA leaves the release choice and its record
     with the user.
2. **Keep the API path when `--uniprot` is not given.**
3. **A gene with no Swiss-Prot entry falls back to the GTF canonical** (with
   `--uniprot` only; the API path is unchanged).
4. **The merged database contains `psq_T1` and `psq_T1_alt_*`**, with
   TrEMBL (`tr|`) records dropped. A later option may add T2 and T4.
5. **Joining by coordinates is deferred.**
