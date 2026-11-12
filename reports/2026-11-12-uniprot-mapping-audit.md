# Audit: how `--uniprot` maps genes to UniProt entries, against the UniProt API

| | |
|---|---|
| **Date** | 2026-11-12 |
| **Status** | Audit. Options in §5; decisions requested in §6 and approved 2026-11-12 (§7). No code changed. Measured at `8b845aa`. |
| **Question** | Under `--uniprot`, where does the Ensembl-to-UniProt mapping come from, and how complete is it next to the UniProt API that JCAST uses without it? Raised with the decisions on the [PR 6–8 audit](2026-11-12-thresholds-runs-canonical-audit.md#7-decisions-2026-11-12). Extends the chr15 comparison of the [UniProt FASTA evaluation §3.1](2026-11-11-uniprot-fasta-evaluation.md#31-mapping-a-swiss-prot-entry-to-each-gene) to a whole genome. |
| **Reproduce** | `python reports/mapping_audit/audit_mapping.py --jcast_out <run> --uniprot <its FASTA> --gtf_cached <gtf>.cached --xrefs <file>.tsv.gz --gtf <gtf> --genome <fa>` (≈ 5 min; downloads the cross-references if the file is missing) |
| **Inputs** | <ul><li>The user's whole-genome run of 2026-11-12: `-m -c --uniprot -p 8`, 17,068 genes, GENCODE v50, Swiss-Prot human with isoforms 2026_04.</li><li>UniProt's Ensembl cross-references for every human UniProtKB entry that has one (262,883 entries, 23,750 genes), release 2026_04, from `rest.uniprot.org`. These are what the API looks a gene up by.</li><li>The local API cache: the API's actual answers for 6,519 of these genes, from an earlier run.</li></ul> None is bundled. |
| **Evidence** | [mapping_audit/genome/summary.txt](mapping_audit/genome/summary.txt) · [genes.tsv](mapping_audit/genome/genes.tsv) (one row per gene) |

## 1. Summary

- **`--uniprot` maps no identifiers.**
  - GENCODE and Ensembl GTFs carry no UniProt accessions, and the UniProt
    FASTA carries no Ensembl ids.
  - JCAST assigns each gene, in this order
    ([jcast/uniprot.py](../jcast/uniprot.py)):
    1. the UniProt canonical whose **sequence** equals a protein
       translated from one of the gene's GTF transcripts;
    2. else the entry whose **gene name** (`GN=`) is the gene's symbol;
    3. else the GTF canonical.

    Each step needs exactly one match, and Swiss-Prot entries come before
    TrEMBL ones.
  - The API instead looks the gene up by **UniProt's own Ensembl
    cross-references,** and takes the first entry returned.
- **The two agree for 99.35% of the genes that matter.**
  - Of the 13,365 protein-coding genes that UniProt cross-references to a
    Swiss-Prot entry, `--uniprot` picks the same entry for 13,278: 13,066
    by sequence and 212 by name.
  - Against the API's actual cached answers (6,519 genes), 6,374 agree
    (97.8%). Most of the rest are genes with only TrEMBL entries, which a
    Swiss-Prot FASTA cannot provide.
- **Where they differ, `--uniprot` is mostly the more defensible of the
  two:**
  - **86 genes:** `--uniprot` declines a match where the API picks one,
    often arbitrarily. These are readthrough genes such as NME1-NME2,
    which the API maps to NME2; loci with several Swiss-Prot entries
    (GNAS, CDKN2A); and unnamed novel genes.
  - **55 protein-coding genes:** `--uniprot` finds the Swiss-Prot entry
    where UniProt cross-references only TrEMBL (POLR2A, SHANK3, NLRP3,
    C1R).
- **The real gap is genes with no canonical at all.**
  - 384 protein-coding genes get neither a Swiss-Prot match nor a GTF
    canonical, and their 361 records are all orphans. That is 2.6% of the
    run's 13,718 orphans, and 0.4% of all 86,516 records.
  - The API would give 351 of them a TrEMBL entry, and 33 a Swiss-Prot
    entry.
  - **Remedies:**
    - A FASTA that includes TrEMBL would give 203 of the 384 a canonical,
      covering 110 of those records.
    - Choosing the GTF canonical by GENCODE's `Ensembl_canonical` tag
      when no transcript passes the current filter would cover 173 of the
      384.

## 2. The two mappings

| | `--uniprot FASTA` | UniProt API (no `--uniprot`) |
|---|---|---|
| Key | Protein sequence, then gene name | The gene's Ensembl id, through UniProt's cross-references |
| Ambiguity | Needs exactly one match; otherwise the next step | Takes the first entry returned (`size=1`) |
| TrEMBL | Only if the FASTA contains it | Yes, when the gene has no Swiss-Prot entry |
| Release | Whatever was downloaded (checksum logged) | Whatever the API served at the time, cached indefinitely |
| Gene without a match | The GTF canonical (longest protein-coding transcript with TSL ≤ 2 and annotated start and stop codons), else none | None |

Cross-references and the FASTA here are the same release (2026_04).
UniProt builds its Ensembl cross-references largely by matching
sequences, so the two methods are less different than they look.

## 3. Results

**Protein-coding genes (13,823)**, by the run's choice and UniProt's
cross-references for the gene. Records are the run's output records of
those genes ([summary.txt](mapping_audit/genome/summary.txt)):

| Genes | `--uniprot` | Cross-references | Records (joined / orphan) |
|---|---|---|---|
| 13,066 | Swiss-Prot by sequence | the same entry | 70,827 / 12,834 |
| 212 | Swiss-Prot by name | the same entry | 1,469 / 357 |
| 1 | Swiss-Prot by name | another entry (PINX1, §4) | 0 / 0 |
| 53 | GTF canonical | Swiss-Prot | 228 / 68 |
| 33 | none | Swiss-Prot | 0 / 23 |
| 55 | Swiss-Prot (49 by name, 6 by sequence) | TrEMBL only | 158 / 78 |
| 4 | Swiss-Prot by name | none | 0 / 0 |
| 48 | GTF canonical | TrEMBL only, **same sequence** as the GTF canonical | 116 / 20 |
| 351 | none | TrEMBL only | 0 / 338 |

**Other genes (3,245):**
- **3,061** have nothing on either side.
- **183** lncRNA and pseudogene loci get a Swiss-Prot entry by name from
  `--uniprot`, with no cross-reference. Examples are PMS2P3, IGF2-AS and
  MIR9-1HG (Protein CROC-4), each a named Swiss-Prot entry. They have no coding transcripts, so they produce no
  records. Their canonical goes to `psq_canonical` under `-c`.
- **1** has a TrEMBL cross-reference only.

**Against the API's cached answers** (6,519 genes from an earlier API run
of the same data):

| `--uniprot` | API | Genes |
|---|---|---|
| same entry | | 6,358 |
| another entry with the same sequence | | 16 |
| none or GTF | TrEMBL | 98 |
| GTF or none | Swiss-Prot | 27 |
| Swiss-Prot | TrEMBL | 17 |
| another Swiss-Prot entry | Swiss-Prot | 3 |

## 4. Where they differ, and why

**86 protein-coding genes cross-referenced to Swiss-Prot, where
`--uniprot` uses the GTF or has nothing:**
- **39 readthrough genes** (NME1-NME2, STON1-GTF2A1L, CORO7-PAM16,
  LY75-CD302, …).
  - UniProt cross-references the readthrough gene to one component's
    entry, because they share transcripts.
  - The API would make, for example, NME2 the canonical of NME1-NME2.
  - `--uniprot` finds no Swiss-Prot sequence equal to the readthrough
    protein, and no entry named NME1-NME2, so it uses the readthrough
    gene's own GTF protein.
- **21 genes cross-referenced to several Swiss-Prot entries:**
  - GNAS (3 entries);
  - CDKN2A (p16 and ARF, from alternative reading frames);
  - MOCS2 (A and B);
  - NRXN1 and NRXN2 (α and β);
  - TMPO, BBC3, ERCC6, GPR89A and GPR89B (identical proteins), …

  `--uniprot` needs a unique match, so it falls back. The API returns
  whichever entry comes first.
- **25 genes without a name,** whose rMATS symbol is the gene id. These
  are mostly novel loci overlapping a named gene, which UniProt
  cross-references to that gene's entry. Examples: ENSG00000285130 →
  FXN, ENSG00000266086 → SRSF1. With no name to match and no identical
  sequence, `--uniprot` falls back.
- **1:** MACF1. No Swiss-Prot sequence equals one of its GENCODE
  translations, and its name matches two entries (isoforms 1–5 and
  isoforms 6/7).

**55 protein-coding genes where `--uniprot` finds a Swiss-Prot entry and
UniProt cross-references only TrEMBL** (POLR2A, SHANK3, NLRP3, C1R, MICA,
…). For these genes, UniProt links the Ensembl translations to TrEMBL
entries, presumably because the reference genome's translations differ
from the Swiss-Prot sequence. For the 49 matched by name, no single
Swiss-Prot sequence equals one of their GENCODE translations, which is
why the sequence step did not match.
- The API would return that TrEMBL entry.
- `--uniprot` returns the Swiss-Prot entry: by name for 49 genes, by
  sequence (from another transcript) for 6.
- For proteomics, the reviewed entry is the better canonical.

**PINX1:** GENCODE v50 has two genes named PINX1. The second,
ENSG00000258724, a locus that UniProt links to SOX7, also took PINX1's
entry by name. It has no records. Gene names are not unique in GENCODE,
so mapping by name can assign another gene's protein.

**384 protein-coding genes with no canonical:**
- **Why the GTF gives none:**
  - 211 have no protein-coding transcript at all in GENCODE v50, despite
    the gene type;
  - 149 have only transcripts of TSL 3–5 or NA;
  - 24 have TSL ≤ 2 transcripts without annotated start and stop codons.
- **With TrEMBL entries in the FASTA as well,** the same rules would give
  203 of them a canonical: 199 by sequence and 4 by name. That covers
  110 of their 361 orphan records. The other 181 still have none.
- **173 of the 384 have a protein-coding transcript tagged
  `Ensembl_canonical`.** That includes 86 of the 181 left by TrEMBL.

## 5. Options

1. **Keep the mapping as it is** (recommended). It agrees with UniProt's
   own cross-references for 99.35% of the coding genes that have a
   Swiss-Prot entry. Most of its differences are matches the API makes
   arbitrarily, or Swiss-Prot entries the API misses.
2. **Document a Swiss-Prot + TrEMBL FASTA** for users who want the
   fullest coverage.
   - 203 more genes would get a canonical, covering at most 110 more
     joinable records.
   - Records joined to a TrEMBL canonical are named `tr|…`, and are left
     out of `jcast_uniprot_merged.fasta` by design.
   - The human TrEMBL set is large: about ten times Swiss-Prot.
   - This needs only documentation; the code already handles TrEMBL.
3. **The GTF canonical from GENCODE's `Ensembl_canonical` tag** when no
   transcript passes the current filter. This would cover 173 of the
   384, without TrEMBL. It is the deferred "tag-based canonical
   selection" item (REVIEW.md §9, deferred to v0.5.x). It changes output
   in GTF-only and `--uniprot` runs, so it needs its own audit.
4. **Map by name only when the name is unique in the GTF.** This fixes
   PINX1-like cases: 1 gene, no records. Small.
5. **Not recommended:** a cross-reference step, with UniProt's Ensembl
   cross-references as an extra input. Its additional answers are mostly
   the readthrough and multi-entry assignments in §4.

## 6. Decisions requested

1. Keep the mapping (recommended)?
2. Document the Swiss-Prot + TrEMBL option in the README (recommended)?
3. Schedule the `Ensembl_canonical` fallback for the GTF canonical: in
   v0.4.0 after PR 16, or leave it in v0.5.x?
4. Map by name only for gene names unique in the GTF (recommended; small,
   could join PR 8)?

## 7. Decisions (2026-11-12)

1. **Keep the mapping** as it is. Approved.
2. **Document the Swiss-Prot + TrEMBL option** in the README. Approved,
   with PR 8.
   - For human runs, TrEMBL is not expected to be used, and the 203 genes
     it would add are accepted as missing.
   - For other organisms, a FASTA with TrEMBL may be the default.
3. **The `Ensembl_canonical` fallback for the GTF canonical** is scheduled
   in v0.4.0, after PR 16, as PR 17. It changes output, so it gets its
   own audit first.
4. **Map by gene name only when the name is unique in the GTF.**
   Approved, with PR 8.

