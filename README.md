# JCAST — Junction-Centric Alternative Splicing Translator

JCAST translates rMATS splice junctions into custom protein FASTA databases
for shotgun mass-spectrometry-based proteomics. By restricting the database to
the splice isoforms that are likely present in a specific RNA-seq sample,
JCAST avoids the database-size inflation and false-positive identification
issues that occur when searching against a full transcriptome translation.

- Documentation: <https://ed-lau.github.io/jcast/>
- Source: <https://github.com/ed-lau/jcast>
- Changelog: [CHANGELOG.md](CHANGELOG.md)
- Developer review & v0.4.0 roadmap: [REVIEW.md](REVIEW.md)

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Command-line options](#command-line-options)
- [Output files](#output-files)
- [Worked example: ENCODE heart](#worked-example-encode-heart)
- [Known issues](#known-issues)
- [FAQ](#faq)
- [Citations](#citations)
- [Contributors](#contributors)
- [License](#license)

## Installation

JCAST requires Python 3.7+. (For v0.4.0 the floor will move to 3.10; see
[REVIEW.md](REVIEW.md).) We recommend installing into a virtual environment:

```bash
pip install jcast
```

To install the optional Qt-based GUI (PySide6):

```bash
pip install 'jcast[gui]'
```

To install from source for development:

```bash
git clone https://github.com/ed-lau/jcast
cd jcast
pip install -e '.[gui]'
```

## Quick start

After installing, launch JCAST as a module to see usage:

```bash
python -m jcast
# or, equivalently,
jcast
```

A typical analysis run:

```bash
python -m jcast \
    /data/encode_human_pancreas/rmats \
    /data/gtf/Homo_sapiens.GRCh38.89.gtf \
    /data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -o encode_human_pancreas \
    -q 0 1 \
    -r 1 -m -c
```

To launch the (optional) graphical interface:

```bash
python -m jcast --gui    # or: jcast-gui
```

To run the bundled smoke test (rMATS sample + chromosome-15 genome):

```bash
pip install tox
tox
```

To exercise JCAST end-to-end on the bundled fixtures:

```bash
python -m jcast \
    tests/data/rmats \
    tests/data/genome/Homo_sapiens.GRCh38.89.chromosome.15.gtf \
    tests/data/genome/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz \
    -o ~/Desktop/jcast_smoke
```

## Command-line options

```text
usage: jcast [-h] [-o OUT] [-r READ] [-m] [-c] [--mask MASK]
             [-s {MXE,SE,RI,A5SS,A3SS} [...]] [-q q_lo q_hi]
             [--g_or_ln G_OR_LN]
             rmats_folder gtf_file genome

positional arguments:
  rmats_folder          path to folder storing rMATS output
  gtf_file              path to Ensembl/Gencode GTF file
  genome                path to genome FASTA (optionally .gz)

optional arguments:
  -h, --help            show this help message and exit
  -o, --out OUT         output directory [default: out]
  -r, --read READ       minimum skipped-junction read count for a junction to
                        be translated [default: 1]
  -m, --model           model the read-count cutoff via a mixture model
                        [default: False]
  -c, --canonical       write canonical protein sequences even when a slice
                        is untranslatable [default: False]
  --mask MASK           masked-nucleotide characters; slices containing any
                        of them are skipped [default: NnXx]
  -s, --splice_type     splice types to translate [default: MXE SE RI A5SS A3SS]
  -q, --qvalue q_lo q_hi
                        keep junctions whose rMATS FDR is in this range
                        [default: 0 1]
  --g_or_ln G_OR_LN     experimental: low-end distribution selector
                        (0 = Gamma, otherwise LogNormal)
```

The `--g_or_ln` flag is a development-only switch and will be removed in
v0.4.0; see [REVIEW.md](REVIEW.md).

## Output files

JCAST writes a small number of FASTA files into the output directory. Each
entry's header encodes the source rMATS row, slice number, junction
coordinates, translation strand, frame, and tier (see below).

| File | Contents |
|---|---|
| `psq_canonical.fasta` | Sequences identical to the canonical SwissProt/GTF protein for the gene. Headers follow UniProt convention. |
| `psq_T1.fasta` | **Tier 1.** Both slices translate in their annotated reading frames with no premature stop codon (PTC), and the event does not shift the reading frame: the slices' coding lengths, measured from each slice's frame, differ by a multiple of 3. Both slices are stitched into the canonical full-length protein. |
| `psq_T2.fasta` | **Tier 2.** As T1, but the coding lengths differ by a non-multiple of 3 (frameshift). A slice that cannot be stitched at both ends may be joined at its N-terminus only. |
| `psq_T4.fasta` | **Tier 4.** Only one slice translates; the other hits a PTC. The PTC slice is translated *up to* its stop codon, giving a truncated product, which is joined to the canonical at its N-terminus. If it is shorter than `params.ptc_threshold` × the other slice, or cannot be joined, it goes to `psq_T4_orphan.fasta`. **Treat T4 with caution — exclude from MS searches unless you have a specific reason to include it.** |
| `psq_T#_orphan.fasta` | Slices of that tier that translated but could not be stitched back to the canonical sequence. **Treat orphans with caution.** |

There is no Tier 3 from v0.4.0. It re-translated in a guessed frame when the
annotated frame failed, and was removed together with the three-frame
fallback (see [CHANGELOG.md](CHANGELOG.md)).

How a junction is processed (v0.4.0):

1. **Trimming.** Every exon is trimmed to the gene's annotated coding region.
   This runs from the most upstream start codon to the most distal stop
   codon of protein-coding transcripts with transcript support level ≤
   `params.tsl_threshold`. Exons entirely in the UTR are removed.
2. **Skip if outside the coding region.** If the splice event lies entirely
   outside the coding region, both slices are identical after trimming. The
   junction is skipped (`SKIPPED. Splice event lies outside the annotated
   coding region.`), and the canonical is still written with `-c`.
3. **Reading frame.** Each slice is read in its own frame, taken from the
   annotated coding sequence at the slice's 5′ end. A slice with no annotated
   frame is not translated, for example in a lncRNA or a pseudogene. Frames
   are never guessed.

Skipped and failed junctions write no isoform sequences. Each junction's
outcome (fate) is written to the log.

For a typical proteomics workflow, append `psq_T1.fasta` and
`psq_canonical.fasta` to a base SwissProt database and search as usual.
Including T2, T4 and orphans is possible but increases false-discovery risk.

### Header format

Non-canonical FASTA entries are named as follows:

```
>sp|Q91VW5|GOGA4_MOUSE|ENSMUSG00000038708|MXE1|0|chr9|118560742:118560872|118565557:118565667|+2|r521|T1 sp|Q91VW5|...
```

The `|`-separated fields are:

1. Knowledgebase code (`sp` for SwissProt).
2. UniProt accession of the canonical entry.
3. UniProt entry name.
4. Ensembl gene id.
5. rMATS splice type and slice number (e.g. `MXE1`, `MXE2`).
6. rMATS row id.
7. Chromosome.
8. Anchor exon `start:end`.
9. Alternative exon `start:end`.
10. Translated strand and that slice's phase (e.g. `+2`).
11. Skipped-junction read count, prefixed with `r`.
12. Tier (`T1` / `T2` / `T4`).

## Worked example: ENCODE heart

A representative full pipeline from FASTQ to JCAST database:

```bash
# 1. Align with STAR
STAR --runThreadN 10 \
     --genomeDir /ref/GRCh38/STARindex \
     --sjdbGTFfile /ref/Homo_sapiens.gtf \
     --sjdbOverhang 100 \
     --readFilesIn /raw/heart_b1t1_R1.fastq.gz /raw/heart_b1t1_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix /work/star/b1t1/

# 2. Identify splice junctions with rMATS-turbo (in Docker, as one example)
sudo docker run -v /work:/data rmats:turbo01 \
     --b1 /data/b1.txt --b2 /data/b2.txt \
     --gtf /data/GRCh38.gtf --od /data/output \
     -t paired --nthread 4 --readLength 101 --anchorLength 1

# 3. Translate splice events with JCAST
python -m jcast /work/output \
                /ref/Homo_sapiens.gtf \
                /ref/GRCh38.primary_assembly.fa \
                -o /work/jcast_heart -m -c
```

## Known issues

- **Soft-masked (lowercase) genome sequences cause translation failure.** Use
  an unmasked or hard-masked genome until this is fixed. (Tracked for v0.4.0;
  see [REVIEW.md](REVIEW.md) HIGH-2/HIGH-3.)
- **rMATS rows with `NA` gene names can fail.** A defensive substitution was
  added in v0.3.5a; further fixes are planned.
- **Tier 4 contains products truncated at a stop codon.** Some of these are
  genuine alternative C-termini rather than premature stops; separating them
  is planned for v0.4.0 (REVIEW.md, PR 4c). Treat T4 outputs as candidates
  for follow-up validation, not as confirmed identifications.
- **UniProt and GTF `.cached` files are written to the current working
  directory** and the directory of the input GTF respectively — this fails on
  read-only volumes. Workaround: `cd` to a writable directory before running.
  Permanent fix is part of the v0.4.0 hardening release.

## FAQ

**What do I do with the output files?**
For a typical analysis, append `psq_T1.fasta` to either `psq_canonical.fasta`
(if you used `-c`) or your existing SwissProt canonical database, and search
that combined file in your MS pipeline. Lower-tier (T2/T4) and orphan
outputs contain frameshifted or truncated sequences and should be
interpreted with caution.

**What files are required in the rMATS directory input?**
JCAST expects an unaltered rMATS output folder with the five per-splice-type
files named exactly `MXE.MATS.JC.txt`, `SE.MATS.JC.txt`, `RI.MATS.JC.txt`,
`A3SS.MATS.JC.txt`, `A5SS.MATS.JC.txt`.

**Can I use a Gencode GTF?**
Yes. JCAST silently aliases `transcript_type` to `transcript_biotype` when
reading a Gencode-style GTF.

**Why does my first run write a `.cached` file next to my GTF?**
JCAST caches the parsed GTF as a TSV to speed re-runs. The cache lives next
to the input file. If your reference directory is read-only, work around by
copying the GTF to a writable location. (Targeted fix in v0.4.0.)

**How can I contribute?**
Open an issue or a pull request on
<https://github.com/ed-lau/jcast>, or reach out to
[edward.lau@cuanschutz.edu](mailto:edward.lau@cuanschutz.edu). The
[REVIEW.md](REVIEW.md) document enumerates the current backlog and is a good
starting point for new contributors.

## Citations

If you use JCAST in published work, please cite:

1. **JCAST: Sample-Specific Protein Isoform Databases for Mass
   Spectrometry-based Proteomics Experiments.**
   RW Ludwig, E Lau. *Software Impacts* 10, 100163.
   [Link](https://www.sciencedirect.com/science/article/pii/S2665963821000695)
2. **Determining Alternative Protein Isoform Expression Using RNA Sequencing
   and Mass Spectrometry.**
   Y Han, JM Wright, E Lau, MPY Lam. *STAR Protocols* 1 (3), 100138.
   [Link](https://pubmed.ncbi.nlm.nih.gov/33377032/)
3. **Splice-Junction-based Mapping of Alternative Isoforms in the Human
   Proteome.**
   E Lau, Y Han, DR Williams, CT Thomas, R Shrestha, JC Wu, MPY Lam.
   *Cell Reports* 29 (11), 3751-3765.e5.
   [Link](https://www.sciencedirect.com/science/article/pii/S2211124719314974)

## Contributors

- **Edward Lau, PhD** — code & design — [ed-lau](https://github.com/ed-lau)
- **Maggie Lam, PhD** — code & design — [Maggie-Lam](https://maggielamlab.org)
- **Robert Wes Ludwig, BSc** — modeling — [WesLudwig](https://github.com/WesLudwig)

## License

MIT — see [LICENSE.md](LICENSE.md).
