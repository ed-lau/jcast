# UniProt test data

`uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz`: 711 entries from
UniProtKB/Swiss-Prot human with isoforms, release 2026_04 (4 Nov 2026). It
holds the 275 canonical entries that the chr15 test genes map to with
`--uniprot`, or that the UniProt API returned for them, plus all their
isoforms.

Downloaded from
`https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=(proteome:UP000005640)+AND+(reviewed:true)`
(42,545 entries) and subset with
`reports/uniprot_eval/genes.tsv`; see
[the UniProt FASTA evaluation](../../../reports/2026-11-11-uniprot-fasta-evaluation.md).

UniProt data are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). © The UniProt
Consortium.
