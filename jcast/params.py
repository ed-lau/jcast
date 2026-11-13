# -*- coding: utf-8 -*-

""" Experimental parameters """

"""
Annotated transcript filtering

Setting TSL threshold to 1 excludes some Uniprot canonical transcripts, e.g.,
DDP8_HUMAN with the first 18 amino acids.

"""
tsl_threshold = 2  # the transcript levels below which (lower is better) to consider


"""
Stop codon forced translation

This threshold may be set to allow some stop codon transcripts to be translated to Tier 4

"""
ptc_threshold = 0.33  # the PTC slice should be at least this portion of the long slice to be included


"""
Canonical transcript recognition

Canonical can be chosen from the GTF (longest protein coding sequence) alone 
which speeds up the analysis or through an API call to Uniprot. 
This has the advantage of allowing splice graphs to determine which junctions
map to which transcripts.
However, note that Uniprot has manual annotation on what is the canonical sequence
through prevalence and homology to other species that are not apparent in the GTF file without
calling to other resources like APPRIS (e.g. see DET1_HUMAN). 
Because of the prominence of Uniprot in proteomics work
we have chosen to use Uniprot for now.

2026-11-11: with --uniprot, canonicals come from a supplied UniProt FASTA instead of the API (by sequence, then gene
name), and genes without an entry fall back to the GTF; use_gtf_only is then ignored. See
reports/2026-11-11-uniprot-fasta-evaluation.md

"""
use_gtf_only = False  # use GTF but not Uniprot to get canonical transcripts
uniprot_max_retries = 10  # max number of retries if retrieving sequences from Uniprot
uniprot_timeout = 30  # seconds to wait for each UniProt API request (connect, and between bytes read)
aa_stitch_length = 7  # anchor length (residues) for joining a translated slice to the canonical
aa_terminus_length = 10  # residues compared at each end to decide whether a product's terminus differs from the canonical

"""
2026-11-11 the anchor length was 10. Seven is the shortest length at which, across the MANE Select proteome, anchors
repeated within their own protein and chance matches from other proteins level off; longer anchors only miss slices
with short shared ends. See reports/2026-11-11-stitch-anchor-audit.md
"""
