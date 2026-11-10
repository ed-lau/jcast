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
Three-frame fallback (PR 4b audit switch; to be removed or kept according to that audit)

When a slice has no annotated reading frame, or does not translate in it, try all three frames (Tier 3).
"""
three_frame_fallback = True

# PR 4b audit switch: read a frame for a slice that lacks the gene's start codon when the other slice contains it
frame_for_start_less_slice = True


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

"""
use_gtf_only = False  # use GTF but not Uniprot to get canonical transcripts
uniprot_max_retries = 10  # max number of retries if retrieving sequences from Uniprot
aa_stitch_length = 10  # amino acid joint to stitch to canonical when using aa for joining
