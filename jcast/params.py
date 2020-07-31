# -*- coding: utf-8 -*-

""" Parameters """

"""
Annotated transcript filtering
Setting TSL threshold to 1 excludes some Uniprot canonical transcripts, e.g.,
Human DDP8 with the first 18 amino acids.

"""
tsl_threshold = 2  # the transcript levels below which (lower is better) to consider


"""
Stop codon forced translation
This threshold may be set to allow some stop codon transcripts to be translated to Tier 4

"""
ptc_threshold = 0.33   # the PTC slice should be at least this portion of the long slice to be included


"""
Canonical transcript recognition

"""
use_gtf_only = True  # use GTF but not Uniprot to get canonical transcripts
uniprot_max_retries = 10  # max number of retries if retrieving sequences from Uniprot
aa_stitch_length = 9  # amino acid joint to stitch to canonical when using aa for joining
