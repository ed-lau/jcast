# -*- coding: utf-8 -*-

""" Parameters """

# transcript filtering
tsl_threshold = 1  # the transcript levels below which (lower is better) to consider

# stop codon forced translation
ptc_threshold = 0.33   # the PTC slice should be at least this portion of the long slice to be included

# canonical stitching
max_retries = 10  # max number of retries to Uniprot
stitch_length = 9  # amino acid joint to stitch to canonical
