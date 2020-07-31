# -*- coding: utf-8 -*-

""" Junction callback messages for different tiers """

skipped_low = 'SKIPPED. Sequence discarded due to low coverage. \n\n'
skipped_p = 'SKIPPED. Sequence discarded due to difference across replicates. \n\n'

tier1 = 'SUCCESS TIER 1. Retrieved phase: {0}. Used phase: {1}. No frameshift. \n\n'
tier2 = 'SUCCESS TIER 2. Retrieved phase: {0}. Used phase: {1}. Frameshift. \n\n'
tier3 = 'SUCCESS TIER 3. GTF phase mismatch. Retrieved phase: {0}. Used phase: {1}. \n\n'
tier4 = 'PARTIAL TIER 4. Slice {0} hit a premature termination codon. Translated fragment. \n\n'

fail = 'FAILURE.  No alternative translation. \n\n'