# -*- coding: utf-8 -*-

""" Junction callback messages for different tiers """

SKIPPED_LOW_COUNT = 'SKIPPED. Sequence discarded due to low coverage. \n\n'
SKIPPED_P_VALUE = 'SKIPPED. Sequence discarded due to difference across replicates. \n\n'
SKIPPED_MASKED = 'SKIPPED. Sequence discarded due to masked regions. \n\n'

TIER1 = 'SUCCESS TIER 1. Retrieved phase: {0}. Used phase: {1}. No frameshift. \n\n'
TIER2 = 'SUCCESS TIER 2. Retrieved phase: {0}. Used phase: {1}. Frameshift. \n\n'
TIER3 = 'SUCCESS TIER 3. GTF phase mismatch. Retrieved phase: {0}. Used phase: {1}. \n\n'
TIER4 = 'PARTIAL TIER 4. Slice {0} hit a premature termination codon. Translated fragment. \n\n'

FAIL = 'FAILURE.  No alternative translation. \n\n'