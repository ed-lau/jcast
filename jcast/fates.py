# -*- coding: utf-8 -*-

""" Junction callback messages for different tiers """

SKIPPED_LOW_COUNT = 'SKIPPED. Sequence discarded due to low coverage. \n\n'
SKIPPED_P_VALUE = 'SKIPPED. Sequence discarded due to difference across replicates. \n\n'
SKIPPED_MASKED = 'SKIPPED. Sequence discarded due to masked or ambiguous bases. \n\n'
SKIPPED_OUTSIDE_CDS = 'SKIPPED. Splice event lies outside the annotated coding region. \n\n'

TIER1 = 'SUCCESS TIER 1. Retrieved phase: {0}. Used phase: {1}. No frameshift. \n\n'
TIER1_ALT_TERMINUS = 'SUCCESS TIER 1 ALT TERMINUS. Retrieved phase: {0}. Used phase: {1}. Files: {2}. \n\n'
TIER2 = 'SUCCESS TIER 2. Retrieved phase: {0}. Used phase: {1}. Frameshift. \n\n'
TIER4 = 'PARTIAL TIER 4. Slice {0} hit a premature termination codon. Translated fragment. \n\n'

FAIL = 'FAILURE.  No alternative translation. \n\n'

# an exception while translating the junction: it writes nothing, and the run continues
# (reports/2026-11-12-thresholds-runs-canonical-audit.md §3.3)
EXCEPTION = 'FAILURE (exception). {0}: {1}. \n\n'
