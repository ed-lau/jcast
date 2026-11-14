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

# Each fate's short name in junctions.tsv and the end-of-run summary, in summary order
# (reports/2026-11-12-log-audit.md §4)
NAMES = [(TIER1, 'T1'),
         (TIER1_ALT_TERMINUS, 'T1_alt_terminus'),
         (TIER2, 'T2'),
         (TIER4, 'T4'),
         (FAIL, 'fail'),
         (EXCEPTION, 'exception'),
         (SKIPPED_LOW_COUNT, 'skipped_low_count'),
         (SKIPPED_P_VALUE, 'skipped_fdr'),
         (SKIPPED_MASKED, 'skipped_masked'),
         (SKIPPED_OUTSIDE_CDS, 'skipped_outside_cds'),
         ]


def filtered(min_jc, fdr, cutoff: int, qvalue) -> str:
    """
    The fate of a junction that the read count and FDR filters skip, or None if they keep it: it is translated. The
    run applies it to each junction, and the GUI's preview of the -m model to every junction at once
    (reports/2026-11-13-gui-cli-scope.md §4)

    :param min_jc: the junction's read count (Junction.min_jc)
    :param fdr: its rMATS FDR
    :param cutoff: the read count cutoff (-r, or the -m model's): a junction with fewer reads is skipped
    :param qvalue: (q_lo, q_hi) (-q): a junction whose FDR is outside is skipped
    :return: SKIPPED_LOW_COUNT, SKIPPED_P_VALUE or None
    """
    q_lo, q_hi = qvalue
    if min_jc < cutoff:
        return SKIPPED_LOW_COUNT
    if not q_lo <= fdr <= q_hi:
        return SKIPPED_P_VALUE
    return None


def name(fate: str) -> tuple:
    """
    The fate's short name, and the detail not recorded elsewhere in junctions.tsv: the slice that hit a premature
    termination codon (T4), or the exception (exception)

    :param fate: a fate message
    :return: (short name, detail)
    """
    for message, short in NAMES:
        head = message.split('{')[0]
        if fate.startswith(head):
            if message == TIER4:
                return short, 'slice {0}'.format(fate[len(head):].split(' ')[0])
            if message == EXCEPTION:
                detail = fate[len(head):].rstrip(' \n')
                return short, detail[:-1] if detail.endswith('.') else detail
            return short, ''
    raise ValueError(f'unknown fate: {fate!r}')
