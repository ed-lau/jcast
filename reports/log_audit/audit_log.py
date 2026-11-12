# -*- coding: utf-8 -*-

"""Audit of the run log's size and content (reports/2026-11-12-log-audit.md).

Reads a JCAST logfile.log in one pass and counts lines and bytes per kind of message. A message's bytes include its
line prefix (time stamp and level), the `~module:…~function:…` suffix the formatter appends, and the blank lines
that fate messages end with. Writes into --out:

    kinds.tsv     per kind of message: messages, bytes, share of the file, bytes per junction
    summary.txt   the totals quoted in the report

Usage (from the repository root):
    python reports/log_audit/audit_log.py --log <output dir>/logfile.log
"""

import argparse
import collections
import re
from pathlib import Path

# message kinds, by the start of the message; the first match wins
KINDS = [
    ('junction start (>>>>>> Doing …)', '>>>>>> Doing '),
    ('coding region (Transcription start … end …)', 'Transcription start: '),
    ('slice frames (Retrieved phases)', 'Retrieved phases: '),
    ('slice nucleotides (Retrieved nucleotide …)', 'Retrieved nucleotide for '),
    ('slice peptides (Translated AA …)', 'Translated AA for '),
    ('PTC slice peptide (Force-Translated AA …)', 'Force-Translated AA for '),
    ('orphan notice (Slice … does not stitch …)', 'Slice '),
    ('canonical chosen (Canonical for …)', 'Canonical for '),
    ('fate', ('SKIPPED', 'SUCCESS', 'PARTIAL', 'FAILURE')),
]
LINE = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}~(\w+)~')


def kind_of(message):
    for name, start in KINDS:
        if message.startswith(start):
            return name
    return 'other'


def main_(a):
    lines, size = collections.Counter(), collections.Counter()
    overhead = 0
    junctions = 0
    current = None
    with open(a.log, 'rb') as fh:
        for raw in fh:
            text = raw.decode('utf-8', 'replace')
            m = LINE.match(text)
            if m:
                message = text[m.end():]
                current = kind_of(message)
                lines[current] += 1
                junctions += current.startswith('junction start')
                suffix = message.rfind('~module:')
                overhead += m.end() + (len(message.rstrip('\n')) - suffix if suffix >= 0 else 0)
            size[current or 'other'] += len(raw)
    total = sum(size.values())
    rows = sorted(size, key=lambda k: -size[k])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'kinds.tsv').write_text('kind\tmessages\tbytes\tshare\tbytes_per_junction\n' + ''.join(
        f'{k}\t{lines[k]}\t{size[k]}\t{size[k] / total:.3f}\t{size[k] / max(junctions, 1):.0f}\n' for k in rows))
    summary = [f'{a.log}: {total / 1e9:.2f} GB, {sum(lines.values())} messages, {junctions} junctions, '
               f'{total / max(junctions, 1):.0f} bytes per junction',
               f'line prefix and suffix (time stamp, level, module, function): {overhead / 1e9:.2f} GB '
               f'({overhead / total:.1%})']
    summary += [f'  {k}: {lines[k]} messages, {size[k] / 1e6:.0f} MB ({size[k] / total:.1%})' for k in rows]
    (out / 'summary.txt').write_text('\n'.join(summary) + '\n')
    print('\n'.join(summary))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--log', required=True)
    p.add_argument('--out', default=str(Path(__file__).resolve().parent / 'genome'))
    main_(p.parse_args())
