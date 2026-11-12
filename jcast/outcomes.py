# -*- coding: utf-8 -*-

""" The run record: junctions.tsv, one row per junction, and the end-of-run summary (reports/2026-11-12-log-audit.md) """

import collections
import csv
import glob
import os

from Bio import SeqIO

from jcast import fates

JUNCTIONS_FILE = 'junctions.tsv'

SPLICE_TYPES = ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']

COLUMNS = ['junction_type',         # MXE, SE, RI, A5SS or A3SS
           'junction_id',           # the rMATS ID
           'gene_id',
           'gene_symbol',
           'chr',
           'strand',
           'min_jc',                # the read count: the smaller of the summed inclusion and skipping junction counts
           'fdr',                   # the rMATS FDR
           'fate',                  # the fate's short name (fates.NAMES)
           'fate_detail',           # T4: the slice that hit a premature stop; exception: the exception
           'slice1_phase',          # each slice's annotated frame (0, 1 or 2), empty if none or not looked up
           'slice2_phase',
           'cds_start',             # the coding region the slices were trimmed to, empty if not looked up
           'cds_end',
           'canonical_source',      # where the gene's canonical came from, empty if not looked up
           'canonical_accession',
           'files',                 # the psq_<name>.fasta files of the junction's records, as <name>, ;-separated
           'new_records',           # how many of its records were written: a record already in its file is not
           ]


def outcome(junction) -> dict:
    """ What translating the junction found, for its row: read from the junction after _translate_one """
    return dict(slice_phases=junction.slice_phases, cds=(junction.tx0, junction.tx1), canonical=junction.canonical)


class JunctionTable(object):
    """ junctions.tsv, written row by row in junction order, and the counts for the summary """

    def __init__(self, out_dir):
        self.path = os.path.join(out_dir, JUNCTIONS_FILE)
        self._fh = open(self.path, 'w', newline='')
        self._writer = csv.writer(self._fh, delimiter='\t', lineterminator='\n')
        self._writer.writerow(COLUMNS)
        self.counts = collections.Counter()     # (fate name, splice type) -> junctions

    def add(self, junction, fate, found, suffixes, new_records) -> None:
        """
        Add the junction's row

        :param junction: Junction
        :param fate: its fate
        :param found: outcome(junction), or {} if it failed with an exception
        :param suffixes: the file suffix of each of its records
        :param new_records: how many of its records were written
        """
        short, detail = fates.name(fate)
        self.counts[short, junction.junction_type] += 1
        phases = found.get('slice_phases', (None, None))
        cds = found.get('cds', (None, None))
        canonical = found.get('canonical') or (None, None)
        self._writer.writerow(['' if v is None else v for v in [
            junction.junction_type, junction.name, junction.gene_id, junction.gene_symbol, junction.chr,
            junction.strand, junction.min_jc, junction.fdr, short, detail, *phases, *cds, *canonical,
            ';'.join(dict.fromkeys(suffixes)), new_records]])

    def close(self) -> None:
        self._fh.close()

    def summary(self, splice_types, out_dir) -> list:
        """
        The end-of-run summary: junctions by fate and splice type, and the records in each output FASTA

        :param splice_types: the splice types run
        :param out_dir: the output directory
        :return: lines
        """
        types = [t for t in SPLICE_TYPES if t in splice_types]
        rows = [['fate'] + types + ['total']]
        for short in [s for _, s in fates.NAMES] + ['total']:
            n = [sum(c for (f, t), c in self.counts.items() if t == type_ and short in (f, 'total')) for type_ in types]
            rows.append([short] + [f'{c:,}' for c in n + [sum(n)]])
        width = max(len(r[0]) for r in rows)
        widths = [max(len(r[i]) for r in rows) for i in range(1, len(rows[0]))]
        lines = [r[0].ljust(width) + ''.join(v.rjust(w + 2) for v, w in zip(r[1:], widths)) for r in rows]

        records = []
        for path in sorted(glob.glob(os.path.join(glob.escape(out_dir), 'psq_*.fasta'))):
            with open(path) as fh:
                records.append(f'{os.path.basename(path)} {sum(1 for _ in SeqIO.parse(fh, "fasta")):,}')

        return (['Junctions by fate and splice type (each junction is a row of {0}):'.format(JUNCTIONS_FILE)] + lines
                + ['Records per file: ' + (', '.join(records) if records else 'none')])
