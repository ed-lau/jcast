# -*- coding: utf-8 -*-

""" Per-gene GTF index and in-memory de-duplication of FASTA writes (PR 9; reports/2026-11-12-parallel-audit.md) """

import types

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import helpers as h
from jcast.annots import gene_rows


def test_gene_rows_equal_the_filter(gtf):
    """ the index gives exactly the rows, in file order, that filtering the whole GTF gives """
    for gene_id in list(gtf.by_gene)[::20]:
        assert gene_rows(gtf, gene_id).equals(gtf.annot.query('gene_id == @gene_id'))
    assert len(gtf.by_gene) == gtf.annot['gene_id'].nunique()


def test_gene_rows_missing_gene(gtf):
    rows = gene_rows(gtf, 'ENSG_NOT_IN_GTF')
    assert rows.empty and list(rows.columns) == list(gtf.annot.columns)


def test_gene_rows_without_index(gtf):
    """ an annotation without the index (e.g. a test double) is filtered as before """
    plain = types.SimpleNamespace(annot=gtf.annot)
    gene_id = list(gtf.by_gene)[100]
    assert gene_rows(plain, gene_id).equals(gtf.by_gene[gene_id])


def _record(seq, name):
    return SeqRecord(Seq(seq), id=name, description='')


def _sequences(path):
    return [str(r.seq) for r in SeqIO.parse(path, 'fasta')]


def test_duplicate_sequence_not_written(tmp_path):
    for seq, name in [('MAAA', 'a'), ('MBBB', 'b'), ('MAAA', 'c')]:
        h.write_seqrecord_to_fasta(_record(seq, name), str(tmp_path), 'T1')
    assert _sequences(tmp_path / 'psq_T1.fasta') == ['MAAA', 'MBBB']


def test_file_from_an_earlier_run_is_read(tmp_path):
    """ an existing output file is read before the first write, as before """
    (tmp_path / 'psq_T1.fasta').write_text('>old\nMOLD\n')
    h.write_seqrecord_to_fasta(_record('MOLD', 'new'), str(tmp_path), 'T1')
    h.write_seqrecord_to_fasta(_record('MNEW', 'new'), str(tmp_path), 'T1')
    assert _sequences(tmp_path / 'psq_T1.fasta') == ['MOLD', 'MNEW']


def test_file_changed_or_removed_between_writes(tmp_path):
    """ the cache is not trusted when the file is no longer the one we wrote """
    h.write_seqrecord_to_fasta(_record('MAAA', 'a'), str(tmp_path), 'T1')
    (tmp_path / 'psq_T1.fasta').unlink()
    h.write_seqrecord_to_fasta(_record('MAAA', 'a'), str(tmp_path), 'T1')
    assert _sequences(tmp_path / 'psq_T1.fasta') == ['MAAA']

    with open(tmp_path / 'psq_T1.fasta', 'a') as fh:
        fh.write('>other\nMOTHER\n')
    h.write_seqrecord_to_fasta(_record('MOTHER', 'b'), str(tmp_path), 'T1')
    assert _sequences(tmp_path / 'psq_T1.fasta') == ['MAAA', 'MOTHER']
