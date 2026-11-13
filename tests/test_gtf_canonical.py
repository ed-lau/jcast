# -*- coding: utf-8 -*-

""" The GTF canonical from GENCODE's Ensembl_canonical tag when no transcript passes the rule (PR 17;
reports/2026-11-12-ensembl-canonical-audit.md) """

import logging
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import cache, params
from jcast.annots import ReadAnnotations, gene_rows
from jcast.sequences import Sequence, _ensembl_canonical_cds

# chr1: a short CDS, 101-112 and 201-209 (M A A K | W P E), stop codon 210-212; a longer one adds 301-309 (G G G),
# stop codon 310-312
CHR1 = ['C'] * 400
for start, bases in ((101, 'ATGGCTGCAAAA'), (201, 'TGGCCCGAATAA'), (301, 'GGTGGTGGTTAA')):
    CHR1[start - 1:start - 1 + len(bases)] = bases
GENOME = {'chr1': SeqRecord(Seq(''.join(CHR1)), id='chr1')}

SHORT = [(101, 112), (201, 209)]
LONG = [(101, 112), (201, 209), (301, 309)]


def _transcript(gene, tx, cds, stop, tsl='NA', tx_type='protein_coding', tags=(), start_codon=True):
    """ GENCODE-style rows of one transcript: transcript, CDS, start and stop codons """
    attrs = (f'gene_id "{gene}"; transcript_id "{tx}"; gene_type "protein_coding"; gene_name "{gene}"; '
             f'transcript_type "{tx_type}"; transcript_name "{tx}"; transcript_support_level "{tsl}";'
             + ''.join(f' tag "{t}";' for t in tags))
    rows = [('transcript', cds[0][0], stop[1] if stop else cds[-1][1], '.', attrs)]
    rows += [('CDS', s, e, '0', f'{attrs} exon_number "{i}"; protein_id "{tx}P";') for i, (s, e) in enumerate(cds, 1)]
    if start_codon:
        rows += [('start_codon', cds[0][0], cds[0][0] + 2, '0', f'{attrs} exon_number "1";')]
    if stop:
        rows += [('stop_codon', stop[0], stop[1], '0', f'{attrs} exon_number "{len(cds)}";')]
    return [f'chr1\tHAVANA\t{f}\t{s}\t{e}\t.\t+\t{frame}\t{a}\n' for f, s, e, frame, a in rows]


GTF = (
    # an ordinary gene with transcripts of TSL NA only: the tagged one, not the longest
    _transcript('ORDINARY', 'ORDINARY-201', SHORT, (210, 212), tags=('basic', 'Ensembl_canonical', 'MANE_Select')) +
    _transcript('ORDINARY', 'ORDINARY-202', LONG, (310, 312), tags=('basic',)) +
    # tagged, but a readthrough transcript
    _transcript('READTHROUGH', 'READTHROUGH-201', SHORT, (210, 212),
                tags=('basic', 'Ensembl_canonical', 'readthrough_transcript')) +
    # tags are matched whole
    _transcript('EXTENDED', 'EXTENDED-201', SHORT, (210, 212), tags=('basic', 'Ensembl_canonical_extended')) +
    # tagged, without a stop codon, or a start codon
    _transcript('NO_STOP', 'NO_STOP-201', SHORT, None, tags=('Ensembl_canonical', 'cds_end_NF')) +
    _transcript('NO_START', 'NO_START-201', SHORT, (210, 212), tags=('Ensembl_canonical',), start_codon=False) +
    # tagged, but not protein_coding
    _transcript('NMD', 'NMD-201', SHORT, (210, 212), tx_type='nonsense_mediated_decay', tags=('Ensembl_canonical',)) +
    # a transcript passes the rule: the rule's choice, not the tagged transcript
    _transcript('RULE', 'RULE-201', SHORT, (210, 212), tsl='1') +
    _transcript('RULE', 'RULE-202', LONG, (310, 312), tags=('basic', 'Ensembl_canonical'))
)


@pytest.fixture(scope='module')
def annotation(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('gtf_canonical')
    patch = pytest.MonkeyPatch()
    patch.setattr(cache, '_directory', str(tmp / 'cache'))
    (tmp / 'a.gtf').write_text(''.join(GTF))
    try:
        return ReadAnnotations(logging.getLogger('jcast-test'), tmp / 'a.gtf')
    finally:
        patch.undo()


def _sequence(gene):
    return Sequence(SimpleNamespace(gene_id=gene, gene_symbol=gene, chr='chr1', strand='+', name=gene))


@pytest.mark.parametrize('gene, chosen', [('ORDINARY', 'ORDINARY-201'), ('READTHROUGH', None), ('EXTENDED', None),
                                          ('NO_STOP', None), ('NO_START', None), ('NMD', None),
                                          ('RULE', 'RULE-201')])
def test_canonical_transcript(annotation, gene, chosen):
    sequence = _sequence(gene)
    sequence.get_annotated_transcripts_gtf(annotation)
    transcript = sequence.gtf_canonical_transcript
    assert (str(transcript) if transcript else None) == chosen


def test_canonical_sequence(annotation, monkeypatch):
    """ the tagged transcript is translated from its start codon, as the rule's choice is """
    monkeypatch.setattr(params, 'use_gtf_only', True)
    for gene, protein in (('ORDINARY', 'MAAKWPE'), ('RULE', 'MAAKWPE'), ('READTHROUGH', '')):
        sequence = _sequence(gene)
        sequence.get_canonical_aa(annotation, GENOME)
        assert (str(sequence.canonical_aa.seq), sequence.canonical_source) == \
               (protein, ('gtf', '') if protein else ('none', ''))


def test_gtf_without_tags(annotation):
    rows = gene_rows(annotation, 'ORDINARY')
    assert sorted(_ensembl_canonical_cds(rows)['transcript_name'].unique()) == ['ORDINARY-201']
    assert _ensembl_canonical_cds(rows.drop(columns='tag')).empty


def test_bundled_gtfs_have_no_ensembl_canonical(gtf, rat_gtf):
    """ the golden tests' GTFs (Ensembl 89 and 100) do not carry the tag, so the fallback never applies to them """
    for annotation in (gtf, rat_gtf):
        tags = annotation.annot['tag'].dropna().str.split(',') if 'tag' in annotation.annot.columns else []
        assert not any('Ensembl_canonical' in t for t in tags)
