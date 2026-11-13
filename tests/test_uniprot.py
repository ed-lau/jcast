# -*- coding: utf-8 -*-

""" Canonicals from a supplied UniProt FASTA (--uniprot) and the merged database (PR 14;
reports/2026-11-11-uniprot-fasta-evaluation.md) """

import argparse
import collections
import csv
import hashlib
import logging

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from jcast import main, params
from jcast.sequences import Sequence
from jcast.uniprot import UniProtFasta

from tests.conftest import DATA, HUMAN_GENOME, HUMAN_GTF, RMATS
from tests.test_junction import _junction

# Swiss-Prot human with isoforms, release 2026_04, restricted to the chr15 test genes (tests/data/uniprot/README.md)
UNIPROT = DATA / 'uniprot' / 'uniprot_sprot_human_isoforms_chr15_2026_04.fasta.gz'


def _no_api(self, reviewed='true'):
    raise AssertionError('UniProt API called with --uniprot')


def _write_fasta(path, entries):
    with open(path, 'w') as fh:
        for header, seq in entries:
            fh.write(f'>{header}\n{seq}\n')
    return path


def _read(path):
    return [(r.description, str(r.seq)) for r in SeqIO.parse(path, 'fasta')]


@pytest.fixture(scope='module')
def uniprot():
    return UniProtFasta(logger=logging.getLogger('jcast-test'), path=UNIPROT)


def test_renamed_gene_is_mapped_by_sequence(gtf, genome, rmats_results, uniprot, monkeypatch):
    """ FAM214A in Ensembl 89 is ATOSA in UniProt: no gene name match, but its canonical is a GTF protein """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)

    junction = _junction(rmats_results.rmats_se, 41203)
    sequence = Sequence(junction, uniprot=uniprot)
    sequence.get_canonical_aa(gtf, genome.genome)

    assert junction.gene_symbol == 'FAM214A'
    assert sequence.canonical_aa.id == 'sp|Q32MH5|ATOSA_HUMAN'
    assert uniprot.chosen[junction.gene_id][1] == 'sp by sequence'


def test_gene_without_entry_uses_gtf(gtf, genome, rmats_results, uniprot, monkeypatch):
    """ JMJD7-PLA2G4B (readthrough) has no Swiss-Prot entry: the GTF canonical is used """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)

    rma = rmats_results.rmats_se
    junction = _junction(rma, rma[rma['gene_symbol'] == 'JMJD7-PLA2G4B']['id'].iloc[0])
    sequence = Sequence(junction, uniprot=uniprot)
    sequence.get_canonical_aa(gtf, genome.genome)

    assert uniprot.chosen[junction.gene_id][1:3] == ('gtf', '')
    assert str(sequence.canonical_aa.seq) == str(sequence.get_canonical_aa_gtf(gtf, genome.genome).seq)
    assert len(sequence.canonical_aa) > 0


def test_choice_order(tmp_path):
    """ by sequence before by gene name, reviewed before unreviewed, a unique match only, canonical entries only """
    fasta = _write_fasta(tmp_path / 'uniprot.fasta', [
        ('sp|P00001|A_HUMAN Protein A OS=Homo sapiens GN=GENEA', 'MAAAA'),
        ('sp|P00001-2|A_HUMAN Isoform 2 of Protein A OS=Homo sapiens GN=GENEA', 'MAAAB'),
        ('sp|P00002|B_HUMAN Protein B OS=Homo sapiens GN=GENEB', 'MBBBB'),
        ('sp|P00003|C1_HUMAN Protein C1 OS=Homo sapiens GN=GENEC', 'MCCC1'),
        ('sp|P00004|C2_HUMAN Protein C2 OS=Homo sapiens GN=GENEC', 'MCCC2'),
        ('tr|A00005|A00005_HUMAN Protein C OS=Homo sapiens GN=GENEC', 'MCCC3'),
    ])
    uniprot = UniProtFasta(logger=logging.getLogger('jcast-test'), path=fasta)
    gtf = SeqRecord(Seq('MGTF'), id='gtf')

    def choose(gene, proteins):
        record = uniprot.canonical(gene, gene, lambda: proteins, lambda: gtf)
        return record.id, uniprot.chosen[gene][1]

    # sequence wins over a different gene name
    assert choose('GENEB', {'MAAAA'}) == ('sp|P00001|A_HUMAN', 'sp by sequence')
    # an isoform's sequence does not make it the canonical; the gene name does
    assert choose('GENEA', {'MAAAB'}) == ('sp|P00001|A_HUMAN', 'sp by gene name')
    # two reviewed entries share the gene name: not used; the unreviewed entry has a unique name
    assert choose('GENEC', set()) == ('tr|A00005|A00005_HUMAN', 'tr by gene name')
    # nothing matches: the GTF canonical
    assert choose('GENED', {'MDDDD'}) == ('gtf', 'gtf')


def test_merged_database(tmp_path):
    """ Swiss-Prot entries, then psq_T1 and psq_T1_alt_* records joined to a Swiss-Prot canonical that are new """
    fasta = _write_fasta(tmp_path / 'uniprot.fasta', [
        ('sp|P00001|A_HUMAN Protein A OS=Homo sapiens GN=GENEA', 'MAAAAAAAA'),
        ('sp|P00001-2|A_HUMAN Isoform 2 of Protein A OS=Homo sapiens GN=GENEA', 'MAAAAKAAA'),
    ])
    out = tmp_path / 'out'
    out.mkdir()
    _write_fasta(out / 'psq_T1.fasta', [
        ('sp|P00001|A_HUMAN|G1|SE1|1|r9|T1', 'MAAAAKAAA'),             # identical to isoform 2: not added
        ('sp|P00001|A_HUMAN|G1|SE2|1|r9|T1', 'MAAAANEWA'),             # P00001-J1
        ('tr|A00005|A00005_HUMAN|G5|SE1|5|r9|T1', 'MTTTTNEWT'),        # TrEMBL: not added
    ])
    _write_fasta(out / 'psq_T1_alt_cterm.fasta', [
        ('sp|P00001|A_HUMAN|G1|A3SS1|2|r9|T1_alt_cterm', 'MAAAANEWA'),  # same as a psq_T1 record: not added
        ('sp|P00001|A_HUMAN|G1|A3SS2|3|r9|T1_alt_cterm', 'MAAAANEWB'),  # P00001-J2
    ])
    _write_fasta(out / 'psq_T1_orphan.fasta', [('xx|ORPHN|GENEA|G1|SE1|4|r9|T1', 'NEWNEW')])
    _write_fasta(out / 'psq_T4.fasta', [('sp|P00001|A_HUMAN|G1|SE1|6|r9|T4', 'MAAAANEWC')])     # not a merged tier
    _write_fasta(out / 'psq_T2.fasta', [('sp|P00001|A_HUMAN|G1|SE1|7|r9|T2', 'MAAAAAAAA')])     # identical: listed

    uniprot = UniProtFasta(logger=logging.getLogger('jcast-test'), path=fasta)
    uniprot.write_outputs(out)

    assert _read(out / 'jcast_uniprot_merged.fasta') == [
        ('sp|P00001|A_HUMAN Protein A OS=Homo sapiens GN=GENEA', 'MAAAAAAAA'),
        ('sp|P00001-2|A_HUMAN Isoform 2 of Protein A OS=Homo sapiens GN=GENEA', 'MAAAAKAAA'),
        ('sp|P00001-J1|A_HUMAN sp|P00001|A_HUMAN|G1|SE2|1|r9|T1', 'MAAAANEWA'),
        ('sp|P00001-J2|A_HUMAN sp|P00001|A_HUMAN|G1|A3SS2|3|r9|T1_alt_cterm', 'MAAAANEWB'),
    ]
    rows = list(csv.DictReader(open(out / 'psq_isoforms.tsv'), delimiter='\t'))
    assert [(r['file'], r['accessions'], r['entry']) for r in rows] == [
        ('psq_T1.fasta', 'P00001-2', 'isoform'),
        ('psq_T2.fasta', 'P00001', 'canonical'),
    ]


def test_cli_rejects_missing_uniprot(monkeypatch, capsys):
    monkeypatch.setattr('sys.argv', ['jcast', str(RMATS), str(HUMAN_GTF), str(HUMAN_GENOME),
                                     '--uniprot', 'missing.fasta'])

    with pytest.raises(SystemExit):
        main.main()
    assert 'file not found: missing.fasta' in capsys.readouterr().err


EXPECTED_SOURCES = {'sp by sequence': 267, 'sp by gene name': 8, 'gtf': 1, 'none': 18}

# file name -> (number of records, sha256 of the sorted sequences). PR 14 pinned these at -r 1 keeping counts >= 2, when
# the psq_* counts equalled those with the full Swiss-Prot release (reports/uniprot_eval/final/summary.txt); since PR 6,
# -r 1 keeps every junction (reports/2026-11-12-thresholds-runs-canonical-audit.md §7)
EXPECTED_OUTPUT = {
    'jcast_uniprot_merged.fasta': (893, '9c31b68b076a7d7c'),
    'psq_T1.fasta': (253, '7daf17dd26a5b455'),
    'psq_T1_alt_cterm.fasta': (18, 'd94f4fd7852b9ba4'),
    'psq_T1_alt_cterm_orphan.fasta': (2, 'af4d8be95cb5dd05'),
    'psq_T1_alt_nterm.fasta': (6, '4ecca3ae6061a2d6'),
    'psq_T1_orphan.fasta': (28, 'bf2cbed6e747a093'),
    'psq_T2.fasta': (47, '068c0166f5e601b9'),
    'psq_T2_orphan.fasta': (6, '105037b254205e97'),
    'psq_T4.fasta': (351, '5057eeedae64959e'),
    'psq_T4_orphan.fasta': (77, 'ace037b31473c649'),
    'psq_canonical.fasta': (276, '9f547139112f68ae'),
}


def _digest(path):
    seqs = sorted(str(r.seq) for r in SeqIO.parse(path, 'fasta'))
    return len(seqs), hashlib.sha256('\n'.join(seqs).encode()).hexdigest()[:16]


@pytest.mark.slow
def test_chr15_uniprot(tmp_path, monkeypatch, cache_dir):
    """ the whole chr15 test set through run_jcast with --uniprot: no API call, pinned canonicals and outputs """
    monkeypatch.setattr(Sequence, 'get_canonical_aa_uniprot', _no_api)
    monkeypatch.setattr(params, 'use_gtf_only', False)
    out = tmp_path / 'out'
    args = argparse.Namespace(rmats_folder=str(RMATS), gtf_file=str(HUMAN_GTF), genome=str(HUMAN_GENOME),
                              out=str(out), read=1, model=False, canonical=True, mask=None, qvalue=[0, 1],
                              splice_type={'MXE', 'SE', 'RI', 'A5SS', 'A3SS'}, uniprot=str(UNIPROT),
                              cache_dir=str(cache_dir))

    main.run_jcast(args)

    sources = collections.Counter(r['source'] for r in
                                  csv.DictReader(open(out / 'uniprot_canonicals.tsv'), delimiter='\t'))
    assert dict(sources) == EXPECTED_SOURCES
    assert {p.name: _digest(p) for p in sorted(out.glob('*.fasta'))} == EXPECTED_OUTPUT
