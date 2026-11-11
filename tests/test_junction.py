# -*- coding: utf-8 -*-

""" Tests for splice junction related functions """

import argparse
import shutil
import types

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from jcast import fates, main, params
from jcast.junctions import Junction, RmatsResults
from jcast.sequences import Sequence

from tests.conftest import RMATS


@pytest.mark.network
def test_pkm_manual(gtf, genome, rmats_results):
    """ test for translating PKM, a known AS protein """

    rma = rmats_results.rmats_mxe
    # TODO: test negative strand and trimming
    i = 71  # PKM

    junction = Junction(**rma.iloc[i].to_dict())

    junction.trim_cds(gtf)
    junction.get_translated_phase(gtf)

    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)
    sequence.get_canonical_aa(gtf, genome.genome)
    sequence.translate()

    for i in [1, 2]:
        sequence.stitch_to_canonical_aa(slice_to_stitch=i)

    # Compare with Uniprot retrieved sequences
    retrieved_pkm2 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                      'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                      'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                      'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                      'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                      'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                      'LSGETAKGDYPLEAVRMQHLIAREAEAAIYHLQLFEELRRLAPITSDPTEATAVGAVEAS'
                      'FKCCSGAIIVLTKSGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                      'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

    retrieved_pkm1 = ('MSKPHSEAGTAFIQTQQLHAAMADTFLEHMCRLDIDSPPITARNTGIICTIGPASRSVET'
                      'LKEMIKSGMNVARLNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIR'
                      'TGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVEVGSKIYVDDGL'
                      'ISLQVKQKGADFLVTEVENGGSLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMV'
                      'FASFIRKASDVHEVRKVLGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIE'
                      'IPAEKVFLAQKMMIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIM'
                      'LSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELVRASSHSTDLMEAMAMGSVEAS'
                      'YKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNPQTARQAHLYRGIFPVLCKDPVQE'
                      'AWAEDVDLRVNFAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP')

    assert sequence.slice1_stitched.seq == retrieved_pkm2
    assert sequence.slice2_stitched.seq == retrieved_pkm1


def test_trim_start(gtf, rmats_results):
    """ test for trimming junctions """

    rma = rmats_results.rmats_a5ss
    # TODO: test negative strand and trimming
    i = 6  # SEMA4B A5SS 742

    junction = Junction(**rma.iloc[i].to_dict())
    junction.trim_cds(gtf)

    assert junction.tx0 == 90201579


def _junction(rma, rmats_id):
    return Junction(**rma[rma['id'] == rmats_id].iloc[0].to_dict())


def test_c_terminal_bound_plus_strand(gtf, rmats_results):
    """ tx1 is the end of the stop codon on the + strand (REVIEW.md §5.1) """

    junction = _junction(rmats_results.rmats_se, 6426)  # CASC4
    junction.trim_cds(gtf)

    assert (junction.tx0, junction.tx1) == (44289030, 44413406)
    assert junction.down_ee == 44413406


def test_c_terminal_bound_minus_strand(gtf, rmats_results):
    """ tx0 is the start of the stop codon on the - strand (REVIEW.md §5.1) """

    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8
    junction.trim_cds(gtf)

    assert (junction.tx0, junction.tx1) == (88899395, 88913319)
    assert junction.anc_es == 88899395


def test_c_terminal_event_translates(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ an event spanning the stop-codon exon is no longer read as a premature stop """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=False,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8

    fate = main._translate_one(junction, gtf, genome, args, pred_bound=args.read)

    assert fate.startswith('SUCCESS TIER 1')


def test_bounds_without_tsl_column(gtf, rmats_results):
    """ a GTF without transcript_support_level still yields bounds (REVIEW.md MED-7) """

    no_tsl = types.SimpleNamespace(annot=gtf.annot.drop(columns='transcript_support_level'))
    junction = _junction(rmats_results.rmats_se, 6426)  # CASC4
    junction.trim_cds(no_tsl)

    assert (junction.tx0, junction.tx1) == (44289030, 44413406)


def _exon(junction, exon):
    return getattr(junction, f'{exon}_es'), getattr(junction, f'{exon}_ee')


def test_trim_drops_utr_exon(gtf, rmats_results):
    """ exons lying entirely outside the coding region are removed (REVIEW.md PR 4a) """

    junction = _junction(rmats_results.rmats_se, 6443)  # TMEM62, +: anchor is 5' UTR
    junction.trim_cds(gtf)
    assert _exon(junction, 'anc') == (-1, -1)
    assert _exon(junction, 'alt1') == (43134257, 43134368)

    junction = _junction(rmats_results.rmats_se, 6654)  # UBE3A, -: downstream exon is 5' UTR
    junction.trim_cds(gtf)
    assert _exon(junction, 'down') == (-1, -1)
    assert _exon(junction, 'anc') == (25375465, 25375763)


def test_trim_minus_strand_alternative_sites(gtf, rmats_results):
    """ A5SS/A3SS on the - strand are trimmed at the bound their exons actually contain """

    junction = _junction(rmats_results.rmats_a5ss, 3645)  # WDR73: start codon in the alternative exons
    junction.trim_cds(gtf)
    assert _exon(junction, 'alt1') == (84654227, 84654274)
    assert _exon(junction, 'alt2') == (84654234, 84654274)

    junction = _junction(rmats_results.rmats_a3ss, 1238)  # HDDC3: stop codon in the alternative exons
    junction.trim_cds(gtf)
    assert _exon(junction, 'alt1') == (90931275, 90931556)
    assert _exon(junction, 'alt2') == (90931275, 90931405)


def test_trim_retained_intron_region(gtf, rmats_results):
    """ the retained region of an RI event is trimmed along with its flanking exons """

    junction = _junction(rmats_results.rmats_ri, 3355)  # B2M, +: stop codon in the downstream exon
    junction.trim_cds(gtf)
    assert _exon(junction, 'down') == (44716329, 44716342)
    assert _exon(junction, 'alt1') == (44715423, 44716342)


def test_event_outside_cds_is_skipped(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a splice event lying entirely in the UTR is skipped; the canonical is still written with -c """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=True,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_mxe, 433)  # ZNF592: both alternative exons are 5' UTR

    fate = main._translate_one(junction, gtf, genome, args, pred_bound=args.read)

    assert fate == fates.SKIPPED_OUTSIDE_CDS
    assert [p.name for p in tmp_path.glob('psq_*.fasta')] == ['psq_canonical.fasta']


def test_retained_intron_slices(genome, rmats_results):
    """ RI slice 1 is the retained exon alone; slice 2 splices the flanking exons (REVIEW.md HIGH-11) """

    junction = _junction(rmats_results.rmats_ri, 916)  # UNC45A, +
    sequence = Sequence(junction)
    sequence.make_slice_localgenome(genome.genome)

    chrom = genome.genome['15'].seq
    retained = chrom[junction.alt1_es - 1:junction.alt1_ee]
    spliced = chrom[junction.anc_es - 1:junction.anc_ee] + chrom[junction.down_es - 1:junction.down_ee]

    assert (len(sequence.slice1_nt), len(sequence.slice2_nt)) == (515, 37 + 324)
    assert sequence.slice1_nt.seq == retained
    assert sequence.slice2_nt.seq == spliced


def test_retained_intron_in_utr_is_skipped(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ an intron retained upstream of the start codon is outside the CDS; it used to yield a doubled slice """

    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(tmp_path), read=1, model=False, canonical=False,
                              mask='NnXx', qvalue=[0, 1])
    junction = _junction(rmats_results.rmats_ri, 943)  # FAH, +

    assert main._translate_one(junction, gtf, genome, args, pred_bound=args.read) == fates.SKIPPED_OUTSIDE_CDS


def _annotated_protein(gtf, genome, transcript_name):
    """ the protein encoded by an annotated transcript, translated from its CDS rows """

    cds = gtf.annot.query('transcript_name == @transcript_name & feature == "CDS"').sort_values('start')
    chrom = genome.genome[str(cds['seqname'].iloc[0])].seq
    nt = ''.join(str(chrom[s - 1:e]) for s, e in zip(cds['start'], cds['end']))
    if cds['strand'].iloc[0] == '-':
        nt, frame = str(Seq(nt).reverse_complement()), int(cds['frame'].iloc[-1])
    else:
        frame = int(cds['frame'].iloc[0])
    nt = nt[frame:]
    return str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip('*')


def _translate(junction, gtf, genome, out_dir, monkeypatch):
    monkeypatch.setattr(params, 'use_gtf_only', True)
    args = argparse.Namespace(out=str(out_dir), read=1, model=False, canonical=False,
                              mask='NnXx', qvalue=[0, 1])
    return main._translate_one(junction, gtf, genome, args, pred_bound=args.read)


def test_phase_unchanged_where_it_was_found(gtf, rmats_results):
    """ the per-slice lookup returns the frame the exact-match lookup found (PR 4b audit §5.1) """

    junction = _junction(rmats_results.rmats_se, 5739)  # MFGE8, -
    junction.trim_cds(gtf)
    junction.get_translated_phase(gtf)

    assert junction.slice_phases == (2, 2)
    assert junction.phase == 2


def test_per_slice_phase(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ slices that start in different exons are each read in their own frame (REVIEW.md PR 4b) """

    junction = _junction(rmats_results.rmats_se, 11945)  # GTF2A2, -: 5' UTR exon removed, slices start apart
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert junction.slice_phases == (0, 2)
    assert fate.startswith('SUCCESS TIER 1')  # in frame once each slice's frame is accounted for (PR 4b2)


def test_per_slice_phase_reconstructs_annotated_isoform(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ LINS1 MXE 449 slice 2, read in its own frame, is exactly the annotated LINS1-207 protein (T1 since PR 4b2) """

    junction = _junction(rmats_results.rmats_mxe, 449)
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert junction.slice_phases == (2, 0)
    assert fate.startswith('SUCCESS TIER 1')
    written = {str(r.seq) for r in SeqIO.parse(tmp_path / 'psq_T1.fasta', 'fasta')}
    assert _annotated_protein(gtf, genome, 'LINS1-207') in written


def test_slice_without_frame_is_not_guessed(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ only the slice with an annotated frame is translated; the other is not guessed (REVIEW.md PR 4b) """

    junction = _junction(rmats_results.rmats_se, 12446)  # ANXA2, -: no CDS covers slice 1's start
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert junction.slice_phases == (None, 2)
    assert fate.startswith('PARTIAL TIER 4')


def test_no_coding_annotation_is_not_translated(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a gene without any protein-coding CDS gets no frame and fails, rather than a guessed translation """

    junction = _junction(rmats_results.rmats_se, 5796)  # GOLGA2P10, a pseudogene
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert junction.slice_phases == (None, None)
    assert fate == fates.FAIL


def test_annotated_frame_is_not_overridden(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ a slice that hits a stop in its annotated frame is T4, not re-translated in another frame (former T3) """

    junction = _junction(rmats_results.rmats_se, 6181)  # USP3
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert fate.startswith('PARTIAL TIER 4')
    assert not list(tmp_path.glob('psq_T3*.fasta'))


def test_frameshift_counts_each_slice_frame(gtf, genome, rmats_results, tmp_path, monkeypatch):
    """ slices whose raw lengths differ by a multiple of 3 are still a frameshift if their frames differ (PR 4b2) """

    junction = _junction(rmats_results.rmats_se, 11950)  # GTF2A2, -: frames (0, 2)
    fate = _translate(junction, gtf, genome, tmp_path, monkeypatch)

    assert junction.slice_phases == (0, 2)
    assert fate.startswith('SUCCESS TIER 2')


@pytest.mark.parametrize('na', ['NA', '"NA"', ''])
def test_missing_gene_symbol(na, gtf, genome, logger, tmp_path, monkeypatch):
    """ rMATS writes NA when a gene has no symbol; pandas reads it as NaN (REVIEW.md, PR 2) """

    rmats_dir = tmp_path / 'rmats'
    shutil.copytree(RMATS, rmats_dir)
    se = rmats_dir / 'SE.MATS.JC.txt'
    lines = se.read_text().splitlines()
    for i, line in enumerate(lines):
        fields = line.split('\t')
        if fields[0] == '5739':  # MFGE8, ENSG00000140545
            lines[i] = '\t'.join(fields[:2] + [na] + fields[3:])
    se.write_text('\n'.join(lines) + '\n')

    junction = _junction(RmatsResults(logger=logger, rmats_dir=rmats_dir).rmats_se, 5739)
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    _translate(junction, gtf, genome, out_dir, monkeypatch)

    assert junction.gene_symbol == 'NoGeneName_ENSG00000140545'
    assert 'NoGeneName_ENSG00000140545' in str(junction)
    headers = [r.id for r in SeqIO.parse(out_dir / 'psq_T1.fasta', 'fasta')]
    assert headers and all('|NoGeneName_ENSG00000140545|' in h for h in headers)
