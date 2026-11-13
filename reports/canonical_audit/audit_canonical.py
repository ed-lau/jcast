# -*- coding: utf-8 -*-

"""Audit for PR 17: the GTF canonical from GENCODE's Ensembl_canonical tag
(reports/2026-11-12-ensembl-canonical-audit.md).

Two parts:

    genes     From the GTF, and a whole-genome run with -c --uniprot (its uniprot_canonicals.tsv and psq_*.fasta):
                - for every protein-coding gene in the run without a canonical, why the GTF rule gives none, what its
                  Ensembl_canonical transcript is, and what each candidate rule would give it;
                - for the genes the rule does give a GTF canonical, whether that is the Ensembl_canonical
                  transcript's coding sequence.
              With --rmats, also: whether rMATS reports an event under more than one gene, and how many SE events
              of the readthrough genes without a canonical skip an exon annotated in another protein-coding gene.
              Writes genes.tsv and summary.txt.
    compare   Two output directories of the same run, before and with a prototype rule (--before, --after): the files
              that differ, canonical sources changed, and records by gene, joined or orphan. Writes compare_<label>.txt.

The current rule (Sequence.get_annotated_transcripts_gtf): the longest coding sequence among the gene's transcripts
of type protein_coding with CDS rows, transcript support level <= params.tsl_threshold, and annotated start and stop
codons. A tag is matched as a whole: GENCODE also has Ensembl_canonical_extended.

Usage (from the repository root):
    python reports/canonical_audit/audit_canonical.py genes --gtf <gtf> --run <output dir> [--rmats <dir>] [--out <dir>]
    python reports/canonical_audit/audit_canonical.py compare --before <dir> --after <dir> --label <name> [--out <dir>]
"""

import argparse
import collections
import csv
import logging
import sys
from pathlib import Path

from Bio import SeqIO

REPO = Path(__file__).resolve().parents[2]
TAG = 'Ensembl_canonical'


def _tags(value) -> set:
    return set(value.split(',')) if isinstance(value, str) and value else set()


def transcripts(annot, tsl_threshold):
    """ one row per transcript, with what the rules look at """
    cols = ['gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_name', 'transcript_type',
            'transcript_support_level', 'tag']
    tx = annot.loc[annot['feature'] == 'transcript', cols].copy()
    tx['tags'] = tx['tag'].map(_tags)
    cds = annot.loc[annot['feature'] == 'CDS', ['transcript_id', 'start', 'end']]
    tx['cds_length'] = tx['transcript_id'].map((cds['end'] - cds['start'] + 1).groupby(cds['transcript_id']).sum()) \
        .fillna(0).astype(int)
    structure = cds.sort_values('start').groupby('transcript_id').apply(
        lambda t: tuple(zip(t['start'], t['end'])), include_groups=False)
    tx['cds_structure'] = tx['transcript_id'].map(structure)
    for feature in ['start_codon', 'stop_codon']:
        tx[feature] = tx['transcript_id'].isin(set(annot.loc[annot['feature'] == feature, 'transcript_id']))
    tx['coding'] = (tx['transcript_type'] == 'protein_coding') & (tx['cds_length'] > 0)
    tx['supported'] = tx['transcript_support_level'] <= tsl_threshold
    tx['rule'] = tx['coding'] & tx['supported'] & tx['start_codon'] & tx['stop_codon']
    tx['ensembl_canonical'] = tx['tags'].map(lambda t: TAG in t)
    return tx


def why_none(g) -> str:
    """ why the current rule gives the gene no transcript """
    if not (g['cds_length'] > 0).any():
        return 'no transcript with a CDS'
    if not g['coding'].any():
        return 'CDS only in transcripts not of type protein_coding'
    if not (g['coding'] & g['supported']).any():
        return 'protein_coding transcripts of TSL 3-5 or NA only'
    if not g['rule'].any():
        return 'TSL <= 2 protein_coding transcripts lack a start or stop codon'
    return 'the rule finds a transcript, but its translation is empty'


def candidate(ec) -> dict:
    """ the gene's Ensembl_canonical transcript, and what each candidate rule would take """
    if ec is None:
        return dict(ec_transcript='', ec_type='', ec_tsl='', ec_cds='', ec_start_stop='', ec_tags='',
                    v1_complete='no', v2_any_cds='no')
    complete = bool(ec['start_codon'] and ec['stop_codon'])
    shown = sorted(t for t in ec['tags'] if t in {'MANE_Select', 'basic', 'CCDS', 'cds_start_NF', 'cds_end_NF',
                                                   'mRNA_start_NF', 'mRNA_end_NF', 'readthrough_transcript'})
    tsl = ec['transcript_support_level']
    return dict(ec_transcript=ec['transcript_name'], ec_type=ec['transcript_type'],
                ec_tsl='NA' if tsl != tsl else str(int(tsl)), ec_cds=ec['cds_length'],
                ec_start_stop=('start ' if ec['start_codon'] else '') + ('stop' if ec['stop_codon'] else ''),
                ec_tags=','.join(shown),
                v1_complete='yes' if ec['coding'] and complete else 'no',
                v2_any_cds='yes' if ec['cds_length'] > 0 else 'no')


EVENT = __import__('re').compile(r'^(MXE|SE|RI|A5SS|A3SS)[12]$')


def events(out_dir):
    """
    (splice type and slice, chromosome, anchor, alternative exon) -> [(gene id, joined or orphan)] for every record in
    the run's psq files but psq_canonical, from its header
    """
    found = collections.defaultdict(list)
    for f in sorted(Path(out_dir).glob('psq_*.fasta')):
        if f.name == 'psq_canonical.fasta':
            continue
        for r in SeqIO.parse(f, 'fasta'):
            fields = r.id.split('|')
            i = next(i for i, x in enumerate(fields) if EVENT.match(x))
            found[fields[i], fields[i + 2], fields[i + 3], fields[i + 4]].append(
                (fields[i - 1], 'orphan' if 'orphan' in f.name else 'joined'))
    return found


def output_records(out_dir, gene_ids):
    """ gene id -> Counter of records in the run's psq files but psq_canonical, as 'joined' or 'orphan' """
    records = collections.defaultdict(collections.Counter)
    for f in sorted(Path(out_dir).glob('psq_*.fasta')):
        if f.name == 'psq_canonical.fasta':
            continue
        for r in SeqIO.parse(f, 'fasta'):
            gene = next((x for x in r.id.split('|') if x in gene_ids), None)
            if gene:
                records[gene]['orphan' if 'orphan' in f.name else 'joined'] += 1
    return records


def genes(a):
    sys.path.insert(0, str(REPO))
    from jcast import params
    from jcast.annots import ReadAnnotations

    annot = ReadAnnotations(logging.getLogger('jcast-audit'), a.gtf).annot
    tx = transcripts(annot, params.tsl_threshold)
    by_gene = dict(list(tx.groupby('gene_id')))
    with open(Path(a.run) / 'uniprot_canonicals.tsv') as fh:
        chosen = {r['gene_id']: r for r in csv.DictReader(fh, delimiter='\t')}
    records = output_records(a.run, set(chosen))

    s = [f'GTF {a.gtf}; run {a.run}', '']
    coding = tx.loc[tx['gene_type'] == 'protein_coding', 'gene_id'].drop_duplicates()
    s.append(f'GTF: {tx["gene_id"].nunique():,} genes, {len(coding):,} of type protein_coding')
    has_rule = tx.groupby('gene_id')['rule'].any()
    has_ec = tx[tx['ensembl_canonical']].groupby('gene_id').size()
    s.append(f'  protein-coding genes with a transcript the rule takes: {has_rule[coding].sum():,}; with an '
             f'{TAG} transcript: {coding.isin(has_ec.index).sum():,}')

    # ---- the genes with a GTF canonical by the rule: is it the Ensembl_canonical transcript's CDS? ----
    agree = collections.Counter()
    for gene_id in coding:
        g = by_gene[gene_id]
        taken = g[g['rule']]
        ec = g[g['ensembl_canonical']]
        if taken.empty:
            continue
        # the rule's choice: the longest CDS; among equals, the first transcript name
        best = taken.sort_values(['cds_length', 'transcript_name'], ascending=[False, True]).iloc[0]
        if ec.empty:
            agree['no Ensembl_canonical transcript'] += 1
        elif ec.iloc[0]['cds_structure'] == best['cds_structure']:
            agree['same coding sequence'] += 1
        elif ec.iloc[0]['rule']:
            agree['other: Ensembl_canonical passes the rule but is shorter'] += 1
        else:
            agree['other: Ensembl_canonical does not pass the rule'] += 1
    s.append('  protein-coding genes the rule gives a transcript, its choice against the Ensembl_canonical '
             'transcript: ' + ', '.join(f'{k} {v:,}' for k, v in agree.most_common()))
    s.append('')

    # ---- the run's genes without a canonical ----
    rows = []
    for gene_id, r in sorted(chosen.items()):
        g = by_gene.get(gene_id)
        if g is None or r['source'] != 'none':
            continue
        gene_type = g['gene_type'].iloc[0]
        ec = g[g['ensembl_canonical']]
        row = dict(gene_id=gene_id, gene_symbol=r['gene_symbol'], gene_type=gene_type, why=why_none(g),
                   **candidate(ec.iloc[0] if len(ec) else None),
                   readthrough='yes' if len(ec) and 'readthrough_transcript' in ec.iloc[0]['tags'] else 'no',
                   joined_records=records[gene_id]['joined'], orphan_records=records[gene_id]['orphan'])
        rows.append(row)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / 'genes.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    s.append(f'run: {len(chosen):,} genes; without a canonical: {len(rows):,}')
    for coding_gene in (True, False):
        sub = [r for r in rows if (r['gene_type'] == 'protein_coding') == coding_gene]
        s.append(f'{"protein-coding" if coding_gene else "other"} genes without a canonical: {len(sub):,}; '
                 f'their orphan records: {sum(r["orphan_records"] for r in sub):,}')
        for why, n in collections.Counter(r['why'] for r in sub).most_common():
            of = [r for r in sub if r['why'] == why]
            s.append(f'  {n:6,}  {why}; with an {TAG} transcript: {sum(1 for r in of if r["ec_transcript"]):,} '
                     f'(of type protein_coding with a CDS: {sum(1 for r in of if r["ec_type"] == "protein_coding" and r["ec_cds"]):,})')
        for rule, label in [('v1_complete', 'V1, Ensembl_canonical of type protein_coding with start and stop codons'),
                            ('v2_any_cds', 'V2, Ensembl_canonical with any CDS')]:
            taken = [r for r in sub if r[rule] == 'yes']
            s.append(f'  {label}: {len(taken):,} genes, {sum(r["orphan_records"] for r in taken):,} orphan records')
        for rt in ('yes', 'no'):
            of = [r for r in sub if r['readthrough'] == rt]
            s.append(f'  Ensembl_canonical tagged readthrough_transcript: {rt}: {len(of):,} genes '
                     f'({sum(1 for r in of if r["v1_complete"] == "yes"):,} with V1; MANE_Select '
                     f'{sum(1 for r in of if "MANE_Select" in r["ec_tags"]):,}), '
                     f'{sum(r["orphan_records"] for r in of):,} orphan records '
                     f'({sum(r["orphan_records"] for r in of if r["v1_complete"] == "yes"):,} of genes with V1)')
        ec = [r for r in sub if r['ec_transcript']]
        s.append('  their Ensembl_canonical transcripts: type ' +
                 str(dict(collections.Counter(r['ec_type'] for r in ec).most_common())) +
                 '; TSL ' + str(dict(sorted(collections.Counter(r['ec_tsl'] for r in ec).items()))) +
                 '; codons ' + str(dict(collections.Counter(r['ec_start_stop'] or 'neither' for r in ec
                                                           if r['ec_cds']).most_common())) +
                 '; MANE_Select ' + str(sum(1 for r in ec if 'MANE_Select' in r['ec_tags'])))
    # ---- the orphan records of the protein-coding genes without a canonical: the same event under another gene ----
    without = {r['gene_id']: r for r in rows if r['gene_type'] == 'protein_coding'}
    found = events(a.run)
    tally = collections.Counter()
    for key, genes_ in found.items():
        for gene_id, kind in genes_:
            if gene_id in without and kind == 'orphan':
                others = {k for g, k in genes_ if g != gene_id}
                where = 'joined under another gene' if 'joined' in others else \
                        'orphan under another gene' if others else 'under no other gene'
                tally[without[gene_id]['readthrough'], where] += 1
    s.append('')
    s.append('orphan records of protein-coding genes without a canonical, by whether the same event (splice type, '
             'slice, chromosome, anchor and alternative exon) has a record under another gene:')
    for (rt, where), n in sorted(tally.items()):
        s.append(f'  readthrough {rt}: {where}: {n:,}')

    # ---- the run's genes with a GTF canonical: the rule's choice against the Ensembl_canonical transcript ----
    gtf_genes = collections.Counter()
    for gene_id, r in chosen.items():
        if r['source'] != 'gtf' or gene_id not in by_gene:
            continue
        g = by_gene[gene_id]
        best = g[g['rule']].sort_values(['cds_length', 'transcript_name'], ascending=[False, True]).iloc[0]
        ec = g[g['ensembl_canonical']]
        gtf_genes['no Ensembl_canonical' if ec.empty else 'same coding sequence'
                  if ec.iloc[0]['cds_structure'] == best['cds_structure'] else 'other'] += 1
    s.append(f'the run\'s genes with a GTF canonical ({sum(gtf_genes.values()):,}): the rule\'s choice against the '
             f'Ensembl_canonical transcript: {dict(gtf_genes.most_common())}')
    if a.rmats:
        readthrough_events(a.rmats, annot, {g for g, r in without.items() if r['readthrough'] == 'yes'}, s)
    (out / 'summary.txt').write_text('\n'.join(s) + '\n')
    print('\n'.join(s))


def readthrough_events(rmats, annot, readthrough, s):
    """ rMATS SE events: unique across genes? those of the readthrough genes: is the skipped exon another gene's? """
    import pandas as pd
    key = ['chr', 'strand', 'exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
    se = pd.read_table(Path(rmats) / 'SE.MATS.JC.txt', usecols=['GeneID'] + key)
    ex = annot.loc[(annot['feature'] == 'exon') & (annot['gene_type'] == 'protein_coding') &
                   ~annot['gene_id'].isin(readthrough), ['seqname', 'strand', 'start', 'end']]
    exons = set(zip(ex['seqname'], ex['strand'], ex['start'], ex['end']))
    r = se[se['GeneID'].isin(readthrough)]
    hit = sum((c, st, b + 1, e) in exons for c, st, b, e in zip(r['chr'], r['strand'], r['exonStart_0base'],
                                                               r['exonEnd']))
    s.append(f'rMATS SE: {len(se):,} events, {len(se.drop_duplicates(key)):,} distinct by coordinates (an event is '
             f'reported under one gene only)')
    s.append(f'  SE events of the {len(readthrough):,} readthrough genes without a canonical: {len(r):,}; the skipped '
             f'exon is an annotated exon of another protein-coding gene: {hit:,} ({hit / max(len(r), 1):.0%})')


def compare(a):
    before, after = Path(a.before), Path(a.after)
    s = [f'{a.label}: {before} -> {after}', '']
    files = sorted({p.name for d in (before, after) for p in d.iterdir()} - {'logfile.log', 'model.pdf'})
    for f in files:
        b, n = before / f, after / f
        if not b.exists() or not n.exists() or b.read_bytes() != n.read_bytes():
            count = (lambda p: sum(1 for line in open(p) if line.startswith('>')) if p.exists() else 0) \
                if f.endswith('.fasta') else (lambda p: sum(1 for _ in open(p)) - 1 if p.exists() else 0)
            s.append(f'  {f}: {count(b):,} -> {count(n):,} {"records" if f.endswith(".fasta") else "rows"}')
    if len(s) == 2:
        s.append('  no file differs')

    sources = {}
    for label, d in (('before', before), ('after', after)):
        with open(d / 'uniprot_canonicals.tsv') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                sources.setdefault(r['gene_id'], {})[label] = r['source']
    changed = collections.Counter((v.get('before'), v.get('after')) for v in sources.values()
                                  if v.get('before') != v.get('after'))
    s.append(f'  canonical sources changed: {dict(changed)}')
    genes_changed = {g for g, v in sources.items() if v.get('before') != v.get('after')}
    rb, ra = output_records(before, genes_changed), output_records(after, genes_changed)
    totals = collections.Counter()
    for g in genes_changed:
        for k in ('joined', 'orphan'):
            totals[k, 'before'] += rb[g][k]
            totals[k, 'after'] += ra[g][k]
    s.append(f'  their records: joined {totals["joined", "before"]:,} -> {totals["joined", "after"]:,}; '
             f'orphan {totals["orphan", "before"]:,} -> {totals["orphan", "after"]:,}')
    if len(genes_changed) <= 30:
        symbols = {}
        with open(after / 'uniprot_canonicals.tsv') as fh:
            symbols = {r['gene_id']: r['gene_symbol'] for r in csv.DictReader(fh, delimiter='\t')}
        for g in sorted(genes_changed):
            s.append(f'    {g} {symbols.get(g, "")}: {sources[g].get("before")} -> {sources[g].get("after")}; '
                     f'joined {rb[g]["joined"]} -> {ra[g]["joined"]}, orphan {rb[g]["orphan"]} -> {ra[g]["orphan"]}')
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f'compare_{a.label}.txt').write_text('\n'.join(s) + '\n')
    print('\n'.join(s))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = p.add_subparsers(dest='part', required=True)
    g = sub.add_parser('genes')
    g.add_argument('--gtf', required=True)
    g.add_argument('--run', required=True)
    g.add_argument('--rmats', default=None)
    c = sub.add_parser('compare')
    c.add_argument('--before', required=True)
    c.add_argument('--after', required=True)
    c.add_argument('--label', required=True)
    for x in (g, c):
        x.add_argument('--out', default=str(Path(__file__).resolve().parent / 'genome'))
    a = p.parse_args()
    genes(a) if a.part == 'genes' else compare(a)
