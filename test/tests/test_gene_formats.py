"""
Tests for gene files that aren't prodigal .fna files (GFF3 and genbank)

The prodigal genes in the test data are converted to the other formats, so the results can be compared
directly against the prodigal results
"""

import glob
import importlib
import logging
import os
import shutil

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import inStrain
import inStrain.argumentParser
import inStrain.controller
import inStrain.GeneProfile
import inStrain.SNVprofile
from tests.test_utils import BTO


def _prodigal_genes(genes_loc):
    """
    Return a list of (gene, scaffold, start, end, strand) from a prodigal .fna file (1-based coordinates)
    """
    genes = []
    for record in SeqIO.parse(genes_loc, 'fasta'):
        fields = record.description.split('#')
        gene = record.id
        scaffold = '_'.join(gene.split('_')[:-1])
        genes.append((gene, scaffold, int(fields[1]), int(fields[2]), '+' if fields[3].strip() == '1' else '-'))
    return genes


def write_gff(genes_loc, fasta_loc, out_loc, embed_fasta=False):
    """
    Write the prodigal genes as a GFF3 file, in the style of Bakta / Prokka (gene + CDS lines, URL-encoded
    attributes)
    """
    with open(out_loc, 'w') as o:
        o.write('##gff-version 3\n')
        for gene, scaffold, start, end, strand in _prodigal_genes(genes_loc):
            o.write(f'{scaffold}\tProdigal\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene}_gene;locus_tag={gene}\n')
            o.write(f'{scaffold}\tProdigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t'
                    f'ID={gene};Parent={gene}_gene;locus_tag={gene};product=hypothetical protein%3B putative\n')
        if embed_fasta:
            o.write('##FASTA\n')
            with open(fasta_loc) as f:
                o.write(f.read())


def write_genbank(genes_loc, fasta_loc, out_loc):
    """
    Write the prodigal genes as a genbank file, naming genes with the locus_tag qualifier
    """
    scaff2genes = {}
    for gene, scaffold, start, end, strand in _prodigal_genes(genes_loc):
        scaff2genes.setdefault(scaffold, []).append((gene, start, end, strand))

    records = []
    for record in SeqIO.parse(fasta_loc, 'fasta'):
        if record.id not in scaff2genes:
            continue
        rec = SeqRecord(Seq(str(record.seq)), id=record.id, name=record.id[:16], description='test')
        rec.annotations['molecule_type'] = 'DNA'
        for gene, start, end, strand in scaff2genes[record.id]:
            rec.features.append(SeqFeature(FeatureLocation(start - 1, end, strand=1 if strand == '+' else -1),
                                           type='CDS', qualifiers={'locus_tag': [gene]}))
        records.append(rec)
    SeqIO.write(records, out_loc, 'genbank')


def _combined(scaff2geneinfo):
    db = pd.concat(scaff2geneinfo.values()).reset_index(drop=True)
    return db.sort_values('gene').reset_index(drop=True)[['gene', 'scaffold', 'direction', 'partial', 'start', 'end']]


def _assert_same_genes(parsed, expected):
    gdb, seqs = parsed
    egdb, eseqs = expected

    a = _combined(gdb)
    b = _combined(egdb)
    b['partial'] = b['partial'].astype(bool)
    a['partial'] = a['partial'].astype(bool)
    pd.testing.assert_frame_equal(a, b, check_dtype=False)

    for scaff, g2s in eseqs.items():
        for gene, seq in g2s.items():
            assert str(seqs[scaff][gene]).upper() == str(seq).upper(), gene


def test_gff_matches_prodigal(BTO):
    expected = inStrain.GeneProfile.parse_genes(BTO.genes)

    # Sequences from the .fasta file
    gff = os.path.join(BTO.test_dir, 'genes.gff')
    write_gff(BTO.genes, BTO.fasta, gff)
    _assert_same_genes(inStrain.GeneProfile.parse_genes(gff, fasta=BTO.fasta), expected)

    # Sequences from a ##FASTA section in the GFF
    gff3 = os.path.join(BTO.test_dir, 'genes_with_fasta.gff3')
    write_gff(BTO.genes, BTO.fasta, gff3, embed_fasta=True)
    _assert_same_genes(inStrain.GeneProfile.parse_genes(gff3), expected)

    # No way to get sequences
    with pytest.raises(Exception):
        inStrain.GeneProfile.parse_genes(gff)


def test_gff_edge_cases(BTO):
    gff = os.path.join(BTO.test_dir, 'edge.gff')
    seq = 'ATGAAACCCGGGTTTTAA' * 20
    with open(gff, 'w') as o:
        o.write('##gff-version 3\n')
        o.write('s1\tx\tCDS\t1\t18\t.\t+\t0\tID=cds-1;locus_tag=A_1;partial=true\n')
        o.write('s1\tx\tCDS\t19\t36\t.\t-\t0\tlocus_tag=A_2\n')
        o.write('s1\tx\tCDS\t40\t60\t.\t+\t0\tID=split\n')
        o.write('s1\tx\tCDS\t62\t80\t.\t+\t0\tID=split\n')
        o.write('s1\tx\ttRNA\t100\t170\t.\t+\t.\tID=trna\n')
        o.write('missing_scaffold\tx\tCDS\t1\t18\t.\t+\t0\tID=gone\n')
        o.write(f'##FASTA\n>s1 some description\n{seq}\n')

    gdb, seqs = inStrain.GeneProfile.parse_genes(gff)
    db = _combined(gdb).set_index('gene')

    # ID is used when present, then locus_tag; split CDSs, other features and unknown scaffolds are skipped
    assert sorted(db.index) == ['A_2', 'cds-1']
    assert db.loc['cds-1', 'partial'] == True
    assert db.loc['A_2', 'partial'] == False
    assert db.loc['cds-1', 'direction'] == '1'
    assert db.loc['A_2', 'direction'] == '-1'
    assert (db.loc['A_2', 'start'], db.loc['A_2', 'end']) == (18, 35)
    assert str(seqs['s1']['cds-1']) == seq[0:18]
    assert str(seqs['s1']['A_2']) == str(Seq(seq[18:36]).reverse_complement())


def test_genbank_matches_prodigal(BTO):
    """
    Genes on both strands, named by locus_tag; strands must come out as '1' / '-1' like prodigal
    """
    expected = inStrain.GeneProfile.parse_genes(BTO.genes)
    gbk = os.path.join(BTO.test_dir, 'genes.gbk')
    write_genbank(BTO.genes, BTO.fasta, gbk)
    _assert_same_genes(inStrain.GeneProfile.parse_genes(gbk), expected)


def _profile_genes(BTO, gene_file, name):
    importlib.reload(logging)
    location = os.path.join(BTO.test_dir, name)
    shutil.copytree(BTO.IS_nogenes, location)
    cmd = f"profile_genes -i {location} -g {gene_file}"
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split()))
    return inStrain.SNVprofile.SNVprofile(location)


@pytest.mark.parametrize('fmt', ['gff', 'gbk'])
def test_profile_genes_with_other_formats(BTO, fmt):
    """
    Running profile_genes with a GFF3 / genbank file gives the same results as the prodigal file
    """
    if fmt == 'gff':
        gene_file = os.path.join(BTO.test_dir, 'genes.gff')
        write_gff(BTO.genes, BTO.fasta, gene_file, embed_fasta=True)
    else:
        gene_file = os.path.join(BTO.test_dir, 'genes.gbk')
        write_genbank(BTO.genes, BTO.fasta, gene_file)

    expected = _profile_genes(BTO, BTO.genes, 'prodigal.IS')
    got = _profile_genes(BTO, gene_file, f'{fmt}.IS')

    for thing in ['genes_table', 'genes_coverage', 'genes_clonality', 'genes_SNP_count', 'SNP_mutation_types']:
        e = expected.get(thing)
        g = got.get(thing)
        sort_cols = [c for c in ['gene', 'scaffold', 'position', 'mm'] if c in e.columns]
        e = e.sort_values(sort_cols).reset_index(drop=True)
        g = g.sort_values(sort_cols).reset_index(drop=True)[list(e.columns)]
        if 'partial' in e.columns:
            e['partial'] = e['partial'].astype(bool)
            g['partial'] = g['partial'].astype(bool)
        pd.testing.assert_frame_equal(e, g, check_dtype=False, obj=thing)
