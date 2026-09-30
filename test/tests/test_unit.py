"""
Small unit tests that don't need any test data
"""

import os

import pandas as pd
from Bio.Seq import Seq

import inStrain.GeneProfile
import inStrain.SNVprofile
import inStrain.compare_controller


def test_multiallelic_snvs_get_genes():
    """
    SNVs with > 2 alleles should still be linked to their gene, using the highest mm
    """
    gene2sequence = {'g1': Seq("ATGAAACCCGGGTAA")}
    gdb = pd.DataFrame({'gene': ['g1'], 'scaffold': ['s'], 'start': [100], 'end': [114],
                        'direction': ['1'], 'partial': [False]})
    Ldb = pd.DataFrame({'scaffold': ['s'] * 4, 'position': [103, 103, 106, 200], 'mm': [1, 0, 0, 0],
                        'allele_count': [3, 2, 2, 2], 'con_base': ['A', 'A', 'C', 'A'],
                        'var_base': ['G', 'G', 'A', 'G'], 'ref_base': ['A', 'A', 'C', 'A']})

    sdb = inStrain.GeneProfile.Characterize_SNPs_wrapper(Ldb, gdb, gene2sequence)
    sdb = sdb.set_index('position')
    assert sdb.loc[103, 'allele_count'] == 3
    assert sdb.loc[103, 'gene'] == 'g1'
    assert sdb.loc[103, 'mutation'] == 'N:K3E'
    assert sdb.loc[200, 'mutation_type'] == 'I'

    GGdb, _ = inStrain.GeneProfile.calc_gene_snp_counts(gdb, Ldb, sdb.reset_index(), gene2sequence, scaffold='s')
    row = GGdb[GGdb['mm'] == 1].iloc[0]
    assert row['SNS_count'] + row['SNV_count'] == row['divergent_site_count'] == 2
    assert row['SNV_N_count'] == 2


def test_reorder_columns_is_stable():
    db = pd.DataFrame({'b': [1], 'z': [2], 'a': [3], 'y': [4], 'x': [5]})
    assert list(inStrain.SNVprofile.reorder_columns(db, ['a', 'b'])) == ['a', 'b', 'z', 'y', 'x']


def test_numeric_genome_names(tmp_path):
    loc = str(tmp_path / 'test.IS')
    IS = inStrain.SNVprofile.SNVprofile(loc)
    IS.store('genome_level_info', pd.DataFrame({'genome': ['243', 'genome_1'], 'breadth_minCov': [0.9, 0.9]}),
             'pandas', 'test')
    genomes = inStrain.SNVprofile.SNVprofile(loc).get('genome_level_info')['genome'].tolist()
    assert genomes == ['243', 'genome_1']


def test_compare_input_list(tmp_path):
    loc = str(tmp_path / 'inputs.txt')
    with open(loc, 'w') as o:
        o.write("a.IS\n\n# comment\nb.IS\n")
    assert inStrain.compare_controller.parse_input_list([loc, 'c.IS']) == ['a.IS', 'b.IS', 'c.IS']
