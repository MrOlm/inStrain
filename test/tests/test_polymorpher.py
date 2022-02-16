"""
Run tests on inStrain compare
"""

import glob
import os
import shutil
from collections import defaultdict
from subprocess import call
import importlib
import logging
import pytest

import numpy as np
import pandas as pd
import pysam

import inStrain
import inStrain.SNVprofile
import inStrain.argumentParser
import inStrain.polymorpher
import inStrain.controller
import inStrain.profile.fasta
import inStrain.profile.profile_utilities
import inStrain.profile.snv_utilities
import inStrain.readComparer
import inStrain.compare_utils
import tests.test_utils as test_utils

from test_utils import BTO

def test_polymorpher_unit_0(BTO):
    """
    Test the ability to extract SNPs from a .bam file
    """
    bam_loc = BTO.bam1
    IS_loc = BTO.IS1

    # Load a list of positions to extract
    ISP = inStrain.SNVprofile.SNVprofile(IS_loc)
    SSdb = ISP.get_nonredundant_snv_table()

    # Load R2M
    R2M = ISP.get('Rdic')
    #print(ISP)

    for scaffold, Sdb in SSdb.groupby('scaffold'):
        Sdb = Sdb[(Sdb['scaffold'] == scaffold) & (Sdb['cryptic'] == False)][['scaffold', 'position', 'A', 'C', 'T', 'G', 'cryptic']]

        # Re_extract these
        position2counts = inStrain.polymorpher.extract_SNVS_from_bam(bam_loc, R2M[scaffold], Sdb['position'].tolist(), scaffold)

        # Compare
        for i, row in Sdb.iterrows():
            position = row['position']
            table_counts = [row['A'], row['C'], row['T'], row['G']]

            if position not in position2counts:
                print(scaffold)
                print(position)
                print(Sdb)

            calcd_counts = position2counts[position]

            if set(table_counts == calcd_counts) != set([True]):
                print(scaffold)
                print(row)
                print(table_counts)
                print(calcd_counts)

                assert False

def test_polymorpher_unit_1(BTO):
    """
    Test "extract_SNV_positions"
    """
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)

    tables = [S1.get('cumulative_snv_table').rename(
                columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                         'baseCoverage': 'position_coverage'}),
             S2.get('cumulative_snv_table').rename(
             columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                     'baseCoverage': 'position_coverage'})]

    names = [os.path.basename(S1.get('bam_loc')), os.path.basename(S2.get('bam_loc'))]

    scaffolds = set(tables[0]['scaffold']).intersection(set(tables[1]['scaffold']))

    for scaff in scaffolds:
        scaff_tables = [x[(x['scaffold'] == scaff) & (x['cryptic'] == False)] for x in tables]
        name2locs, all_SNVs = inStrain.polymorpher.extract_SNV_positions(scaff_tables, names)

        assert set(scaff_tables[0]['position']).union(set(scaff_tables[1]['position'])) \
               - set(scaff_tables[0]['position']).intersection(set(scaff_tables[1]['position'])) \
               == set(name2locs[names[0]]).union(set(name2locs[names[1]]))

def test_polymorpher_unit_2(BTO):
    """
    Test DDST and PST generation
    """
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)

    tables = [S1.get('cumulative_snv_table').rename(
        columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                 'baseCoverage': 'position_coverage'}),
        S2.get('cumulative_snv_table').rename(
            columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                     'baseCoverage': 'position_coverage'})]

    tables = [x.sort_values('mm')\
            .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
            .sort_index().drop(columns=['mm']) for x in tables]

    names = [os.path.basename(S1.get('bam_loc')), os.path.basename(S2.get('bam_loc'))]
    scaffolds = set(tables[0]['scaffold']).intersection(set(tables[1]['scaffold']))

    name2Rdic = {names[0]:S1.get('Rdic'),
                 names[1]: S2.get('Rdic')}
    name2bam_loc = {names[0]:BTO.bam1,
                  names[1]: BTO.bam2}

    for scaff in scaffolds:

        scaff_tables = [x[(x['scaffold'] == scaff) & (x['cryptic'] == False)] for x in tables]

        PM = inStrain.polymorpher.PoolController(scaff_tables, names, name2Rdic, name2bam_loc, scaff)
        PM.main()

        DSTdb = PM.DDST
        assert set(DSTdb.index.get_level_values(0)) == set(names)
        assert set(DSTdb.index.get_level_values(1)) == set(scaff_tables[0]['position']).union(set(scaff_tables[1]['position']))
        assert len(DSTdb) == (len(set(scaff_tables[0]['position']).union(set(scaff_tables[1]['position']))) * len(names))

        PMdb = PM.PST
        assert set(PMdb.index) == set(scaff_tables[0]['position']).union(set(scaff_tables[1]['position']))
        assert len(PMdb) == len(set(scaff_tables[0]['position']).union(set(scaff_tables[1]['position'])))