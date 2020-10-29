"""
Run tests on inStrain quick_profile
"""

import os
import numpy as np
import shutil
from subprocess import call

import pandas as pd

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.profile.fasta
import tests.test_utils as test_utils

from test_utils import BTO

# class test_quickProfile:
#     def __init__(BTO):
#         BTO.test_dir = test_utils.load_random_test_dir()
#         BTO.bam1 = test_utils.load_data_loc() + \
#                    'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
#         BTO.fasta = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_min1000.fa'
#         BTO.script = test_utils.get_script_loc('quickProfile')
#         BTO.stb = test_utils.load_data_loc() + \
#                    'N5_271_010G1.maxbin2.stb'
# 
#     def setUp(BTO):
# 
#         if os.path.isdir(BTO.test_dir):
#             shutil.rmtree(BTO.test_dir)
#         os.mkdir(BTO.test_dir)
# 
#     def tearDown(BTO):
#         if os.path.isdir(BTO.test_dir):
#             shutil.rmtree(BTO.test_dir)
# 
#     def run(BTO, min=0, max=2, tests='all'):
#         if tests == 'all':
#             tests = np.arange(min, max + 1)
# 
#         for test_num in tests:
#             BTO.setUp()
#             print("\n*** Running {1} test {0} ***\n".format(test_num, BTO.__class__))
#             eval('BTO.test{0}()'.format(test_num))
#             BTO.tearDown()

def test_quick_profile_0(BTO):
    """
    Just run the program
    """
    # Run program
    base = BTO.test_dir + 'test'

    cmd = "inStrain quick_profile {1} {2} -o {3} -s {4}".format(True, BTO.bam1, BTO.fasta,
                                                                base, BTO.stb)
    print(cmd)
    call(cmd, shell=True)

    # Load output

    # Load true scaffolds
    tscaffs = []
    with open(BTO.stb, 'r') as o:
        for line in o.readlines():
            tscaffs.append(line.split()[0].strip())

    # Load reported scaffolds
    scaffs = inStrain.profile.fasta.load_scaff_list(base + '/scaffolds.txt')

    # Compare
    assert set(scaffs) == set(tscaffs)

def test_quick_profile_1(BTO):
    """
    Test the awk filter
    """
    # Run program
    base = BTO.test_dir + 'test'

    cmd = "inStrain quick_profile {1} {2} -o {3} -s {4} --stringent_breadth_cutoff {5}".format(
        True, BTO.bam1, BTO.fasta, base, BTO.stb, 0.5)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    Cdb = pd.read_csv(base + '/coverm_raw.tsv', sep='\t')
    Cdb['breadth'] = [x / y for x, y in zip(Cdb['Covered Bases'], Cdb['Length'])]
    assert abs(Cdb['breadth'].max() - 1) < 0.0001
    assert Cdb['breadth'].min() >= 0.5
    assert Cdb['breadth'].max() <= 1

def test_quick_profile_2(BTO):
    """
    Test no stb
    """
    # Run program
    base = BTO.test_dir + 'test'

    cmd = "inStrain quick_profile {1} {2} -o {3}".format(True, BTO.bam1, BTO.fasta,
                                                                base, BTO.stb)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    scaffs = inStrain.profile.fasta.load_scaff_list(base + '/scaffolds.txt')

    # Compare
    assert len(scaffs) == 178