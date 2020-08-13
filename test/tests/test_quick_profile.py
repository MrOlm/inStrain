"""
Run tests on inStrain quick_profile
"""

import os
import shutil
from subprocess import call

import pandas as pd

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.profile.fasta
import tests.test_utils as test_utils


class test_quickProfile:
    def __init__(self):
        self.test_dir = test_utils.load_random_test_dir()
        self.bam = test_utils.load_data_loc() + \
                   'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa'
        self.script = test_utils.get_script_loc('quickProfile')
        self.stb = test_utils.load_data_loc() + \
                   'N5_271_010G1.maxbin2.stb'

    def setUp(self):

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test0()
        self.tearDown()

        self.setUp()
        self.test1()
        self.tearDown()

    def test0(self):
        """
        Just run the program
        """
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain quick_profile {1} {2} -o {3} -s {4}".format(self.script, self.bam, self.fasta,
                                                                    base, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output

        # Load true scaffolds
        tscaffs = []
        with open(self.stb, 'r') as o:
            for line in o.readlines():
                tscaffs.append(line.split()[0].strip())

        # Load reported scaffolds
        scaffs = inStrain.profile.fasta.load_scaff_list(base + '/scaffolds.txt')

        # Compare
        assert set(scaffs) == set(tscaffs)

    def test1(self):
        """
        Test the awk filter
        """
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain quick_profile {1} {2} -o {3} -s {4} --stringent_breadth_cutoff {5}".format(
            self.script, self.bam, self.fasta, base, self.stb, 0.5)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        Cdb = pd.read_csv(base + '/coverm_raw.tsv', sep='\t')
        Cdb['breadth'] = [x / y for x, y in zip(Cdb['Covered Bases'], Cdb['Length'])]
        assert abs(Cdb['breadth'].max() - 1) < 0.0001
        assert Cdb['breadth'].min() >= 0.5
        assert Cdb['breadth'].max() <= 1
