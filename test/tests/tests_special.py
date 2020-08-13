"""
Run tests on inStrain plotting
"""

import glob
import os
import shutil
from subprocess import call

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.irep_utilities
import inStrain.profile.fasta
import tests.test_utils as test_utils


class test_special:
    def __init__(self):
        self.test_dir = test_utils.load_random_test_dir()
        self.sorted_bam = test_utils.load_data_loc() + \
                          'N5_271_010G1_scaffold_963_Ns.fasta.sorted.bam'
        self.fasta = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_963_Ns.fasta'
        self.bbsam = test_utils.load_data_loc() + \
                     'bbmap_N5_271_010G1_scaffold_963.fasta.sorted.bam'
        self.bbfasta = test_utils.load_data_loc() + \
                       'N5_271_010G1_scaffold_963.fasta'
        self.IS = test_utils.load_data_loc() + \
                  'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
        self.small_fasta = test_utils.load_data_loc() + \
                           'SmallScaffold.fa'
        self.small_bam = test_utils.load_data_loc() + \
                         'SmallScaffold.fa.sorted.bam'
        self.genes = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.stb = test_utils.load_data_loc() + \
                   'GenomeCoverages.stb'

    def setUp(self):
        self.tearDown()

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        # self.setUp()
        # self.test1()
        # self.tearDown()
        #
        # self.setUp()
        # self.test2()
        # self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()
        #
        # self.setUp()
        # self.test4()
        # self.tearDown()

    def test1(self):
        """
        Yell if the minor version is not correct
        """
        assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.3')

        assert not inStrain.SNVprofile.same_versions('1.0.0', '0.8.3')

        assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.4')

        assert not inStrain.SNVprofile.same_versions('0.9.0', '0.8.4')

        assert not inStrain.SNVprofile.same_versions('1.1.0', '0.1.0')

    def test2(self):
        """
        Make sure it works with Ns in the input
        """
        # Run base
        base = self.test_dir + 'testNs'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(
            self.sorted_bam,
            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object

    def test3(self):
        """
        Try other mappers
        """
        # Run bbmap and let it crash
        base = self.test_dir + 'bbMap'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(self.bbsam,
                                                                                                         self.bbfasta,
                                                                                                         base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()
        assert len(MCLdb) == 0

        # Run bbmap and make it work
        base = self.test_dir + 'bbMap'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation --use_full_fasta_header".format(
            self.bbsam,
            self.bbfasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()
        assert len(MCLdb) > 0

    def test4(self):
        """
        Test other --run_statistics
        """
        # Super quick profile run
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        cmd = "inStrain profile {0} {1} -o {2} -s {3} -g {4} --skip_plot_generation".format(
            self.small_bam, self.small_fasta, location, self.stb, self.genes)
        print(cmd)
        call(cmd, shell=True)

        cmd = "inStrain other --run_statistics {0} --debug".format(location)
        print(cmd)
        call(cmd, shell=True)
        assert len(glob.glob(location + '/log/*')) == 3, location
