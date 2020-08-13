"""
Run tests on inStrain plotting
"""

import glob
import importlib
import logging
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


class test_plot:
    def __init__(self):
        self.test_dir = test_utils.load_random_test_dir()
        self.IS = test_utils.load_data_loc() + \
                  'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.ForPlotting.IS'
        self.RC_Loc = test_utils.load_data_loc() + \
                      'readComparer_vCurrent.RC'
        self.stb = test_utils.load_data_loc() + \
                   'N5_271_010G1.maxbin2.stb'
        self.genes = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.sorted_bam = test_utils.load_data_loc() + \
                          'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa'

    def setUp(self):

        self.tearDown()
        importlib.reload(logging)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test1(view=False)
        self.tearDown()

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

    def test1(self, view=False):
        """
        Make sure all plots are made
        """
        # Run the IS plots
        FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
                'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
                'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
                'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
                'GeneHistogram_plot.pdf']

        # !! Maybe when stable you can do this; for now you need to re-run it !!
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)
        for f in glob.glob(location + '/log/*'):
            os.remove(f)

        # Run program
        # location = self.test_dir + 'test'
        # cmd = "inStrain profile {1} {2} -o {3} -g {4} -s {5} --skip_plot_generation -p 6 -d".format(
        #         'junk', self.sorted_bam, self.fasta, location, self.genes, self.stb)
        # print(cmd)
        # call(cmd, shell=True)
        # os.remove(location + '/log/log.log')

        cmd = "inStrain plot -i {0} -d".format(location)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        figs = glob.glob(IS.get_location('figures') + '*')

        # Make sure all figures are made
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000, fig

        # Make sure logging works
        rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
        got = False
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'Plot' in line:
                    got += 1
                # print(line)
        assert got == 11, got  # Its in there twice for random reasons

        if view:
            assert False, "Get on the figures here: " + IS.get_location('figures')

    def test2(self):
        """
        Make sure all RC plots are made
        """
        # Run the IS plots
        FIGS = ['inStrainCompare_dendrograms.pdf']

        location = os.path.join(self.test_dir, os.path.basename(self.RC_Loc))
        shutil.copytree(self.RC_Loc, location)

        cmd = "inStrain genome_wide -i {0} -s {1}".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        cmd = "inStrain plot -i {0} -d".format(location)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        figs = glob.glob(IS.get_location('figures') + '*')

        # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000

    def test3(self):
        """
        Test the breadth cutoff
        """
        # Run the IS plots
        FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
                'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
                'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
                'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
                'GeneHistogram_plot.pdf']

        # FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf']#, 'genomeWide_microdiveristy_metrics.pdf',
        # #         'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
        # #         'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
        # #         'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
        # #         'GeneHistogram_plot.pdf']

        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        cmd = "inStrain plot -i {0} -mb 0.97".format(location)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        figs = glob.glob(IS.get_location('figures') + '*')

        # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000

    def test4(self):
        """
        Test the genome name
        """
        # Run the IS plots
        FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
                'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
                'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
                'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
                'GeneHistogram_plot.pdf']

        # FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf']#, 'genomeWide_microdiveristy_metrics.pdf',
        # #         'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
        # #         'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
        # #         'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
        # #         'GeneHistogram_plot.pdf']

        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        cmd = "inStrain plot -i {0} -g maxbin2.maxbin.001.fasta".format(location)
        # cmd = "inStrain plot -i {0} -g poop".format(location)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        figs = glob.glob(IS.get_location('figures') + '*')

        # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000
