#!/usr/bin/env python
'''
Run tests
'''

import os
import sys
import glob
import shutil
import pickle
import warnings
import logging
import importlib
import numpy as np
warnings.filterwarnings("ignore")
import pandas as pd
from subprocess import call
from Bio import SeqIO
from collections import defaultdict

#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import inStrain
import inStrain.SNVprofile
import inStrain.controller
import inStrain.deprecated
import inStrain.filter_reads
import inStrain.deprecated_filter_reads
import inStrain.profileUtilities
import inStrain.readComparer
import inStrain.argumentParser
import inStrain.logUtils
import inStrain.irep_utilities
import inStrain.genomeUtilities

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)

# Name translations going from version 1.2 to 1.3
twelve2thirteen = {
"baseCoverage":"position_coverage",
"refBase":"ref_base",
"conBase":'con_base',
"varBase":'var_base',
"varFreq":"var_freq",
"conFreq":"con_freq",
'refFreq':'ref_freq',
'coverage_median':'coverage_median',
'std_cov':'coverage_std',
'microdiversity':'nucl_diversity',
'mean_microdiversity':'nucl_diversity',
'median_microdiversity':'nucl_diversity_median',
'masked_breadth':'breadth_minCov',
'unmaskedBreadth':'breadth_minCov',
'rarefied_breadth':'breadth_rarefied',
'expected_breadth':'breadth_expected',
'rarefied_mean_microdiversity':'nucl_diversity_rarefied',
'rarefied_median_microdiversity':'nucl_diversity_rarefied_median',
'median_cov':'coverage_median',
'Reference_SNPs':'SNS_count',
'conANI':'conANI_reference',
'popANI':'popANI_reference',
'SNPs':'divergent_site_count',
'true_length':'length',
'total_SNPs':'divergent_site_count',
'population_SNPs':'population_divergent_sites',
'consensus_SNPs':'consensus_divergent_sites',
'filtered_pairs':'filtered_read_pair_count',
'SNPs':'divergent_site_count',
'SNS_S_count':'',
}

# Removed from version 1.3
del_thirteen = {'bases_w_0_coverage',
                'mean_clonality',
                'median_clonality',
                'BiAllelic_SNPs',
                'SNPs_per_bp',
                'MultiAllelic_SNPs',
                'clonality',
                'min_ANI'
                }

new_thirteen = {'linked_SNV_count',
                'SNV_distance_mean',
                'coverage_SEM',
                'r2_mean',
                'iRep_GC_corrected',
                'iRep',
                'd_prime_mean',
                'coverage_SEM',
                'SNV_count',
                'SNS_S_count',
                'S_sites',
                'SNV_N_count',
                'SNS_N_count',
                'N_sites',
                'dNdS_substitutions',
                'SNV_S_count',
                'pNpS_variants',
                'gene_length'
                }

def load_data_loc():
    return os.path.join(str(os.getcwd()), \
        'test_data/')

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'test_backend/testdir/')
    return loc

def get_script_loc(script):
    if script == 'inStrain':
        return os.path.join(str(os.getcwd()), \
            '../bin/inStrain')
    if script == 'filter_reads':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/filter_reads.py')
    if script == 'deprecated_gene_statistics':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/deprecated_gene_statistics.py')
    if script == 'SNVprofile':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/SNVprofile.py')
    if script == 'readcomparer':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/readComparer.py')
    if script == 'quickProfile':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/quickProfile.py')
    if script == 'GeneProfile':
        return os.path.join(str(os.getcwd()), \
            '../inStrain/GeneProfile.py')

class test_plot():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.ForPlotting.IS'
        self.RC_Loc = load_data_loc() + \
            'readComparer_vCurrent.RC'
        self.stb = load_data_loc() + \
            'N5_271_010G1.maxbin2.stb'
        self.genes = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa.genes.fna'

        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'

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
        '''
        Make sure all plots are made
        '''
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
        log_log = IS.get_location('log') + 'log.log'
        rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
        got = False
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'Plot' in line:
                    got += 1
                #print(line)
        assert got == 11, got # Its in there twice for random reasons

        if view:
            assert False, "Get on the figures here: " + IS.get_location('figures')

    def test2(self):
        '''
        Make sure all RC plots are made
        '''
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

        #assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000

    def test3(self):
        '''
        Test the breadth cutoff
        '''
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

        #assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000

    def test4(self):
        '''
        Test the genome name
        '''
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
        #cmd = "inStrain plot -i {0} -g poop".format(location)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        figs = glob.glob(IS.get_location('figures') + '*')

        #assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
        for F in FIGS:
            assert len([f for f in figs if F in f]) == 1, F
        for fig in figs:
            assert os.path.getsize(fig) > 1000

class test_genome_wide():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.3.0g/'
        self.stb = load_data_loc() + \
            'N5_271_010G1.maxbin2.stb'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.single_scaff = load_data_loc() + \
            'N5_271_010G1_scaffold_101.fasta'
        self.extra_single_scaff = load_data_loc() + \
            'N5_271_010G1_scaffold_101_extra.fasta'
        self.RC_Loc = load_data_loc() + \
            'readComparer_vCurrent.RC'
        self.iRep_test_set = load_data_loc() + \
            'test_iRep.p'

        self.fasta_extra = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000_extra.fa'
        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.stb2 = load_data_loc() + \
            'GenomeCoverages.stb'


        self.tearDown()

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

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

        self.setUp()
        self.test5()
        self.tearDown()

        self.setUp()
        self.test6()
        self.tearDown()

    def test0(self):
        '''
        Just run the program
        '''
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        cmd = "inStrain genome_wide -i {0} -s {1} -d".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if (('genome_info' in f) & ('testdir_' in f))]
        assert len(files) == 1, files

        # Make sure mm level is False in the output
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert 'mm' not in db.columns
            print(db.head())

        # Make sure mm level is True in the raw data
        db = IS.get('genome_level_info')
        assert 'mm' in db.columns
        #print(IS)

    def test1(self):
        '''
        Test the different ways of importing the .stb file
        '''
        # Setup
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)
        for f in glob.glob(location + '/output/*.tsv'):
            os.remove(f)

        # Run with a normal .stb
        cmd = "inStrain genome_wide -i {0} -s {1}".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1, [len(files), files]
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert len(db['genome'].unique()) == 2, [f, db]
            assert len(db[~db['iRep_GC_corrected'].isna()]) > 0, db

        # Run with a skip mm .stb
        cmd = "inStrain genome_wide -i {0} -s {1} --skip_mm_profiling".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1, [len(files), files]
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert len(db['genome'].unique()) == 2, [f, db]
            assert len(db[~db['iRep_GC_corrected'].isna()]) > 0, db

        # Run with a no .stb
        cmd = "inStrain genome_wide -i {0}".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert len(db['genome'].unique()) == 1, [f, db]

        # Run with a single .fasta
        cmd = "inStrain genome_wide -i {0} -s {1}".format(location, self.single_scaff)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert len(db['genome'].unique()) == 1

        # Run with a bunch of .fastas
        cmd = "inStrain genome_wide -i {0} -s {1} {2} {3}".format(location, self.single_scaff,
                self.fasta, self.extra_single_scaff)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert len(db['genome'].unique()) == 2, [f, db]
            assert 'mm' not in db.columns

    def test2(self):
        '''
        Test running on an RC object
        '''
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.RC_Loc))
        shutil.copytree(self.RC_Loc, location)

        cmd = "inStrain genome_wide -i {0} -s {1}".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genomeWide' in f]
        assert len(files) == 1

    def test3(self):
        '''
        Test running on the mm level
        '''
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)
        for f in glob.glob(location + '/output/*.tsv'):
            os.remove(f)

        cmd = "inStrain genome_wide -i {0} -s {1} --mm_level".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1, files
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert 'mm' in list(db.columns)

    def test4(self):
        '''
        Test running RC object on the mm level
        '''
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.RC_Loc))
        shutil.copytree(self.RC_Loc, location)
        for f in glob.glob(location + '/output/*.tsv'):
            os.remove(f)

        cmd = "inStrain genome_wide -i {0} -s {1} --mm_level".format(location, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genomeWide' in f]
        assert len(files) == 1
        for f in files:
            db = pd.read_csv(f, sep='\t')
            assert 'mm' in list(db.columns)

    def test5(self):
        '''
        Test iRep methods
        '''
        # Make sure it can run on a random thing
        iRep, iRep_junk = inStrain.irep_utilities.calculate_iRep_from_coverage_array(
                            np.array([5]*50000 + [6]*50000 + [7]*50000), 1)
        assert abs(iRep - 1.62) < 0.01

        # Load pickled test data
        Test_sets = pickle.load( open( self.iRep_test_set, "rb" ) )

        # Make sure they all work
        for t_set in Test_sets:
            order, length, gc, windows,  OLTwidnows, LTwidnows, raw_iRep, iRep = t_set

            # Make windows into df
            Idb = pd.DataFrame({'index':windows[0], 'coverage':windows[1]})
            GCdb = pd.DataFrame({'index':gc[0], 'GC_content':gc[1]})
            Idb = pd.merge(Idb, GCdb, on='index')

            # Using Chris' array, filter
            Idb = inStrain.irep_utilities._iRep_filter_windows(Idb, on='coverage')
            Idb = inStrain.irep_utilities._iRep_gc_bias(Idb)

            # Log transform
            Idb['coverage_OLT'] = inStrain.irep_utilities._iRep_log_transform(Idb['coverage'])
            Idb['coverage_LT'] = inStrain.irep_utilities._iRep_log_transform(Idb['corrected_coverage'])

            assert len(LTwidnows[1]) == len(OLTwidnows[1])

            assert len(Idb['coverage_OLT'].tolist()) == len(OLTwidnows[1])
            for x, y in zip(Idb['coverage_OLT'].tolist(), OLTwidnows[1]):
                if abs(x - y) > 0.0001:
                    assert False, [x, y]

            assert len(Idb['coverage_LT'].tolist()) == len(LTwidnows[1])
            for x, y in zip(Idb['coverage_LT'].tolist(), LTwidnows[1]):
                if abs(x - y) > 0.0001:
                    assert False, [x, y]

            # Get raw iRep values
            r = inStrain.irep_utilities._calc_iRep(Idb, length, on='coverage_OLT')
            i = inStrain.irep_utilities._calc_iRep(Idb, length, on='coverage_LT')

            assert abs(r - raw_iRep) < 0.0001, "raw is diff! {0} {1}".format(r, raw_iRep)
            assert abs(i - iRep) < 0.0001, print("iRep is diff! {0} {1}".format(i, iRep))

    def test6(self):
        '''
        Test running on database mode (there is a genome with no reads)
        '''
        location = os.path.join(self.test_dir, os.path.basename(self.IS))

        # Make sure empty scaffolds don't mess it up
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 -s {4} --skip_plot_generation".format(
            'junk', self.sorted_bam, self.fasta_extra, location, self.stb2)
        print(cmd)
        exit_code = call(cmd, shell=True)
        assert exit_code == 0; exit_code

        print(cmd)
        call(cmd, shell=True)

        # Load output
        IS = inStrain.SNVprofile.SNVprofile(location)
        files = glob.glob(IS.get_location('output') + '*')
        files = [f for f in files if 'genome_info' in f]
        assert len(files) == 1, files
        for f in files:
            db = pd.read_csv(f, sep='\t')
            print(len(db))

class test_quickProfile():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.script = get_script_loc('quickProfile')
        self.stb = load_data_loc() + \
            'N5_271_010G1.maxbin2.stb'

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
        '''
        Just run the program
        '''
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain quick_profile {1} {2} -o {3} -s {4}".format(self.script, self.bam, self.fasta, \
            base, self.stb)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        Cdb = pd.read_csv(base + '/genomeCoverage.csv')

        # Load true scaffolds
        tscaffs = []
        with open(self.stb, 'r') as o:
            for line in o.readlines():
                tscaffs.append(line.split()[0].strip())

        # Load reported scaffolds
        scaffs = inStrain.controller.load_scaff_list(base + '/scaffolds.txt')

        # Compare
        assert set(scaffs) == set(tscaffs)

    def test1(self):
        '''
        Test the awk filter
        '''
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain quick_profile {1} {2} -o {3} -s {4} --stringent_breadth_cutoff {5}".format(
            self.script, self.bam, self.fasta, base, self.stb, 0.5)
        print(cmd)
        call(cmd, shell=True)

        # Load output
        Cdb = pd.read_csv(base + '/coverm_raw.tsv', sep='\t')
        Cdb['breadth'] = [x/y for x, y in zip(Cdb['Covered Bases'], Cdb['Length'])]
        assert abs(Cdb['breadth'].max() - 1) < 0.0001
        assert Cdb['breadth'].min() >= 0.5
        assert Cdb['breadth'].max() <= 1

class test_SNVprofile():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.v3.IS.pickle'
        self.script = get_script_loc('SNVprofile')

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

        importlib.reload(logging)

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
        '''
        Test basic things
        '''
        # Make an object
        location = self.test_dir + 'test.IS'
        s = inStrain.SNVprofile.SNVprofile(location)
        assert s.get('location') == os.path.abspath(location)

        # Test overwriting a value
        s.store('version', '1.0.0.testeroni', 'value', 'Version of inStrain')
        assert s.get('version') == '1.0.0.testeroni'

        # Test adding a new value
        s.store('testeroni', 'testeroni', 'value', 'Description of testeroni')

        # Load a new version of this and make sure changes are saved
        s2 = inStrain.SNVprofile.SNVprofile(location)
        assert s2.get('testeroni') == 'testeroni'

        # Make sure the README is there
        assert os.path.exists(location + '/raw_data/_README.txt')

    def test1(self):
        '''
        Test changing from old to new version
        '''
        # Make a copy of this file
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copyfile(self.IS, location)

        # Run the command
        cmd = "inStrain other --old_IS {1}".format(self.script, location)
        print(cmd)
        call(cmd, shell=True)

        # Load the old version
        oIS = inStrain.SNVprofile.SNVprofile_old()
        oIS.load(location)

        # Load the new version
        nIS = inStrain.SNVprofile.SNVprofile(location.replace('.pickle', ''))

        # Compare
        Adb = nIS._get_attributes_file()
        assert len(Adb) > 2
        assert len(Adb) == 16, len(Adb)
        for name, row in Adb.iterrows():

            # Translate names
            oldname = name
            if name == 'location':
                oldname = 'file'

            # Skip some
            if name in ['location', 'version', 'fasta_loc', 'bam_loc', 'old_version']:
                continue

            # Make sure they're the same
            if row['type'] == 'numpy':
                new = nIS.get(name)
                old = getattr(oIS, name)
                assert type(new) == type(old), [type(new), type(old)]
                assert len(new) == len(old)
                for i in range(0, len(new)):
                    assert (new[i] == old[i]).all()

            if row['type'] == 'special':
                if name in ['covT', 'snpsCounted']:
                    new = nIS.get(name)
                    old = getattr(oIS, name)
                    assert type(new) == type(old), [type(new), type(old)]

                    old_scaffs = [x for x, i in old.items() if len(i) > 0]
                    new_scaffs = [x for x, i in new.items() if len(i) > 0]

                    assert set(old_scaffs) == set(new_scaffs), [new, old]

            elif row['type'] == 'pandas':
                new = nIS.get(name)
                old = getattr(oIS, name)
                print(new.head())
                print(old.head())

            else:
                nIS.get(name) == getattr(oIS, name), name

class test_gene_statistics():
    def setUp(self):
        self.script = get_script_loc('GeneProfile')
        self.old_script = get_script_loc('deprecated_gene_statistics')
        self.test_dir = load_random_test_dir()
        self.old_IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.v4.IS'
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.3.0g'
        self.genes = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.IS2 = load_data_loc() + \
            'sars_cov_2_MT039887.1.fasta.bt2-vs-SRR11140750.sam.IS'
        self.genbank = load_data_loc() + \
            'sars_cov_2_MT039887.1.gb'

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

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

    def test0(self):
        '''
        Test deprecated gene_statistic
        '''
        location = os.path.join(self.test_dir, os.path.basename(self.old_IS))
        shutil.copytree(self.old_IS, location)

        # Run program
        base = self.test_dir + 'test'

        cmd = "{0} {1} -g {3}".format(self.old_script, location, \
            base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert len(glob.glob(location + '/output/*')) >= 6

        # Read something
        Rdb = pd.read_csv(location + '/output/aa-SNVs.tsv', sep='\t')
        print(len(Rdb))
        #assert len(Rdb) == 31
        assert len(Rdb) == 350

        # Check gene coverages [THIS DOESN'T WORK]
        # Gdb = pd.read_csv(location + '/output/genes.tsv', sep='\t')
        # print(len(Gdb))
        # assert len(Gdb) > 0

    def test1(self):
        '''
        Test GeneProfile vs deprecated gene_statistic
        '''
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.old_IS))
        shutil.copytree(self.old_IS, location)

        # Run program
        base = self.test_dir + 'testO'

        cmd = "{0} {1} -g {3}".format(self.old_script, location, \
            base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Read the mutations
        Rdb = pd.read_csv(location + '/output/aa-SNVs.tsv', sep='\t')
        assert len(Rdb) == 350

        # Get the mutations in this IS (OLD VERSION)
        Sdb =  inStrain.SNVprofile.SNVprofile(location).get_nonredundant_snv_table()
        Sdb['morphia'] = Sdb['morphia'].astype(int)
        Sdb = Sdb[Sdb['cryptic'] == False]
        if 'morphia' in Sdb.columns:
            Sdb = Sdb[Sdb['morphia'] == 2]
            Sdb = Sdb.drop(columns="morphia")
        Sdb = Sdb.drop(columns="cryptic")
        Sdb['position'] = Sdb['position'].astype(int)
        SOdb = Sdb

        assert len(Rdb) == len(SOdb)

        # Run my version
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3}".format(self.script, location, \
            base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Get the mutations in this IS
        Sdb =  inStrain.SNVprofile.SNVprofile(location).get_nonredundant_snv_table()
        Sdb['allele_count'] = Sdb['allele_count'].astype(int)
        Sdb = Sdb[Sdb['cryptic'] == False]
        if 'allele_count' in Sdb.columns:
            Sdb = Sdb[Sdb['allele_count'].isin([1,2])]
            Sdb = Sdb.drop(columns="allele_count")
        Sdb = Sdb.drop(columns="cryptic")
        Sdb['position'] = Sdb['position'].astype(int)
        SNdb = Sdb
        SNdb['mut_key'] = ["{0}:{1}".format(x, y) for x, y in zip(SNdb['scaffold'], SNdb['position'])]

        IS = inStrain.SNVprofile.SNVprofile(location)
        RNdb = IS.get('SNP_mutation_types')
        RNdb['mut_key'] = ["{0}:{1}".format(x, y) for x, y in zip(RNdb['scaffold'], RNdb['position'])]

        # Filter out the mutaitons from scaffolds without and genes
        gdb = IS.get('genes_table')
        SNdb = SNdb[SNdb['scaffold'].isin(list(gdb['scaffold'].unique()))]
        Rdb = Rdb[Rdb['scaffold'].isin(list(gdb['scaffold'].unique()))]

        # Make sure you characterized all SNPs
        assert len(RNdb) == len(SNdb), (len(RNdb), len(SNdb))

        # Filter to only those mutations in both lists
        Rdb['mut_key'] = ["{0}:{1}".format(x, y) for x, y in zip(Rdb['scaffold'], Rdb['position'])]
        RNdb['mut_key'] = ["{0}:{1}".format(x, y) for x, y in zip(RNdb['scaffold'], RNdb['position'])]

        # Filter out SNPs that are in multiple genes, which CCs script doesnt handle
        RNdb = RNdb[RNdb['mutation_type'] != 'M']

        # COMPARE
        overlap = set(Rdb['mut_key']).intersection(set(RNdb['mut_key']))
        Rdb = Rdb[Rdb['mut_key'].isin(overlap)]
        RNdb = RNdb[RNdb['mut_key'].isin(overlap)]

        #assert len(RNdb) > 300, len(RNdb)
        assert len(RNdb) > 250, len(RNdb)
        assert len(Rdb) == len(RNdb)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_colwidth', -1)

        Rdb = Rdb[['scaffold', 'position', 'mutation_type', 'mutation', 'gene', 'mut_key']].sort_values(['scaffold', 'position'])
        Rdb['gene'] = [0 if (g == None) | (g == 'None') else g for g in Rdb['gene']]
        RNdb = RNdb[['scaffold', 'position', 'mutation_type', 'mutation', 'gene', 'mut_key']].sort_values(['scaffold', 'position'])

        Rdb = Rdb.fillna(0)
        RNdb = RNdb.fillna(0)

        for i, row in Rdb.iterrows():
            Nrow = RNdb[RNdb['mut_key'] == row['mut_key']]
            for col in Rdb.columns:
                assert row[col] == Nrow[col].tolist()[0], [col, row, Nrow]

        assert compare_dfs(RNdb, Rdb, verbose=True)

        # Also compare gene clonalities
        Gdb = inStrain.GeneProfile.get_gene_info(IS)

    def test2(self):
        '''
        Make sure it produces output
        '''
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3} --store_everything -d".format(self.script, location, \
            base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(location + '/output/*')
        assert len(rawfiles) == 8, [len(rawfiles), rawfiles, location + '/output/*']

        # Internally verify
        IS = inStrain.SNVprofile.SNVprofile(location)
        Gdb = IS.get('genes_table')
        db = IS.get('genes_coverage')
        gene2sequence = IS.get('gene2sequence')
        RGdb = pd.merge(Gdb, db, on='gene', how='left')

        for scaff, d in RGdb.groupby('gene'):
            d = d.sort_values('mm')
            for thing in ['coverage', 'breadth']:
                assert d[thing].tolist() == sorted(d[thing].tolist()), [d, thing]
            for i, row in d.iterrows():
                seq = gene2sequence[row['gene']]
                assert len(seq) == (row['end'] - row['start'] + 1), [len(seq), row['end'] - row['start'] + 1]

        # Make sure the values make a little bit of sense when compared to the scaffolds
        diffs = []
        Sdb = IS.get('cumulative_scaffold_table')
        for scaff, sdb in Sdb.groupby('scaffold'):
            for mm, db in sdb.groupby('mm'):
                ScaffCov = db['coverage'].tolist()[0]
                GeneCov = RGdb[(RGdb['mm'] == mm) & (RGdb['scaffold'] == scaff)]['coverage'].mean()
                if GeneCov != GeneCov:
                    continue
                #print(ScaffCov, GeneCov, scaff, mm)
                #assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov-GeneCov)/ScaffCov))
        assert np.mean(diffs) < 0.15 # The average difference is less than 15%

        # Make sure logs are good
        # Make sure the missing scaffold is reported
        rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
        got = 0
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'calculate_gene_metrics' in line:
                    got += 1
        assert got == 1, got

    def test3(self):
        '''
        Make sure it works with a genbank file
        '''
        location = os.path.join(self.test_dir, os.path.basename(self.IS2))
        shutil.copytree(self.IS2, location)

        # Oh god this is horrible; re-name some variables
        tis = inStrain.SNVprofile.SNVprofile(location)
        db = tis.get('cumulative_snv_table')
        db = db.rename(columns={'varBase':'var_base', 'conBase':'con_base'})
        tis.store('cumulative_snv_table', db, 'pandas',
                    'Cumulative SNP on mm level. Formerly snpLocations.pickle')

        # db = tis.get('raw_snp_table')
        # db = db.rename(columns={'var_base':'varBase', 'con_base':'conBase'})
        # tis.store('raw_snv_table', db, 'pandas',
        #             'Contains raw SNP information on a mm level')

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3} --store_everything".format(self.script, location, \
            base, self.genbank)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(location + '/output/*')
        assert len(rawfiles) == 7, [rawfiles, location + '/output/*', len(rawfiles)]

        # Internally verify
        IS = inStrain.SNVprofile.SNVprofile(location)
        Gdb = IS.get('genes_table')
        db = IS.get('genes_coverage')
        gene2sequence = IS.get('gene2sequence')
        RGdb = pd.merge(Gdb, db, on='gene', how='left')

        for i, row in Gdb.iterrows():
            seq = gene2sequence[row['gene']]
            try:
                assert len(seq) == (row['end'] - row['start'] + 1), [len(seq), row['end'] - row['start'] + 1]
            except:
                print(row['gene'] + ' failed')
                print(len(seq))
                print(row['end'] - row['start'] + 1)
                print(row)

        for scaff, d in RGdb.groupby('gene'):
            d = d.sort_values('mm')
            for thing in ['coverage', 'breadth']:
                assert d[thing].tolist() == sorted(d[thing].tolist()), [d, thing]

        # Make sure the values make a little bit of sense when compared to the scaffolds
        diffs = []
        Sdb = IS.get('cumulative_scaffold_table')
        for scaff, sdb in Sdb.groupby('scaffold'):
            for mm, db in sdb.groupby('mm'):
                ScaffCov = db['coverage'].tolist()[0]
                GeneCov = RGdb[(RGdb['mm'] == mm) & (RGdb['scaffold'] == scaff)]['coverage'].mean()
                if GeneCov != GeneCov:
                    continue
                #print(ScaffCov, GeneCov, scaff, mm)
                #assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov-GeneCov)/ScaffCov))
        assert np.mean(diffs) < 0.16, np.mean(diffs) # The average difference is less than 16%

    def test4(self):
        '''
        Make sure it works when there are no SNVs
        '''
        location = os.path.join(self.test_dir, os.path.basename(self.IS2))
        shutil.copytree(self.IS2, location)

        # Delete SNVs
        s_loc = location + '/output/sars_cov_2_MT039887.1.fasta.bt2-vs-SRR11140750.sam.IS_SNVs.tsv'
        with open(s_loc, 'w'):
            pass
        inStrain.SNVprofile.SNVprofile(location).store('cumulative_snv_table', pd.DataFrame(),
                        'pandas', 'Cumulative SNP on mm level. Formerly snpLocations.pickle')

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3}".format(self.script, location, \
            base, self.genbank)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(location + '/output/*')
        assert len(rawfiles) == 7, [rawfiles, location + '/output/*', len(rawfiles)]

        # Internally verify
        IS = inStrain.SNVprofile.SNVprofile(location)
        Gdb = IS.get('genes_table')
        db = IS.get('genes_coverage')
        RGdb = pd.merge(Gdb, db, on='gene', how='left')

        for scaff, d in RGdb.groupby('gene'):
            d = d.sort_values('mm')
            for thing in ['coverage', 'breadth']:
                assert d[thing].tolist() == sorted(d[thing].tolist()), [d, thing]

        # Make sure the values make a little bit of sense when compared to the scaffolds
        diffs = []
        Sdb = IS.get('cumulative_scaffold_table')
        for scaff, sdb in Sdb.groupby('scaffold'):
            for mm, db in sdb.groupby('mm'):
                ScaffCov = db['coverage'].tolist()[0]
                GeneCov = RGdb[(RGdb['mm'] == mm) & (RGdb['scaffold'] == scaff)]['coverage'].mean()
                if GeneCov != GeneCov:
                    continue
                #print(ScaffCov, GeneCov, scaff, mm)
                #assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov-GeneCov)/ScaffCov))
        assert np.mean(diffs) < 0.16, np.mean(diffs) # The average difference is less than 16%
        rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
        os.remove(rr)

        # # Make sure it can make plots
        # importlib.reload(logging)
        # cmd = "inStrain plot -i {0}".format(location)
        # print(cmd)
        # call(cmd, shell=True)
        #
        # # Run the IS plots
        # FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
        #         'readANI_distribution.pdf',
        #         'ScaffoldInspection_plot.pdf',
        #         'ReadFiltering_plot.pdf',
        #         'GeneHistogram_plot.pdf']
        #
        # # Load output
        # IS = inStrain.SNVprofile.SNVprofile(location)
        # figs = glob.glob(IS.get_location('figures') + '*')
        #
        # # Make sure logging works
        # log_log = IS.get_location('log') + 'log.log'
        # rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
        # got = 0
        # with open(rr, 'r') as o:
        #     for line in o.readlines():
        #         line = line.strip()
        #         if 'Plot' in line:
        #             got += 1
        #         print(line)
        # assert got == 8, got # Its in there twice for random reasons
        #
        # # Make sure all figures are made
        # for F in FIGS:
        #     assert len([f for f in figs if F in f]) == 1, F
        # for fig in figs:
        #     assert os.path.getsize(fig) > 1000, fig

    def test5(self):
        pass

class test_filter_reads():
    def setUp(self):
        self.script = get_script_loc('filter_reads')
        self.test_dir = load_random_test_dir()
        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.readloc = load_data_loc() + \
            'N5_271_010G1.R1.fastq.gz'
        self.readloc2 = load_data_loc() + \
            'N5_271_010G1.R1.reads'

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

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

    def test0(self):
        '''
        Test basic functionality
        '''
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain filter_reads {0} {1} -o {2}".format(self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1)
        assert len(Rdb) == 179, len(Rdb)

        # Make it make the detailed read report
        cmd = "inStrain filter_reads {0} {1} -o {2} --detailed_mapping_info".format(self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        files = glob.glob(base + '/*')
        print(files)
        assert len(files) == 2
        Rdb = pd.read_csv([f for f in files if 'detailed' not in f][0], sep='\t', header=1)
        assert len(Rdb) == 179, len(Rdb)
        RRdb = pd.read_csv([f for f in files if 'detailed' in f][0], sep='\t', header=1)
        assert len(RRdb) == 41257, len(RRdb)


    def test1(self):
        '''
        Compare version 1 and 2 of getting paired reads
        '''
        # Load
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta")) # set up .fasta file
        scaffolds = list(scaff2sequence.keys())

        # Run version 2 to get scaffold to pairs to info
        s2pair2info2 = inStrain.filter_reads.get_paired_reads_multi(self.sorted_bam, scaffolds)
        db2 = inStrain.filter_reads.make_detailed_mapping_info(s2pair2info2)

        # Run version 1 to get scaffold to pairs to info
        s2pair2info, scaff2total = inStrain.deprecated_filter_reads.get_paired_reads_multi(self.sorted_bam, scaffolds, ret_total=True)
        db = inStrain.filter_reads.make_detailed_mapping_info(s2pair2info, version=1)

        # Make sure they get the same base reads
        assert len(db) == len(db2[db2['reads'] == 2])

        # Load priority_reads method 1
        reads = inStrain.filter_reads.load_priority_reads(self.readloc)
        s1 = set(db['read_pair'].tolist())
        assert len(s1.intersection(reads)) > 0

        # Load priority_reads method 2
        reads2 = inStrain.filter_reads.load_priority_reads(self.readloc2)
        assert len(s1.intersection(reads2)) > 0

        assert reads == reads2

        # Filter version two with only pairs
        kwargs = {}
        pair2info2 = inStrain.filter_reads.paired_read_filter(s2pair2info2, **kwargs)
        filtered_2 = set(p for s, p2i in pair2info2.items() for p, i in p2i.items())

        # Get the reads in a weird way and make sure it works
        filtered_1A = set()
        for s, p2i in s2pair2info.items():
            for p, i in p2i.items():
                filtered_1A.add(p)
        filtered_1 = set([p for s, p2i in s2pair2info.items() for p, i in p2i.items()])
        assert filtered_1A == filtered_1

        # Make sure the filtered version 2 is the same as regular version 1
        assert filtered_1 == filtered_2, [len(filtered_1), len(filtered_2)]

        # Make sure that not removing pairs changes this
        kwargs = {'pairing_filter':'non_discordant'}
        pair2info2N = inStrain.filter_reads.paired_read_filter(s2pair2info2, **kwargs)
        #filtered_2N = set(p for p, i in pair2info2N.items())
        filtered_2N = set([p for s, p2i in pair2info2N.items() for p, i in p2i.items()])
        assert filtered_1 != filtered_2N, [len(filtered_1), len(filtered_2N)]

        kwargs = {'pairing_filter':'all_reads'}
        pair2info2N = inStrain.filter_reads.paired_read_filter(s2pair2info2, **kwargs)
        #filtered_2A = set(p for p, i in pair2info2N.items())
        filtered_2A = set([p for s, p2i in pair2info2N.items() for p, i in p2i.items()])
        assert filtered_2N != filtered_2A, [len(filtered_2N), len(filtered_2A)]
        assert filtered_1 != filtered_2A, [len(filtered_1), len(filtered_2A)]

        # Make sure that not priority reads changes this
        priority_reads_loc = self.readloc
        priority_reads = inStrain.filter_reads.load_priority_reads(priority_reads_loc)
        kwargs = {}
        pair2info2P = inStrain.filter_reads.paired_read_filter(s2pair2info2, priority_reads_set=priority_reads, **kwargs)
        filtered_2P = set(p for p, i in pair2info2P.items())
        assert filtered_2 != filtered_2P, [len(filtered_2), len(filtered_2P)]

        # Make the old-style read report
        # Merge into single
        pair2info = {}
        for scaff, p2i in s2pair2info.items():
            for p, i in p2i.items():
                pair2info[p] = i
        RRo = inStrain.deprecated_filter_reads.makeFilterReport(s2pair2info, scaff2total=scaff2total)

        # Test out the read filtering report
        tallys = {}
        s2pair2info2_filtered = inStrain.filter_reads.paired_read_filter(s2pair2info2,
                                                                    tallys=tallys)
        j, RR = inStrain.filter_reads.filter_scaff2pair2info(s2pair2info2_filtered,
            tallys=tallys, **kwargs)

        # This now relies on the tally; that's ok, there's other tests for this
        for col in ['singletons']:
            assert RR['unfiltered_' + col].tolist()[0] > 0
            assert RR['filtered_' + col].tolist()[0] == 0
        assert RR['unfiltered_priority_reads'].tolist()[0] == 0
        assert RR['filtered_priority_reads'].tolist()[0] == 0

        # Make sure they have the right number of things in total
        assert RRo['unfiltered_pairs'].tolist()[0] == RR['unfiltered_pairs'].tolist()[0] == len(filtered_2) == len(filtered_1),\
        [RRo['unfiltered_pairs'].tolist()[0], RR['unfiltered_pairs'].tolist()[0], len(filtered_2), len(filtered_1)]

        o2n = {'pass_min_read_ani':'pass_filter_cutoff'}
        for item in ['pass_min_read_ani', 'pass_min_mapq',
            'pass_max_insert', 'pass_min_insert', 'filtered_pairs']:
            t = RR[item].tolist()[0]

            if item in o2n:
                item = o2n[item]

            o = RRo[item].tolist()[0]

            assert o == t, [item, o, t]

        # Make sure they're the same for a number of random scaffolds
        scaffolds = RR['scaffold'].tolist()[4:8]
        for scaff in scaffolds:
            for item in ['unfiltered_pairs', 'pass_min_read_ani', 'pass_min_mapq',
                'pass_max_insert', 'pass_min_insert', 'filtered_pairs']:
                t = RR[RR['scaffold'] == scaff][item].tolist()[0]

                if item in o2n:
                    item = o2n[item]
                o = RRo[RRo['scaffold'] == scaff][item].tolist()[0]

                assert o == t, [item, o, t]

        # Make sure no singletons when they're filtered out
        for col in ['singletons']:
            assert RR['unfiltered_' + col].tolist()[0] > 0
            assert RR['filtered_' + col].tolist()[0] == 0
        assert RR['unfiltered_priority_reads'].tolist()[0] == 0
        assert RR['filtered_priority_reads'].tolist()[0] == 0

        # Try out allowing singletons
        RRs = inStrain.filter_reads.makeFilterReport2(s2pair2info2, scaffold_level_mapping_info=True, **kwargs)
        for col in ['singletons']:
            assert RRs['unfiltered_' + col].tolist()[0] > 0
        assert RRs['unfiltered_priority_reads'].tolist()[0] == 0
        assert RRs['filtered_priority_reads'].tolist()[0] == 0

        # Try out priority_reads
        RRs = inStrain.filter_reads.makeFilterReport2(s2pair2info2, scaffold_level_mapping_info=True, pairTOinfo=pair2info2, priority_reads_set=priority_reads, **kwargs)
        for col in ['singletons']:
            assert RRs['unfiltered_' + col].tolist()[0] > 0
        assert RRs['unfiltered_priority_reads'].tolist()[0] > 0


    def test2(self):
        '''
        Compare method 1 and 2 of getting filtered reads
        '''
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta")) # set up .fasta file
        scaffolds = list(scaff2sequence.keys())

        # Try new method
        pair2infoF, RR = inStrain.filter_reads.load_paired_reads(self.sorted_bam, scaffolds)
        assert len(set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()]))\
                == int(RR['filtered_pairs'].tolist()[0])

        # Try old method

        # Load paired reads
        scaff2pair2info, scaff2total = inStrain.deprecated_filter_reads.get_paired_reads_multi(self.sorted_bam, scaffolds, ret_total=True)
        # Merge into single
        pair2info = {}
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                pair2info[p] = i
        # Make read report
        logging.info('Making read report')
        RRo = inStrain.deprecated_filter_reads.makeFilterReport(scaff2pair2info, scaff2total, pair2info=pair2info)
        # Filter the dictionary
        logging.info('Filtering reads')
        pair2infoFo = inStrain.deprecated_filter_reads.filter_paired_reads_dict(pair2info)
        assert len(pair2infoFo.keys()) == int(RRo['filtered_pairs'].tolist()[0])

        # Compare
        assert set(pair2infoFo.keys()) == \
                set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()])
        for s, p2i in pair2infoF.items():
            for pair, info in p2i.items():
                assert info == pair2infoFo[pair], pair

        # Try new method with priority_reads
        kwargs = {"priority_reads":self.readloc2}
        pair2infoF, RR = inStrain.filter_reads.load_paired_reads(self.sorted_bam, scaffolds, **kwargs)

        pair2infoF_pairs = set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()])

        PReads = inStrain.filter_reads.load_priority_reads(self.readloc2)
        assert set(pair2infoFo.keys()) != pair2infoF_pairs
        assert len(pair2infoF_pairs - set(pair2infoFo.keys())) > 0
        assert len(pair2infoF_pairs - set(pair2infoFo.keys()) - set(PReads)) == 0,\
            len(set(pair2infoF.keys()) - set(pair2infoFo.keys()) - set(PReads))

    def test3(self):
        '''
        Test read filtering options
        '''
        # Run on default
        base = self.test_dir + 'test'
        cmd = "inStrain filter_reads {0} {1} -o {2}".format(self.sorted_bam, \
            self.fasta, base)
        call(cmd, shell=True)

        # Load base output
        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1).head(1)
        UFB = Rdb['unfiltered_reads'].tolist()[0]
        FB = Rdb['filtered_pairs'].tolist()[0]
        assert Rdb['filtered_singletons'].tolist()[0] == 0
        assert Rdb['filtered_priority_reads'].tolist()[0] == 0

        # Run it keeping all reads
        base = self.test_dir + 'test'
        cmd = "inStrain filter_reads {0} {1} -o {2} --pairing_filter all_reads".format(self.sorted_bam, \
            self.fasta, base)
        call(cmd, shell=True)

        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1).head(1)
        UFN = Rdb['unfiltered_reads'].tolist()[0]
        FN = Rdb['filtered_pairs'].tolist()[0]
        assert Rdb['filtered_singletons'].tolist()[0] > 0
        assert Rdb['filtered_priority_reads'].tolist()[0] == 0
        assert UFN == UFB
        assert FN > FB

        # Run it keeping non_discordant reads
        base = self.test_dir + 'test'
        cmd = "inStrain filter_reads {0} {1} -o {2} --pairing_filter non_discordant".format(self.sorted_bam, \
            self.fasta, base)
        call(cmd, shell=True)

        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1).head(1)
        UFN = Rdb['unfiltered_reads'].tolist()[0]
        FN = Rdb['filtered_pairs'].tolist()[0]
        assert Rdb['filtered_singletons'].tolist()[0] > 0
        assert Rdb['filtered_priority_reads'].tolist()[0] == 0
        assert UFN == UFB
        assert FN > FB

        # Run it with priority_reads
        base = self.test_dir + 'test'
        cmd = "inStrain filter_reads {0} {1} -o {2} --priority_reads {3}".format(self.sorted_bam, \
            self.fasta, base, self.readloc2)
        call(cmd, shell=True)

        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1).head(1)
        UFN = Rdb['unfiltered_reads'].tolist()[0]
        FN = Rdb['filtered_pairs'].tolist()[0]
        assert Rdb['filtered_singletons'].tolist()[0] > 0
        assert Rdb['filtered_priority_reads'].tolist()[0] > 0
        assert UFN == UFB
        assert FN > FB

    def test4(self):
        '''
        Test that chaning the read filtering options actually changes the results
        '''
        # Run base
        base = self.test_dir + 'testMatt'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(self.sorted_bam, \
            self.fasta, base)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()

        # Run with priority_reads
        base = self.test_dir + 'testMatt'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation --priority_reads {3}".format(self.sorted_bam, \
            self.fasta, base, self.readloc2)
        call(cmd, shell=True)

        # Load the object
        Priority_object = inStrain.SNVprofile.SNVprofile(base)
        PCLdb = Priority_object.get_nonredundant_scaffold_table()

        # Compare
        assert set(PCLdb['scaffold'].tolist()) == set(MCLdb['scaffold'].tolist())
        scaffs = set(PCLdb['scaffold'].tolist())

        cl_diffs = 0
        cov_diffs = 0
        for scaffold in scaffs:
            db1 = PCLdb[PCLdb['scaffold'] == scaffold]
            db2 = MCLdb[MCLdb['scaffold'] == scaffold]

            cov1 = db1['coverage'].tolist()[0]
            cov2 = db2['coverage'].tolist()[0]
            if cov1 != cov2:
                cov_diffs += 1

            cl1 = db1['nucl_diversity'].tolist()[0]
            cl2 = db2['nucl_diversity'].tolist()[0]
            if cl1 != cl2:
                cl_diffs += 1

        print("{0} scaffolds, {1} clonality diffs, {2} coverage diffs".format(
            len(scaffs), cl_diffs, cov_diffs))
        # Make sure some, but not all, things are different
        assert cl_diffs > 0
        assert cl_diffs < len(scaffs)
        assert cov_diffs > 0
        assert cov_diffs < len(scaffs)


class test_readcomparer():
    def setUp(self):
        self.script = get_script_loc('readcomparer')
        self.test_dir = load_random_test_dir()

        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.bam1 = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.bam2 = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam'

        # These get updated by UPDATE_COMPARE_TEST_DATA
        self.IS1 = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.forRC.IS'
        self.IS2 = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.forRC.IS'

        # This gets updated by test0; for use by other tests
        self.RC_Loc = load_data_loc() + \
            'readComparer_vCurrent.RC'

        self.v12_solution = load_data_loc() + \
            'readComparer_v1.3.0r.RC'
        self.SIS = load_data_loc() + \
            'Ecoli_ani.100.0.subset.sorted.bam.IS'
        self.SIS2 = load_data_loc() + \
            'Ecoli_ani.99.9.subset.sorted.bam.IS'
        self.SIS3 = load_data_loc() + \
            'Ecoli_ani.98.0.subset.sorted.bam.IS'

        self.scafflist = load_data_loc() + \
            'scaffList.txt'
        self.scafflistF = load_data_loc() + \
            'scaffList.fasta'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.testS()
        self.tearDown()

        # self.setUp()
        # self.UPDATE_COMPARE_TEST_DATA()
        # self.tearDown()

        self.setUp()
        self.test0()
        self.tearDown()

        self.setUp()
        self.test1()
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

        self.setUp()
        self.test5()
        self.tearDown()

        self.setUp()
        self.test6()
        self.tearDown()

        self.setUp()
        self.test7()
        self.tearDown()

        self.setUp()
        self.test8()
        self.tearDown()
        #
        # # # THE .BAM FILES TO MAKE THESE IS FILES ARE DELETED FROM BIOTITE; SHOULD BE RE-GENERATED
        # # # self.setUp()
        # # # self.test9()
        # # # self.tearDown()
        # #
        # # # ### THIS IS GREEDY CLUSTERING! NOT WORKING NOW!
        # # # ### self.setUp()
        # # # ### self.test10()
        # # # ### self.tearDown()
        
        self.setUp()
        self.test11()
        self.tearDown()

        self.setUp()
        self.test12()
        self.tearDown()

        self.setUp()
        self.test13()
        self.tearDown()

    def UPDATE_COMPARE_TEST_DATA(self):
        '''
        Run inStrain on bam1 and bam2, and store the results where IS1 and IS2 are
        '''
        # Run the first one
        base1 = self.test_dir + os.path.basename(self.IS1)
        cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(self.bam1, self.fasta,
                base1)
        print(cmd)
        code = call(cmd, shell=True)
        assert code == 0, code

        # Copy to new location
        if os.path.isdir(self.IS1):
            shutil.rmtree(self.IS1)
        shutil.copytree(base1, self.IS1)

        # Run the second one
        base2 = self.test_dir + os.path.basename(self.IS2)
        cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(self.bam2, self.fasta,
                base2)
        print(cmd)
        code = call(cmd, shell=True)
        assert code == 0, code

        # Copy to new location
        if os.path.isdir(self.IS2):
            shutil.rmtree(self.IS2)
        shutil.copytree(base2, self.IS2)

    def testS(self):
        '''
        Like test0 but quicker
        '''
        # Run program
        base = self.test_dir + 'RC_test'

        cmd = "inStrain compare -i {1} {2} -o {3} -s {4} --include_self_comparisons --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflistF)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        outfiles = glob.glob(base + '/output/*')
        assert len(outfiles) == 1

        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Read the scaffold table
        Rdb = pd.read_csv(glob.glob(base + '/raw_data/*' + 'comparisonsTable.csv.gz')[0])
        # Rdb['name1'] = [x.split('-vs-')[1] for x in Rdb['name1']]
        # Rdb['name2'] = [x.split('-vs-')[1] for x in Rdb['name2']]

        # Read the pickle file
        RIS = inStrain.SNVprofile.SNVprofile(base)
        MSdb = RIS.get('pairwise_SNP_locations')

        # Make sure it's correct
        S1 = inStrain.SNVprofile.SNVprofile(self.IS1)
        S2 = inStrain.SNVprofile.SNVprofile(self.IS2)
        SRdb1 = S1.get_nonredundant_scaffold_table()
        SRdb2 = S2.get_nonredundant_scaffold_table()
        t1 = os.path.basename(S1.get('bam_loc'))
        t2 = os.path.basename(S2.get('bam_loc'))

        for scaff in S2.get('cumulative_scaffold_table')['scaffold'].unique():
            db = Rdb[Rdb['scaffold'] == scaff]
            if len(db) == 0:
                continue

            S1db = SRdb1[SRdb1['scaffold'] == scaff]
            S2db = SRdb2[SRdb2['scaffold'] == scaff]


            # Make sure self comparisons are equal to unmasked breadh
            o = db[(db['name1'] == t1) & (db['name2'] == t1)].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            tt1 = S1db['breadth_minCov'].tolist()[0]
            assert o - tt1 < 0.000001, [o, tt1]

            o = db[(db['name1'] == t2) & (db['name2'] == t2)].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            tt2 = S2db['breadth_minCov'].tolist()[0]
            assert o - tt2 < 0.0000001, [o, tt2]

            assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist()) == set([1])

            # Make sure the mms work; breadth should only go up
            for loc in [t1, t2]:
                dd = db[(db['name1'] == loc) & (db['name2'] == loc)].sort_values('mm')
                for thing in ['percent_genome_compared', 'compared_bases_count']:
                    assert dd[thing].tolist() == sorted(dd[thing].tolist()), [dd[thing].tolist(), sorted(dd[thing].tolist())]

            # Make sure the coverage_overlap is within bounds
            dd = db.sort_values('mm').drop_duplicates(
            subset=['scaffold', 'name1', 'name2'], keep='last')\
            .sort_values('scaffold')
            d = dd[dd['name1'] != dd['name2']]
            assert len(d) == 1, d
            co = d['coverage_overlap'].tolist()[0]
            assert co > (1 - (1-tt1) - (1-tt2)), "scaffold {4}; co is {0}, breadth 1 is {1}, breadh 2 is {2}, calculation is {3}".format(
                co, tt1, tt2, (1 - (1-tt1) - (1-tt2)), scaff)

            # DO SOME BASIC TESTS FOR ANI
            # When you have no ANI, percent_genome_compared must be 0
            x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
            if len(x) > 0:
                assert (x == set([0]) | x == set([]) | x == set([float(0)])), x

            # When you compare yourself, ANI must be 1
            x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
            if len(x) > 0:
                assert (x == set([1]) | x == set([]) | x == set([float(1)])), x

            # Check on the pickle;
            for i, row in db.iterrows():
                msdb = MSdb[(MSdb['name1'] == row['name1']) & (MSdb['name2'] == row['name2'])]
                snp_locs = msdb[msdb['mm'] == 'mm']['position'].tolist()

                cov = row['compared_bases_count']
                ani = row['popANI']
                snps = cov - (float(ani) * cov)
                try:
                    snps = int(round(snps))
                except:
                    snps = 0
                assert len(snp_locs) == snps, [len(snp_locs), snps, cov, ani, msdb]

    def test0(self):
        '''
        Test calling readcomparer. These are two reads from samples from the same DOL to the same assembly.

        All of these checks are on the internals; that is comparing self to self, basically
        '''
        # Run program
        base = self.test_dir + 'RC_test'

        cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        outfiles = glob.glob(base + '/output/*')
        assert len(outfiles) == 1

        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5


        # Read the scaffold table
        Rdb = pd.read_csv(glob.glob(base + '/raw_data/*' + 'comparisonsTable.csv.gz')[0])
        # Rdb['name1'] = [x.split('-vs-')[1] for x in Rdb['name1']]
        # Rdb['name2'] = [x.split('-vs-')[1] for x in Rdb['name2']]

        # Read the pickle file
        RIS = inStrain.SNVprofile.SNVprofile(base)
        MSdb = RIS.get('pairwise_SNP_locations')

        # Make sure it's correct
        S1 = inStrain.SNVprofile.SNVprofile(self.IS1)
        S2 = inStrain.SNVprofile.SNVprofile(self.IS2)
        SRdb1 = S1.get_nonredundant_scaffold_table()
        SRdb2 = S2.get_nonredundant_scaffold_table()
        t1 = os.path.basename(S1.get('bam_loc'))
        t2 = os.path.basename(S2.get('bam_loc'))

        for scaff in S1.get('cumulative_scaffold_table')['scaffold'].unique():
            db = Rdb[Rdb['scaffold'] == scaff]
            if len(db) == 0:
                print("skipping " + scaff)
                continue

            if set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist()) != set([1]):
                assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist()) == set(), set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist())

            S1db = SRdb1[SRdb1['scaffold'] == scaff]
            S2db = SRdb2[SRdb2['scaffold'] == scaff]

            # Make sure self comparisons are equal to unmasked breadh
            o = db[(db['name1'] == t1) & (db['name2'] == t1)].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            tt1 = S1db['breadth_minCov'].tolist()[0]
            assert o - tt1 < 0.000001, [o, tt1]

            o = db[(db['name1'] == t2) & (db['name2'] == t2)].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            tt2 = S2db['breadth_minCov'].tolist()[0]
            assert o - tt2 < 0.0000001, [o, tt2]

            # Make sure the mms work; breadth should only go up
            for loc in [t1, t2]:
                dd = db[(db['name1'] == loc) & (db['name2'] == loc)].sort_values('mm')
                for thing in ['percent_genome_compared', 'compared_bases_count']:
                    assert dd[thing].tolist() == sorted(dd[thing].tolist()), [dd[thing].tolist(), sorted(dd[thing].tolist())]

            # Make sure the coverage_overlap is within bounds
            dd = db.sort_values('mm').drop_duplicates(
            subset=['scaffold', 'name1', 'name2'], keep='last')\
            .sort_values('scaffold')
            d = dd[dd['name1'] != dd['name2']]
            assert len(d) == 1
            co = d['coverage_overlap'].tolist()[0]
            assert co > (1 - (1-tt1) - (1-tt2)), "scaffold {4}; co is {0}, breadth 1 is {1}, breadh 2 is {2}, calculation is {3}".format(
                co, tt1, tt2, (1 - (1-tt1) - (1-tt2)), scaff)

            # DO SOME BASIC TESTS FOR ANI
            # When you have no ANI, percent_genome_compared must be 0
            x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
            if len(x) > 0:
                assert (x == set([0]) | x == set([]) | x == set([float(0)])), x

            # When you compare yourself, ANI must be 1
            x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
            if len(x) > 0:
                assert (x == set([1]) | x == set([]) | x == set([float(1)])), x

            # Check on the pickle;
            for i, row in db.iterrows():
                msdb = MSdb[(MSdb['name1'] == row['name1']) & (MSdb['name2'] == row['name2']) \
                        & (MSdb['mm'] == float(row['mm'])) & (MSdb['scaffold'] == row['scaffold'])]
                snp_locs = msdb[(msdb['population_SNP'] == True)]['position'].tolist()

                cov = row['compared_bases_count']
                ani = row['popANI']
                snps = cov - (float(ani) * cov)
                try:
                    snps = int(round(snps))
                except:
                    snps = 0

                assert snps == int(row['population_SNPs']), [snps, row]
                assert len(snp_locs) == snps, [len(snp_locs), snps, cov, ani, msdb]

        # THIS UPDATES THE RC_LOC!
        if os.path.isdir(self.RC_Loc):
            shutil.rmtree(self.RC_Loc)
        shutil.copytree(base, self.RC_Loc)

    def test1(self):
        '''
        Test readcomparer functions
        '''
        # Load
        S1 = inStrain.SNVprofile.SNVprofile(self.IS1)
        S2 = inStrain.SNVprofile.SNVprofile(self.IS2)

        # Get the scaffolds to compare
        scaffs = set(S1.get('cumulative_scaffold_table')['scaffold'].unique())\
                    .intersection(set(S2.get('cumulative_scaffold_table')['scaffold'].unique()))

        # Get the covTs
        covTs = [S1.get('covT'), S2.get('covT')]

        # Get the SNP tables
        SNP_tables = [S1.get('cumulative_snv_table'), S2.get('cumulative_snv_table')]

        # Get the nonredundant scaffold table
        SRdb1 = S1.get_nonredundant_scaffold_table()
        SRdb2 = S2.get_nonredundant_scaffold_table()

        s2l = S1.get('scaffold2length')

        # Run it
        for scaff in scaffs:
            # if scaff != 'N5_271_010G1_scaffold_6':
            #     continue
            db, pair2mm2SNPlocs, pair2mm2cov, scaffold = inStrain.readComparer.compare_scaffold(scaff, ['t1', 't2'], [S1, S2], s2l[scaff], include_self_comparisons=True)

            # MAKE SURE THE COVERAGES ARE RIGHT
            S1db = SRdb1[SRdb1['scaffold'] == scaff]
            S2db = SRdb2[SRdb2['scaffold'] == scaff]

            o = db[(db['name1'] == 't1') & (db['name2'] == 't1')].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            t = S1db['breadth_minCov'].tolist()[0]
            assert o - t < 0.0001, [o, t]

            o = db[(db['name1'] == 't2') & (db['name2'] == 't2')].sort_values('mm')\
                    ['percent_genome_compared'].tolist()[-1]
            t = S2db['breadth_minCov'].tolist()[0]
            assert o - t < 0.0001, [o, t]

            if set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist()) != set([1]):
                assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist()) == set(), set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist())

            # DO SOME BASIC TESTS FOR ANI
            x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
            if len(x) > 0:
                assert (x == set([0]) | x == set([]) | x == set([float(0)])), x

            x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
            if len(x) > 0:
                assert (x == set([1]) | x == set([]) | x == set([float(1)])), x

    def test2(self):
        '''
        Test the min_coverage function
        '''
        # Run program normally
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons".format(self.script, self.IS1, self.IS2, \
            base)
        print(cmd)
        call(cmd, shell=True)

        # Run program with high min_cov
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3}.2 -c 50 --include_self_comparisons".format(self.script, self.IS1, self.IS2, \
            base)
        print(cmd)
        call(cmd, shell=True)

        # Read the scaffold tables
        Rdb = inStrain.SNVprofile.SNVprofile(base).get_nonredundant_RC_table()
        Rdb2 = inStrain.SNVprofile.SNVprofile(base + '.2').get_nonredundant_RC_table()
        # Rdb = pd.read_csv(base + '.comparisonsTable.csv.gz').sort_values('mm').drop_duplicates(
        #     subset=['scaffold', 'name1', 'name2'], keep='last')\
        #     .sort_values('scaffold')
        # Rdb2 = pd.read_csv(base + '.2.comparisonsTable.csv.gz').sort_values('mm').drop_duplicates(
        #     subset=['scaffold', 'name1', 'name2'], keep='last')\
        #     .sort_values('scaffold')

        # Compare
        assert set(Rdb['scaffold'].tolist()) == set(Rdb2['scaffold'].tolist())

        higher = 0
        total = 0
        for scaff in set(Rdb['scaffold'].tolist()):
            rdb = Rdb[(Rdb['scaffold'] == scaff) & (Rdb['name1'] == Rdb['name2'])]
            rdb2 = Rdb2[(Rdb2['scaffold'] == scaff) & (Rdb2['name1'] == Rdb2['name2'])]
            assert len(rdb) == len(rdb2) == 2

            for name in rdb['name1'].tolist():
                b1 = rdb[rdb['name1'] == name]['percent_genome_compared'].tolist()[0]
                b2 = rdb2[rdb2['name1'] == name]['percent_genome_compared'].tolist()[0]
                if b1 > b2:
                    higher += 1
                total += 1

        assert higher/total >= 0.75, "{0} of {1} comparisons are higher".format(higher, total)

    def test3(self):
        '''
        Test random things
        '''
        P2C = {'A':0, 'C':1, 'T':2, 'G':3, 'X':4}
        S = inStrain.SNVprofile.SNVprofile(self.IS1)
        db = S.get('cumulative_snv_table')

        db.at['ref_base', 0] = np.nan
        db.at['con_base', 0] = np.nan
        db.at['ref_base', 1] = 'W'
        db.at['con_base', 1] = 'Y'

        db['con_base'] = [x if x in P2C else 'X' for x in db['con_base']]
        db['ref_base'] = [x if x in P2C else 'X' for x in db['con_base']]

        db['con_base'] = db['con_base'].map(P2C).astype(int)
        db['ref_base'] = db['ref_base'].map(P2C).astype(int)

    def test4(self):
        '''
        Test providing a .fasta file of scaffolds
        '''
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} -s {4}".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        Rdb = inStrain.SNVprofile.SNVprofile(base).get('comparisonsTable')
        #Rdb = pd.read_csv(base + '.comparisonsTable.csv.gz')
        scaffs = set(inStrain.controller.load_scaff_list(self.scafflist))
        assert set(scaffs) == set(Rdb['scaffold'].tolist())

    def test5(self):
        '''
        Test store_coverage_overlap
        '''
        # Run program
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} -s {4} --store_coverage_overlap".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

    def test6(self):
        '''
        Test --store_mismatch_locations
        '''
        # Run program
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} -s {4} --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

    def test7(self):
        '''
        Test --compare_consensus_bases
        '''
        # Run program
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} -s {4} --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Load snps
        scaff2snps = {}
        MSdb = inStrain.SNVprofile.SNVprofile(base).get('comparisonsTable')
        for i, row in MSdb.iterrows():
            if row['conANI'] != row['conANI']:
                continue
            assert row['conANI'] <= row['popANI'], row


    def test8(self):
        '''
        Test readComparer fdr adjustment
        '''
        # Run program
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Load snps
        scaff2popsnps = {}
        scaff2consnps = {}
        MSdb = inStrain.SNVprofile.SNVprofile(base).get('pairwise_SNP_locations')
        for scaff, msdb in MSdb.groupby('scaffold'):
            if scaff not in scaff2popsnps:
                scaff2popsnps[scaff] = set()
                scaff2consnps[scaff] = set()

            consnps = set(msdb[msdb['consensus_SNP'] == True]['position'].tolist())
            popsnps = set(msdb[msdb['population_SNP'] == True]['position'].tolist())
            scaff2consnps[scaff] = scaff2consnps[scaff].union(consnps)
            scaff2popsnps[scaff] = scaff2popsnps[scaff].union(popsnps)

        # Print the total number
        total_con_snps = sum([len(s) for sc, s in scaff2consnps.items()])
        total_pop_snps = sum([len(s) for sc, s in scaff2popsnps.items()])


        # Run program with compare_consensus_bases
        base = self.test_dir + 'testR2'

        cmd = "inStrain compare -i {1} {2} -o {3} --fdr 0.5 --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Load snps
        scaff2popsnps = {}
        scaff2consnps = {}
        MSdb = inStrain.SNVprofile.SNVprofile(base).get('pairwise_SNP_locations')
        for scaff, msdb in MSdb.groupby('scaffold'):
            if scaff not in scaff2popsnps:
                scaff2popsnps[scaff] = set()
                scaff2consnps[scaff] = set()

            consnps = set(msdb[msdb['consensus_SNP'] == True]['position'].tolist())
            popsnps = set(msdb[msdb['population_SNP'] == True]['position'].tolist())
            scaff2consnps[scaff] = scaff2consnps[scaff].union(consnps)
            scaff2popsnps[scaff] = scaff2popsnps[scaff].union(popsnps)

        # Print the total number
        total_con_snps2 = sum([len(s) for sc, s in scaff2consnps.items()])
        total_pop_snps2 = sum([len(s) for sc, s in scaff2popsnps.items()])

        # Make sure you catch less SNPs when you lower your fdr
        assert total_con_snps2 == total_con_snps, [total_con_snps2, total_con_snps]
        assert total_pop_snps2 < total_pop_snps, [total_pop_snps2, total_pop_snps]

    def test9(self):
        '''
        Test readComparer on synthetic dataset
        '''
        # Run program
        base = self.test_dir + 'testR'

        cmd = "inStrain compare -i {1} {2} -o {3} --store_mismatch_locations".format(self.script, self.SIS, self.SIS2, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Load snps
        scaff2popsnps = {}
        scaff2consnps = {}
        MSdb = inStrain.SNVprofile.SNVprofile(base).get('pairwise_SNP_locations')
        for scaff, msdb in MSdb.groupby('scaffold'):
            if scaff not in scaff2popsnps:
                scaff2popsnps[scaff] = set()
                scaff2consnps[scaff] = set()

            consnps = set(msdb[msdb['consensus_SNP'] == True]['position'].tolist())
            popsnps = set(msdb[msdb['population_SNP'] == True]['position'].tolist())
            scaff2consnps[scaff] = scaff2consnps[scaff].union(consnps)
            scaff2popsnps[scaff] = scaff2popsnps[scaff].union(popsnps)

        total_con_snps = sum([len(s) for sc, s in scaff2consnps.items()])
        total_pop_snps = sum([len(s) for sc, s in scaff2popsnps.items()])

        # Load the SNV tables
        SNdb1 = inStrain.SNVprofile.SNVprofile(self.SIS).get_nonredundant_snv_table()
        SNdb2 = inStrain.SNVprofile.SNVprofile(self.SIS2).get_nonredundant_snv_table()

        # Make sure they're all called
        assert (len(SNdb2) + len(SNdb1)) == total_con_snps == total_pop_snps == 11, \
                [len(SNdb2), len(SNdb1), total_con_snps, total_pop_snps]

        # TRY THE OTHER OPNE

        base = self.test_dir + 'testX'

        cmd = "inStrain compare -i {2} {1} -o {3} --store_mismatch_locations".format(self.script, self.SIS2, self.SIS3, \
            base, self.scafflist)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        rawfiles = glob.glob(base + '/raw_data/*')
        assert len(rawfiles) == 5

        # Load snps
        scaff2popsnps = {}
        scaff2consnps = {}
        MSdb = inStrain.SNVprofile.SNVprofile(base).get('pairwise_SNP_locations')
        for scaff, msdb in MSdb.groupby('scaffold'):
            if scaff not in scaff2popsnps:
                scaff2popsnps[scaff] = set()
                scaff2consnps[scaff] = set()

            consnps = set(msdb[msdb['consensus_SNP'] == True]['position'].tolist())
            popsnps = set(msdb[msdb['population_SNP'] == True]['position'].tolist())
            scaff2consnps[scaff] = scaff2consnps[scaff].union(consnps)
            scaff2popsnps[scaff] = scaff2popsnps[scaff].union(popsnps)
        total_con_snps = sum([len(s) for sc, s in scaff2consnps.items()])
        total_pop_snps = sum([len(s) for sc, s in scaff2popsnps.items()])
        assert total_con_snps == total_pop_snps, [total_con_snps, total_pop_snps]

        # Load the SNV table
        SNdb2 = inStrain.SNVprofile.SNVprofile(self.SIS2).get_nonredundant_snv_table()
        SNdb3 = inStrain.SNVprofile.SNVprofile(self.SIS3).get_nonredundant_snv_table()

        # Find mixed positions
        Positions = set(SNdb2['position']).union(set(SNdb3['position']))

        # Load the coverage
        Cov2 = inStrain.SNVprofile.SNVprofile(self.SIS2).get('covT')['CP011321.1']
        Cov2 = inStrain.profileUtilities._mm_counts_to_counts_shrunk(Cov2)

        Cov3 = inStrain.SNVprofile.SNVprofile(self.SIS3).get('covT')['CP011321.1']
        Cov3 = inStrain.profileUtilities._mm_counts_to_counts_shrunk(Cov3)

        to_rm = []
        for p in sorted(Positions):
            #print("{0} - {1} {2}".format(p, Cov2.loc[p], Cov3.loc[p]))
            if (Cov2.loc[p] < 5) | (Cov3.loc[p] < 5):
                to_rm.append(p)
        for p in to_rm:
            Positions.remove(p)


        # Make sure they're all called
        assert len(Positions) == total_con_snps == 18, [len(Positions) , total_con_snps]

    # def test10(self):
    #     '''
    #     Test greedy clustering
    #     '''
    #     # Run program
    #     base = self.test_dir + 'testR'
    #
    #     # Make a copy of one
    #     location = os.path.join(self.test_dir, os.path.basename(self.IS.replace('N5_271_010G1', 'testeroni')))
    #     shutil.copytree(self.IS, location)
    #     inStrain.SNVprofile.SNVprofile(location).store('bam_loc', 'testeroni', 'value', 'Location of .bam file')
    #
    #     cmd = "inStrain compare -i {1} {2} {5} -o {3} -s {4} --greedy_clustering --g_cov 0.46".format(self.script, self.IS, self.IS2, \
    #         base, self.scafflist, location)
    #     print(cmd)
    #     call(cmd, shell=True)
    #
    #     # Make sure it produced output
    #     outfiles = glob.glob(base + '/output/*')
    #     assert len(outfiles) == 2
    #
    #     # Load the clusters file
    #     Cdb = pd.read_csv(base + '/output/testR_greedyClusters.tsv', sep='\t')
    #     print(Cdb)
    #     assert len(Cdb['name'].unique()) == 3
    #     assert len(Cdb['cluster'].unique()) == 2

    def test11(self):
        '''
        Test skipping scaffolds
        '''
        S1 = inStrain.SNVprofile.SNVprofile(self.IS1)
        S2 = inStrain.SNVprofile.SNVprofile(self.IS2)

        s2l = S1.get('scaffold2length')
        Sprofiles = [S1, S2]
        names = ['name1', 'name2']
        scaffolds_to_compare = list(S1._get_covt_keys())[:3]
        scaffolds_to_compare.append('junkero')
        s2l['junkero'] = 1000

        # Run program
        Cdb, Mdb, scaff2pair2mm2cov = inStrain.readComparer.compare_scaffolds(names, Sprofiles, scaffolds_to_compare, s2l)

        # NO REAL WAY TO TEST THIS; JUST MAKE SURE THAT THE OUTPUT DOESNT SHOW RE-RUNNING ANY SCAFFOLDS

    def test12(self):
        '''
        Test _calc_SNP_count_alternate
        '''
        # Make table1
        table = defaultdict(list)

        # Make a consensus SNP in table1
        table['position'].append(1)
        table['con_base'].append('A')
        table['ref_base'].append('C')
        table['var_base'].append('A')
        table['position_coverage'].append(100)
        table['A'].append(100)
        table['C'].append(0)
        table['G'].append(0)
        table['T'].append(0)
        table['allele_count'].append(1)
        table['mm'].append(0)

        # Make a population SNP in table1
        table['position'].append(2)
        table['con_base'].append('A')
        table['ref_base'].append('C')
        table['var_base'].append('C')
        table['position_coverage'].append(100)
        table['A'].append(60)
        table['C'].append(40)
        table['G'].append(0)
        table['T'].append(0)
        table['allele_count'].append(2)
        table['mm'].append(0)

        SNPtable1 = pd.DataFrame(table)
        SNPtable2 = pd.DataFrame()

        mm2overlap = {0:[1,2]}

        mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                mm2overlap)
        assert len(mdb[mdb['consensus_SNP'] == True]) == 2
        assert len(mdb[mdb['population_SNP'] == True]) == 1

        # Change up table2
        table = defaultdict(list)

        # Make a consensus SNP in table1
        table['position'].append(1)
        table['con_base'].append('T')
        table['ref_base'].append('C')
        table['var_base'].append('T')
        table['position_coverage'].append(100)
        table['A'].append(0)
        table['C'].append(0)
        table['G'].append(0)
        table['T'].append(100)
        table['allele_count'].append(1)
        table['mm'].append(0)

        # Make a population SNP in table1
        table['position'].append(2)
        table['con_base'].append('T')
        table['ref_base'].append('C')
        table['var_base'].append('G')
        table['position_coverage'].append(100)
        table['A'].append(0)
        table['C'].append(0)
        table['G'].append(40)
        table['T'].append(60)
        table['allele_count'].append(2)
        table['mm'].append(0)

        SNPtable2 = pd.DataFrame(table)
        mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                mm2overlap)
        assert len(mdb[mdb['consensus_SNP'] == True]) == 2
        assert len(mdb[mdb['population_SNP'] == True]) == 2

        # Change up table2 again
        table = defaultdict(list)

        table['position'].append(1)
        table['con_base'].append('A')
        table['ref_base'].append('C')
        table['var_base'].append('A')
        table['position_coverage'].append(100)
        table['A'].append(80)
        table['C'].append(20)
        table['G'].append(0)
        table['T'].append(0)
        table['allele_count'].append(2)
        table['mm'].append(0)

        # Make a population SNP in table1
        table['position'].append(2)
        table['con_base'].append('G')
        table['ref_base'].append('C')
        table['var_base'].append('C')
        table['position_coverage'].append(100)
        table['A'].append(0)
        table['C'].append(40)
        table['G'].append(60)
        table['T'].append(0)
        table['allele_count'].append(2)
        table['mm'].append(0)

        SNPtable2 = pd.DataFrame(table)
        mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                mm2overlap)
        assert len(mdb[mdb['consensus_SNP'] == True]) == 1
        assert len(mdb[mdb['population_SNP'] == True]) == 0


        # Try with nothing
        SNPtable1 = pd.DataFrame()
        SNPtable2 = pd.DataFrame()
        mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                mm2overlap)
        assert len(mdb) == 0

    def test13(self):
        '''
        Re-run and ensure that the results are the same as a previous run
        '''
        # Run program
        base = self.test_dir + 'RC_test'

        cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons --store_mismatch_locations".format(self.script, self.IS1, self.IS2, \
            base, self.scafflistF)
        print(cmd)
        call(cmd, shell=True)

        exp_RC = inStrain.SNVprofile.SNVprofile(base)
        sol_RC = inStrain.SNVprofile.SNVprofile(self.v12_solution)

        # Print what the output of the solutions directory looks like
        if True:
            s_out_files = glob.glob(exp_RC.get_location('output') + os.path.basename(
                                    exp_RC.get('location')) + '_*')
            print("The output has {0} tables".format(len(s_out_files)))
            for f in s_out_files:
                name = os.path.basename(f)
                print("{1}\n{0}\n{1}".format(name, '-'*len(name)))

                if 'mapping_info.tsv' in name:
                    s = pd.read_csv(f, sep='\t', header=1)
                else:
                    s = pd.read_csv(f, sep='\t')
                print(s.head())
                print()

        # Make sure log is working
        assert len(glob.glob(base + '/log/*')) == 2, glob.glob(base + '/log/*')
        Ldb = exp_RC.get_parsed_log()
        print(Ldb)

        # Check output files
        e_out_files = glob.glob(exp_RC.get_location('output') + os.path.basename(
                                exp_RC.get('location')) + '_*')
        s_out_files = glob.glob(sol_RC.get_location('output') + '*_*')
        assert len(s_out_files) == 1, sol_RC.get_location('output') + '*_*'

        for s_file in s_out_files:
            name = os.path.basename(s_file).split('RC_test_')[1]
            e_file = [e for e in e_out_files if name in os.path.basename(e)]

            print("checking {0}".format(name))

            if len(e_file) == 1:
                #print("Both have {0}!".format(name))

                e = pd.read_csv(e_file[0], sep='\t')
                s = pd.read_csv(s_file, sep='\t')

                if name == 'comparisonsTable.tsv':
                    e = e.sort_values(['scaffold', 'name1', 'name2']
                                        ).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'name1', 'name2']
                                        ).reset_index(drop=True)

                assert set(s.columns) == set(e.columns), \
                        [set(s.columns) - set(e.columns),
                         set(e.columns) - set(s.columns),]
                s = s[list(e.columns)]
                assert compare_dfs2(e, s, verbose=True), name

            else:
                assert False, name

        # Check attributes
        sAdb = sol_RC._get_attributes_file()
        eAdb = exp_RC._get_attributes_file()

        for i, row in sAdb.iterrows():
            print("checking {0}".format(i))

            if i in ['location', 'version']:
                continue

            s = sol_RC.get(i)
            e = exp_RC.get(i)

            if i in ['comparisonsTable', 'pairwise_SNP_locations']:
                s = s.sort_values(['scaffold', 'name1', 'name2']).reset_index(drop=True)
                e = e.sort_values(['scaffold', 'name1', 'name2']).reset_index(drop=True)

                # Re-arange column order
                assert set(e.columns) == set(s.columns),\
                    [i,
                     set(e.columns) - set(s.columns),\
                     set(s.columns) - set(e.columns)]
                s = s[list(e.columns)]
                assert compare_dfs2(e, s, verbose=True), i

            elif i in ['scaffold2length']:
                assert compare_dicts(e, s, verbose=True), i


class test_strains():
    def setUp(self, destroy=True):
        self.script = get_script_loc('inStrain')

        self.test_dir = load_random_test_dir()
        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'
        self.genes = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa.genes.fna'

        self.failure_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_failureScaffold.sorted.bam'
        self.single_scaff = load_data_loc() + \
            'N5_271_010G1_scaffold_101.fasta'
        self.fasta_extra = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000_extra.fa'
        self.small_fasta = load_data_loc() + \
            'SmallScaffold.fa'
        self.small_bam = load_data_loc() + \
            'SmallScaffold.fa.sorted.bam'
        self.extra_single_scaff = load_data_loc() + \
            'N5_271_010G1_scaffold_101_extra.fasta'
        self.failure_fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_failureScaffold.fa'
        self.failure_genes = load_data_loc() + \
            'N5_271_010G1_scaffold_failureScaffold.fa.genes.fna.fa'
        self.cc_solution = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam.CB'
        self.pp_snp_solution = load_data_loc() + \
            'strainProfiler_v0.3_results/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam_SP_snpLocations.pickle'
        self.cc_snp_solution = load_data_loc() + \
            'v0.4_results/test_0.98.freq'
        self.v12_solution = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.2.14'
        self.sam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam'
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
        self.scafflist = load_data_loc() + \
            'scaffList.txt'
        self.genes = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.stb = load_data_loc() + \
            'GenomeCoverages.stb'

        if destroy:
            if os.path.isdir(self.test_dir):
                shutil.rmtree(self.test_dir)
            os.mkdir(self.test_dir)

        importlib.reload(logging)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        # YOU HAVE TO RUN THIS ONE ON ITS OWN, BECUASE IT MESSES UP FUTURE RUNS
        # self.setUp()
        # self.test0()
        # self.tearDown()

        self.setUp()
        self.test1()
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

        self.setUp()
        self.test5()
        self.tearDown()

        self.setUp()
        self.test6()
        self.tearDown()

        self.setUp()
        self.test7()
        self.tearDown()

        self.setUp()
        self.test8()
        self.tearDown()

        self.setUp()
        self.test9()
        self.tearDown()

        self.setUp()
        self.test10()
        self.tearDown()

        self.setUp()
        self.test11()
        self.tearDown()

        self.setUp()
        self.test12()
        self.tearDown()

        self.setUp()
        self.test13()
        self.tearDown()

        self.setUp()
        self.test14()
        self.tearDown()

        self.setUp()
        self.test15()
        self.tearDown()

        self.setUp(destroy=True)
        self.test16()
        self.tearDown()

        self.setUp()
        self.test17()
        self.tearDown()

        self.setUp()
        self.test18()
        self.tearDown()

    def test0(self):
        '''
        Compare Matts to CCs methodology
        '''
        # Run Matts program
        base = self.test_dir + 'testMatt'
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --min_mapq 2 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        sys_args = cmd.split(' ')
        args = inStrain.argumentParser.parse_args(sys_args[1:])
        inStrain.controller.Controller().main(args)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)

        # Run CCs program
        base = self.test_dir + 'testCC'
        cmd = "{0} {1} {2} -o {3} -l 0.95".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        args = inStrain.deprecated.parse_arguments(cmd.split()[1:])
        inStrain.deprecated.main(args)

        # Load the object
        CC_object = inStrain.deprecated.SNVdata()
        CC_object.load(name=base + '_0.95')

        ## COMPARE SNPS

        # Parse CCs dumb SNV table
        CPdb = CC_object.snv_table
        CPdb['var_base'] = [s.split(':')[1] for s in CPdb['SNV']]
        CPdb['loc'] = [int(s.split(':')[0].split('_')[-1]) for s in CPdb['SNV']]
        CPdb['scaff'] = ['_'.join(s.split(':')[0].split('_')[:-1]) for s in CPdb['SNV']]
        CPdb = CPdb.drop_duplicates(subset=['scaff','loc'])
        CPdb = CPdb.sort_values(['scaff', 'loc'])

        # Load Matts beautiful object
        MPdb = Matt_object.get_nonredundant_snv_table()

        # Allowing for cryptic SNPs, make sure Matt calls everything CC does
        MS = set(["{0}-{1}".format(x, y) for x, y in zip(MPdb['scaffold'], MPdb['position'])])
        CS = set(["{0}-{1}".format(x, y) for x, y in zip(CPdb['scaff'], CPdb['loc'])])
        assert len(MS) > 0
        assert len(CS - MS) == 0, CS-MS

        # Not allowing for cyptic SNPs, make sure CC calls everything Matt does
        MPdb = MPdb[MPdb['cryptic'] == False]
        MPdb = MPdb[MPdb['allele_count'] >= 2]
        MS = set(["{0}-{1}".format(x, y) for x, y in zip(MPdb['scaffold'], MPdb['position'])])
        CS = set(["{0}-{1}".format(x, y) for x, y in zip(CPdb['scaff'], CPdb['loc'])])
        assert len(MS - CS) == 0, MS-CS

        ## COMPARE CLONALITY

        # Parse CCs dumb clonality table
        CLdb = CC_object.clonality_table
        p2c = CLdb.set_index('position')['clonality'].to_dict()

        # Load Matt's beautiful table
        MCLdb = Matt_object.get_clonality_table()

        #print(set(p2c.keys()) - set(["{0}_{1}".format(s, p) for s, p in zip(MCLdb['scaffold'], MCLdb['position'])]))
        assert len(MCLdb) == len(CLdb), (len(MCLdb), len(CLdb))
        for i, row in MCLdb.dropna().iterrows():
            assert (p2c["{0}_{1}".format(row['scaffold'], row['position'])] \
                - row['clonality']) < .001, (row, p2c["{0}_{1}".format(row['scaffold'], row['position'])])

        ## COMPARE LINKAGE
        CLdb = CC_object.r2linkage_table
        CLdb['position_A'] = [eval(str(x))[0].split(':')[0].split('_')[-1] for x in CLdb['total_A']]
        CLdb['position_B'] = [eval(str(x))[0].split(':')[0].split('_')[-1] for x in CLdb['total_B']]
        CLdb['scaffold'] = [x.split(':')[0] for x  in CLdb['Window']]

        MLdb = Matt_object.get_nonredundant_linkage_table()
        # Mark cryptic SNPs
        MPdb = Matt_object.get_nonredundant_snv_table()
        dbs = []
        for scaff, db in MPdb.groupby('scaffold'):
            p2c = db.set_index('position')['cryptic'].to_dict()
            mdb = MLdb[MLdb['scaffold'] == scaff]
            mdb['cryptic'] = [True if ((p2c[a] == True) | (p2c[b] == True)) else False for a, b in zip(
                                mdb['position_A'], mdb['position_B'])]
            dbs.append(mdb)
        MLdb = pd.concat(dbs)

        # # Make sure MLdb and MPdb aggree
        # MLS = set(["{0}-{1}".format(x, y, z) for x, y, z in zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])]).union(
        #       set(["{0}-{2}".format(x, y, z) for x, y, z in zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])]))
        # print([len(MS), len(MLS), len(MS - MLS), len(MLS - MS)])
        # assert MS == MLS

        # Allowing for cryptic SNPs, make sure Matt calls everything CC does
        MS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])])
        CS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in zip(CLdb['scaffold'], CLdb['position_A'], CLdb['position_B'])])
        assert len(CS - MS) <= 1, [CS - MS]
        # At scaffold N5_271_010G1_scaffold_110 from position 525 to 546 you end up in an edge case
        # where you skip it because you have absolutely no minor alleles to counterbalance it. It's fine,
        # CC just reports an r2 of np.nan, and this seems like the easiest way to handle it

        # Not allowing for cyptic SNPs, make sure CC calls everything Matt does
        MLdb = MLdb[MLdb['cryptic'] == False]
        MS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])])
        CS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in zip(CLdb['scaffold'], CLdb['position_A'], CLdb['position_B'])])
        assert len(MS - CS) == 0, [len(MS), len(CS), len(MS - CS), MS - CS]

    def test1(self):
        '''
        Basic test- Make sure whole version doesn't crash when run from the command line
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98".format(self.script, self.sorted_bam, \
            self.fasta_extra, base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert os.path.isdir(base)
        assert len(glob.glob(base + '/output/*')) == 5, glob.glob(base + '/output/*')

        # Make sure the output makes sense
        S1 = inStrain.SNVprofile.SNVprofile(base)
        db = S1.get('cumulative_scaffold_table')
        _internal_verify_Sdb(db)

        # Make sure it doesn't mess up at the lower-case bases (I put a lower-case in scaffold N5_271_010G1_scaffold_0 at a c; make sure it's not there)
        db = S1.get('raw_snp_table')
        assert len(db[db['ref_base'] == 'c']) == 0

    def test2(self):
        '''
        Test filter reads; make sure CCs and Matt's agree
        '''
        # Set up
        positions, total_length = inStrain.deprecated_filter_reads.get_fasta(self.fasta)
        min_read_ani = 0.98

        # Run initial filter_reads
        subset_reads, Rdb = inStrain.deprecated_filter_reads.filter_reads(self.sorted_bam, positions, total_length,
                            filter_cutoff=min_read_ani, max_insert_relative=3, min_insert=50, min_mapq=2)

        # Run Matts filter_reads
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta")) # set up .fasta file
        s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length
        scaffolds = list(s2l.keys())
        subset_reads2 = inStrain.deprecated_filter_reads.filter_paired_reads(self.sorted_bam,
                        scaffolds, filter_cutoff=min_read_ani, max_insert_relative=3,
                        min_insert=50, min_mapq=2)

        # Run Matts filter_reads in a different way
        pair2info = inStrain.deprecated_filter_reads.get_paired_reads(self.sorted_bam, scaffolds)
        pair2infoF = inStrain.deprecated_filter_reads.filter_paired_reads_dict(pair2info,
                        filter_cutoff=min_read_ani, max_insert_relative=3,
                        min_insert=50, min_mapq=2)
        subset_readsF = list(pair2infoF.keys())

        # Run Matts filter_reads in a different way still
        scaff2pair2infoM, Rdb = inStrain.filter_reads.load_paired_reads(
                            self.sorted_bam, scaffolds, min_read_ani=min_read_ani,
                            max_insert_relative=3, min_insert=50, min_mapq=2)
        # pair2infoMF = inStrain.filter_reads.paired_read_filter(scaff2pair2infoM)
        # pair2infoMF = inStrain.filter_reads.filter_scaff2pair2info(pair2infoMF,
        #                 min_read_ani=min_read_ani, max_insert_relative=3,
        #                 min_insert=50, min_mapq=2)

        subset_readsMF = set()
        for scaff, pair2infoC in scaff2pair2infoM.items():
            subset_readsMF = subset_readsMF.union(pair2infoC.keys())
        # subset_readsMF = list(pair2infoMF.keys())

        assert (set(subset_reads2) == set(subset_reads) == set(subset_readsF) == set(subset_readsMF)),\
                [len(subset_reads2), len(subset_reads), len(subset_readsF), len(subset_readsMF)]

        # Make sure the filter report is accurate
        # Rdb = inStrain.filter_reads.makeFilterReport2(scaff2pair2infoM, pairTOinfo=pair2infoMF, min_read_ani=min_read_ani, max_insert_relative=3,
        #                                     min_insert=50, min_mapq=2)
        assert int(Rdb[Rdb['scaffold'] == 'all_scaffolds']\
                    ['unfiltered_pairs'].tolist()[0]) \
                    == len(list(pair2info.keys()))
        assert int(Rdb[Rdb['scaffold'] == 'all_scaffolds']['filtered_pairs'].tolist()[0]) \
                == len(subset_reads)

        # Try another cutuff
        positions, total_length = inStrain.deprecated_filter_reads.get_fasta(self.fasta)
        min_read_ani = 0.90

        # Run initial filter_reads
        subset_reads, Rdb = inStrain.deprecated_filter_reads.filter_reads(self.sorted_bam, positions, total_length,
                            filter_cutoff=min_read_ani, max_insert_relative=3, min_insert=50, min_mapq=2)

        # Run Matts filter_reads
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta")) # set up .fasta file
        s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length
        scaffolds = list(s2l.keys())

        scaff2pair2infoM, Rdb = inStrain.filter_reads.load_paired_reads(
                            self.sorted_bam, scaffolds, min_read_ani=min_read_ani,
                            max_insert_relative=3, min_insert=50, min_mapq=2)
        pair2infoMF_keys = set()
        for scaff, pair2infoC in scaff2pair2infoM.items():
            pair2infoMF_keys = pair2infoMF_keys.union(pair2infoC.keys())
        # Scaff2pair2infoM = inStrain.filter_reads.get_paired_reads_multi(self.sorted_bam, scaffolds)
        # pair2infoMF = inStrain.filter_reads.paired_read_filter(scaff2pair2infoM)
        # pair2infoMF = inStrain.filter_reads.filter_paired_reads_dict2(pair2infoMF,
        #                 min_read_ani=min_read_ani, max_insert_relative=3,
        #                 min_insert=50, min_mapq=2)


        subset_reads2 = pair2infoMF_keys

        assert(set(subset_reads2) == set(subset_reads))

    def test3(self):
        '''
        Testing scaffold table (breadth and coverage) vs. calculate_breadth
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)
        Odb = Sprofile.get('cumulative_scaffold_table')

        # Verify coverage table
        _internal_verify_Sdb(Odb)

        # Compare to calculate_coverage
        Cdb = pd.read_csv(self.cc_solution)
        s2c = Cdb.set_index('scaffold')['coverage'].to_dict()
        s2b = Cdb.set_index('scaffold')['breadth'].to_dict()
        for scaff, db in Odb.groupby('scaffold'):
            db = db.sort_values('mm', ascending=False)
            assert (db['coverage'].tolist()[0] - s2c[scaff]) < .1,\
                            [db['coverage'].tolist()[0], s2c[scaff]]
            assert (db['breadth'].tolist()[0] - s2b[scaff]) < .01,\
                            [db['breadth'].tolist()[0], s2b[scaff]]

        # Verify SNP calls
        Sdb = Sprofile.get('cumulative_snv_table')
        _internal_verify_OdbSdb(Odb, Sdb)

    def test4(self):
        '''
        Test store_everything and database mode
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --skip_plot_generation -s {4}".format(self.script, self.sorted_bam, \
            self.fasta, base, self.stb)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you stored a heavy object
        assert Sprofile.get('testeroni') is  None
        assert Sprofile.get('covT') is not None
        assert Sprofile.get('mm_to_position_graph') is not None

        # Run database mode
        base2 = self.test_dir + 'test2'
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --database_mode --skip_plot_generation -s {4}".format(self.script, self.sorted_bam, \
            self.fasta, base2, self.stb)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile2 = inStrain.SNVprofile.SNVprofile(base2)

        # Make sure you didn't story a heavy object
        assert Sprofile2.get('mm_to_position_graph') is None

        # Make sure you have more reads mapping
        mdb1 = Sprofile.get('mapping_info').sort_values('scaffold').reset_index(
                drop=True)
        mdb2 = Sprofile2.get('mapping_info').sort_values('scaffold').reset_index(
                drop=True)
        assert set(mdb1['scaffold']) == set(mdb2['scaffold'])
        assert not compare_dfs2(mdb1, mdb2, verbose=True)

        # Make sure you have more skip mm level
        mdb1 = Sprofile.get('genome_level_info')
        mdb2 = Sprofile2.get('genome_level_info')
        assert 'mm' in mdb1
        assert 'mm' not in mdb2

        # Make sure you skip junk genomes
        # assert len(set(mdb1['genome']) - set(mdb2['genome'])) > 0

    def test5(self):
        '''
        Test one thread
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 1 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get a result
        Odb = Sprofile.get('cumulative_scaffold_table')
        assert len(Odb['scaffold'].unique()) == 178


    def test6(self):
        '''
        Test the case where only one scaffold is  preset at all in the .bam file AND it has no SNPs
        '''
        # Run program
        base = self.test_dir + 'test'
        cmd = "inStrain profile {1} {2} -o {3} -l 0.99 -p 1 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.single_scaff, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Load output
        Odb = Sprofile.get('cumulative_scaffold_table')
        print(Odb)
        _internal_verify_Sdb(Odb)

    def test7(self):
        '''
        Test the case where a scaffold is not preset at all in the .bam file

        Also test being able to adjust read filtering parameters
        '''
        # Run program
        base = self.test_dir + 'test'
        cmd = "inStrain profile {1} {2} -o {3} -l 0.80 -p 6 --store_everything --skip_genome_wide --skip_plot_generation".format(self.script, self.sorted_bam, \
            self.extra_single_scaff, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure scaffold table is OK
        Odb = Sprofile.get('cumulative_scaffold_table')
        _internal_verify_Sdb(Odb)

        # Check the read report
        rloc = glob.glob(Sprofile.get_location('output') + '*mapping_info.tsv')[0]
        with open(rloc) as f:
            first_line = f.readline()
        assert "min_read_ani:0.8" in first_line

        rdb = pd.read_csv(rloc, sep='\t', header=1)
        total_pairs = rdb[rdb['scaffold'] == 'all_scaffolds']['filtered_pairs'].tolist()[0]
        reads = set()
        for s, rdic in Sprofile.get('Rdic').items():
            reads = reads.union(rdic.keys())
        assert total_pairs == len(reads)

        ORI_READS = len(reads)

        for thing, val in zip(['min_mapq', 'max_insert_relative', 'min_insert'], [10, 1, 100]):
            print("!!!!!")
            print(thing, val)

            cmd = "inStrain profile {1} {2} -o {3} --{4} {5} -p 6 --store_everything --skip_plot_generation".format(self.script, self.sorted_bam, \
                self.extra_single_scaff, base, thing, val)
            inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
            Sprofile = inStrain.SNVprofile.SNVprofile(base)

            rloc = glob.glob(Sprofile.get_location('output') + 'test_mapping_info.tsv')[0]
            with open(rloc) as f:
                first_line = f.readline()
            assert "{0}:{1}".format(thing, val) in first_line, [first_line, thing, val]

            if thing == 'max_insert_relative':
                thing = 'max_insert'

            rdb = pd.read_csv(rloc, sep='\t', header=1)
            passF = rdb[rdb['scaffold'] == 'all_scaffolds']["pass_{0}".format(thing)].tolist()[0]
            print(passF)
            #assert rdb[rdb['scaffold'] == 'all_scaffolds']["pass_{0}".format(thing)].tolist()[0] == 0

            reads = len(Sprofile.get('Rdic').keys())
            print(Sprofile.get('Rdic'))
            assert reads < ORI_READS

    def test8(self):
        '''
        Test the ability to make and sort .bam files from .sam
        '''
        # Copy sam to test dir
        new_sam = os.path.join(self.test_dir, os.path.basename(self.sam))
        shutil.copyfile(self.sam, new_sam)

        # Run program
        base = self.test_dir + 'test'
        cmd = "inStrain profile {1} {2} -o {3} -l 0.80 -p 6 --store_everything --skip_genome_wide -d".format(self.script, new_sam, \
            self.extra_single_scaff, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Load output
        assert len([f for f in glob.glob(base + '/output/*') if '.log' not in f]) == 4, glob.glob(base + '/output/*')

        # Make sure the missing scaffold is reported
        print(glob.glob(base + '/log/*'))
        rr = [f for f in glob.glob(base + '/log/*') if 'runtime' in f][0]
        got = False
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'weird_NAMED_scaff' in line:
                    got = True
        assert got

    def test9(self):
        '''
        Test the debug option

        v0.5.1 - Actually this should happen all the time now...
        v1.2.0 - This test is obsolete now
        '''
        pass
        # # Run Matts program
        # base = self.test_dir + 'testMatt'
        # # cmd = "{0} {1} {2} -o {3} -l 0.95 --store_everything --debug".format(self.script, self.sorted_bam, \
        # #     self.fasta, base)
        # cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --skip_genome_wide".format(self.script, self.sorted_bam, \
        #     self.fasta, base)
        # print(cmd)
        # inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        # Sprofile = inStrain.SNVprofile.SNVprofile(base)
        #
        # # Open the log
        # logfile = Sprofile.get_location('log') + 'log.log'
        #
        # table = defaultdict(list)
        # with open(logfile) as o:
        #     for line in o.readlines():
        #         line = line.strip()
        #         if 'RAM. System has' in line:
        #             linewords = [x.strip() for x in line.split()]
        #             table['scaffold'].append(linewords[0])
        #             table['PID'].append(linewords[2])
        #             table['status'].append(linewords[3])
        #             table['time'].append(linewords[5])
        #             table['process_RAM'].append(linewords[7])
        #             table['system_RAM'].append(linewords[11])
        #             table['total_RAM'].append(linewords[13])
        #             logged = True
        # Ldb = pd.DataFrame(table)
        # assert len(Ldb) > 5

    def test10(self):
        '''
        Test min number of reads filtered and min number of genome coverage
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_scaffold_reads 10 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you actually filtered out the scaffolds
        Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
        Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
        print("{0} of {1} scaffolds have >10 reads".format(len(Rdb[Rdb['filtered_pairs'] >= 10]),
                    len(Rdb)))
        assert len(Sdb['scaffold'].unique()) == len(Rdb[Rdb['filtered_pairs'] >= 10]['scaffold'].unique()) - 1

        # Try min_genome_coverage

        # Make sure it fails with no .stb
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 5 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        exit_code = call(cmd, shell=True)
        assert exit_code == 1

        # Run it with an .stb
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 5  -s {4} --skip_genome_wide".format(
            self.script, self.sorted_bam, self.fasta, base, self.stb)
        exit_code = call(cmd, shell=True)
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you actually filtered out the scaffolds
        Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
        Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
        assert len(Rdb) == 179
        assert len(Sdb) == 42

        # Make sure empty scaffolds don't mess it up
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 5  -s {4} --skip_genome_wide".format(
            self.script, self.sorted_bam, self.fasta_extra, base, self.stb)
        exit_code = call(cmd, shell=True)
        Sprofile = inStrain.SNVprofile.SNVprofile(base)
        assert exit_code == 0; exit_code

        # Make sure you actually filtered out the scaffolds
        Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
        Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
        assert len(Rdb) == 180, len(Rdb)
        assert len(Sdb) == 42, len(Sdb)

        # Run it with an .stb and coverage that cannot be hit
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 100  -s {4} --skip_genome_wide".format(
            self.script, self.sorted_bam, self.fasta, base, self.stb)
        exit_code = call(cmd, shell=True)
        Sprofile = inStrain.SNVprofile.SNVprofile(base)
        assert exit_code == 1

        # Run it with an .stb and coverage that is low
        cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 1.1  -s {4} --skip_genome_wide".format(
            self.script, self.sorted_bam, self.fasta, base, self.stb)
        exit_code = call(cmd, shell=True)
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you actually filtered out the scaffolds
        Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
        Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
        assert len(Rdb) == 179
        assert len(Sdb) > 42


    def test11(self):
        '''
        Test skip mm profiling
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} --skip_mm_profiling --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        IS = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get the same results
        scaffdb = IS.get_nonredundant_scaffold_table().reset_index(drop=True)
        correct_scaffdb = inStrain.SNVprofile.SNVprofile(self.IS).get_nonredundant_scaffold_table().reset_index(drop=True)

        sdb = IS.get('cumulative_scaffold_table')
        scdb = inStrain.SNVprofile.SNVprofile(self.IS).get('cumulative_scaffold_table')

        cols = ['scaffold', 'length', 'breadth', 'coverage']
        assert compare_dfs(scaffdb[cols], correct_scaffdb[cols], verbose=True)

        # Make sure you dont have the raw mm
        sdb = IS.get('cumulative_scaffold_table')
        print(sdb.head())
        assert set(sdb['mm'].tolist())  == set([0]), set(sdb['mm'].tolist())

    def test12(self):
        '''
        Test scaffolds_to_profile
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --min_mapq 2 --scaffolds_to_profile {4} --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base, self.scafflist)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        IS = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get the same results
        scaffdb = IS.get_nonredundant_scaffold_table().reset_index(drop=True)
        assert set(scaffdb['scaffold'].unique()) == set(inStrain.controller.load_scaff_list(self.scafflist))

    def test13(self):
        '''
        Make sure that covT, clonT, and the SNP table agree on coverage
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --min_mapq 2 --scaffolds_to_profile {4} -p 1 --skip_genome_wide".format(self.script, self.sorted_bam, \
            self.fasta, base, self.scafflist)
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        IS = inStrain.SNVprofile.SNVprofile(base)

        SRdb = IS.get('cumulative_snv_table')
        CovT = IS.get('covT')
        ClonT = IS.get('clonT')

        for i, row in SRdb.iterrows():
            cov = 0
            for mm, covs in CovT[row['scaffold']].items():
                if mm <= row['mm']:
                    if row['position'] in covs:
                        cov += covs[row['position']]
            assert row['position_coverage'] == cov, [cov, row['position_coverage'], row]

    def test14(self):
        '''
        Basic test- Make sure genes and genome_wide can be run within the profile option

        Make sure logging produces the run statistics
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98".format(self.script, self.sorted_bam, \
            self.fasta, base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert os.path.isdir(base)
        assert len(glob.glob(base + '/output/*')) == 6, len(glob.glob(base + '/output/*'))
        assert len(glob.glob(base + '/log/*')) == 3

        # Make sure the output makes sense
        S1 = inStrain.SNVprofile.SNVprofile(base)
        db = S1.get('cumulative_scaffold_table')
        _internal_verify_Sdb(db)

        # Read the log
        for l in glob.glob(base + '/log/*'):
            if 'runtime_summary.txt' in l:
                with open(l, 'r') as o:
                    for line in o.readlines():
                        print(line.strip())

    def test15(self):
        '''
        Basic test- Make sure genes and genome_wide can be run within the profile option
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98 --rarefied_coverage 10".format(self.script, self.sorted_bam, \
            self.fasta, base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert os.path.isdir(base)
        assert len(glob.glob(base + '/output/*')) == 6

        # Make sure the output makes sense
        S1 = inStrain.SNVprofile.SNVprofile(base)
        db = S1.get('cumulative_scaffold_table')
        _internal_verify_Sdb(db)
        clontR = S1.get('clonTR')
        counts0 = sum([len(x[2]) if 2 in x else 0 for s, x in clontR.items()])

        # Make sure its in the genome_wide table
        gdb = pd.read_csv(glob.glob(base + '/output/*genome_info*.tsv')[0], sep='\t')
        assert 'nucl_diversity' in gdb.columns, gdb.head()

        # Run again with different rarefied coverage
        base = self.test_dir + 'test2'
        cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98 --rarefied_coverage 50".format(self.script, self.sorted_bam, \
            self.fasta, base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        S1 = inStrain.SNVprofile.SNVprofile(base)
        db = S1.get('cumulative_scaffold_table')
        _internal_verify_Sdb(db)
        clontR = S1.get('clonTR')
        counts2 = sum([len(x[2]) if 2 in x else 0 for s, x in clontR.items()])

        assert counts0 > counts2, [counts0, counts2]

    def test16(self):
        '''
        Make sure the results exactly match a run done with inStrain verion 1.2.14
        '''
        # Run program
        base = self.test_dir + 'test'
        cmd = "inStrain profile {1} {2} -o {3} -g {4} --skip_plot_generation -p 6 -d".format(self.script, self.sorted_bam, \
            self.fasta, base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        exp_IS = inStrain.SNVprofile.SNVprofile(base)
        sol_IS = inStrain.SNVprofile.SNVprofile(self.v12_solution)

        # Print what the output of the solutions directory looks like
        if True:
            s_out_files = glob.glob(exp_IS.get_location('output') + os.path.basename(
                                    exp_IS.get('location')) + '_*')
            print("The output has {0} tables".format(len(s_out_files)))
            for f in s_out_files:
                name = os.path.basename(f)
                print("{1}\n{0}\n{1}".format(name, '-'*len(name)))

                if 'mapping_info.tsv' in name:
                    s = pd.read_csv(f, sep='\t', header=1)
                else:
                    s = pd.read_csv(f, sep='\t')
                print(s.head())
                print()


        # MAKE SURE LOG IS WORKING
        assert len(glob.glob(base + '/log/*')) == 3, base
        Ldb = exp_IS.get_parsed_log()
        rdb = inStrain.logUtils._load_profile_logtable(Ldb)

        LOGGED_SCAFFOLDS = set(rdb[rdb['command'] == 'MergeProfile']['unit'].tolist())
        TRUE_SCAFFOLDS = \
            set(exp_IS.get_nonredundant_scaffold_table()['scaffold'].tolist())
        assert(LOGGED_SCAFFOLDS == TRUE_SCAFFOLDS)

        # CHECK OUTPUT FILES
        e_out_files = glob.glob(exp_IS.get_location('output') + os.path.basename(
                                exp_IS.get('location')) + '_*')
        s_out_files = glob.glob(sol_IS.get_location('output') + os.path.basename(
                                sol_IS.get('location')) + '_*')

        for s_file in s_out_files:
            name = os.path.basename(s_file).split('v1.2.14_')[1]
            e_file = [e for e in e_out_files if name in os.path.basename(e)]

            print("checking {0}".format(name))

            if len(e_file) == 1:
                #print("Both have {0}!".format(name))

                e = pd.read_csv(e_file[0], sep='\t')
                s = pd.read_csv(s_file, sep='\t').rename(columns=twelve2thirteen)

                if name in ['linkage.tsv']:
                    e = e.sort_values(['scaffold', 'position_A', 'position_B']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'position_A', 'position_B']).reset_index(drop=True)

                    # Delete random ones
                    rand = ['r2_normalized', 'd_prime_normalized']
                    for r in rand:
                        del e[r]
                        del s[r]

                if name in ['SNVs.tsv']:
                    e = e.sort_values(['scaffold', 'position']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'position']).reset_index(drop=True)

                    for col in ['mutation_type', 'mutation', 'gene']:
                        del e[col]

                if name in ['scaffold_info.tsv']:
                    e = e.sort_values(['scaffold']).reset_index(drop=True)
                    s = s.sort_values(['scaffold']).reset_index(drop=True)

                    # TRANSLATE THE OLD VERSION
                    s = s.rename(columns=twelve2thirteen)
                    for r in del_thirteen:
                        if r in s.columns:
                            del s[r]
                    for r in new_thirteen:
                        if r in e.columns:
                            del e[r]

                    rand = ['nucl_diversity_rarefied', 'nucl_diversity_rarefied_median']
                    for r in rand:
                        if r in e.columns:
                            del e[r]
                        if r in s.columns:
                            del s[r]

                if name in ['mapping_info.tsv']:
                    e = pd.read_csv(e_file[0], sep='\t', header=1)
                    s = pd.read_csv(s_file, sep='\t', header=1)

                    e = e.sort_values(['scaffold']).reset_index(drop=True)
                    s = s.sort_values(['scaffold']).reset_index(drop=True)

                if name in ['gene_info.tsv']:
                    e = e.sort_values(['scaffold', 'gene']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'gene']).reset_index(drop=True)

                    # TRANSLATE THE OLD VERSION
                    s = s.rename(columns=twelve2thirteen)
                    for r in del_thirteen:
                        if r in s.columns:
                            del s[r]
                    for r in new_thirteen:
                        if r in e.columns:
                            del e[r]

                    rand = ['SNS_count', 'divergent_site_count', 'partial']
                    for r in rand:
                        if r in e.columns:
                            del e[r]
                        if r in s.columns:
                            del s[r]

                # Re-arange column order
                assert set(e.columns) == set(s.columns),\
                    [name,
                     set(e.columns) - set(s.columns),\
                     set(s.columns) - set(e.columns)]
                s = s[list(e.columns)]
                e = e[list(e.columns)]

                assert compare_dfs2(e, s, verbose=True), name

            else:
                #print("Both dont have {0}!".format(name))
                if name in ['read_report.tsv']:
                    e_file = [e for e in e_out_files if 'mapping_info.tsv' in os.path.basename(e)]

                    e = pd.read_csv(e_file[0], sep='\t', header=1)
                    s = pd.read_csv(s_file, sep='\t', header=1).rename(columns={'pass_filter_cutoff':'pass_min_read_ani'})

                    e = e.sort_values(['scaffold']).reset_index(drop=True)
                    s = s.sort_values(['scaffold']).reset_index(drop=True)

                    e = e[e['scaffold'] == 'all_scaffolds']
                    s = s[s['scaffold'] == 'all_scaffolds']

                    for c in list(s.columns):
                        s[c] = s[c].astype(e[c].dtype)

                    for c in ['median_insert']: # calculated in a different way
                        del e[c]
                        del s[c]

                elif name in ['genomeWide_scaffold_info.tsv']:
                    e_file = [e for e in e_out_files if 'genome_info.tsv' in os.path.basename(e)]

                    e = pd.read_csv(e_file[0], sep='\t')
                    s = pd.read_csv(s_file, sep='\t').rename(columns=twelve2thirteen)

                    for r in del_thirteen:
                        if r in s.columns:
                            del s[r]
                    s = s.rename(columns=twelve2thirteen)

                    NEW_HERE = {'coverage_median', 'SNV_count', 'SNS_count',
                                'nucl_diversity', 'filtered_read_pair_count'}
                    for c in new_thirteen.union(NEW_HERE):
                        if c in e.columns:
                            del e[c]

                    # Remove the ones that are gained by the new read filtering
                    for c in e.columns:
                        if c.startswith('reads_'):
                            del e[c]

                    e = e.sort_values(['genome']).reset_index(drop=True)
                    s = s.sort_values(['genome']).reset_index(drop=True)

                    # Re-order columns
                    assert set(e.columns) == set(s.columns),\
                        [set(e.columns) - set(s.columns),\
                         set(s.columns) - set(e.columns)]

                    e = e[list(s.columns)]
                    s = s[list(s.columns)]

                    changed_cols = ['coverage_std',
                                    'nucl_diversity_rarefied']
                    for c in changed_cols:
                        del e[c]
                        del s[c]

                elif name in ['genomeWide_read_report.tsv']:
                    e_file = [e for e in e_out_files if 'genome_info.tsv' in os.path.basename(e)]

                    e = pd.read_csv(e_file[0], sep='\t')
                    s = pd.read_csv(s_file, sep='\t').rename(columns=twelve2thirteen)\
                            .rename(columns={'pass_filter_cutoff':'pass_min_read_ani'})

                    # TRANSLATE THE OLD VERSION
                    s = s.rename(columns=twelve2thirteen)
                    for r in del_thirteen:
                        if r in s.columns:
                            del s[r]
                    for r in new_thirteen:
                        if r in e.columns:
                            del e[r]

                    new_cols = list(set(e.columns) - set(s.columns))
                    for c in new_cols:
                        del e[c]

                    removed_cols = ['unfiltered_reads', 'pass_pairing_filter', 'pass_min_mapq', 'unfiltered_singletons', 'filtered_priority_reads', 'mean_insert_distance', 'mean_pair_length', 'pass_min_insert', 'pass_max_insert', 'pass_min_read_ani', 'mean_PID', 'mean_mistmaches', 'median_insert', 'mean_mapq_score', 'unfiltered_pairs', 'filtered_singletons', 'unfiltered_priority_reads', 'filtered_pairs', 'scaffolds']
                    for r in removed_cols:
                        if r in s.columns:
                            del s[r]

                    e = e.sort_values(['genome']).reset_index(drop=True)
                    s = s.sort_values(['genome']).reset_index(drop=True)
                    assert set(e.columns) == set(s.columns), [set(s.columns) - set(e.columns)]

                elif name in ['SNP_mutation_types.tsv']:
                    e_file = [e for e in e_out_files if 'SNVs.tsv' in os.path.basename(e)]

                    e = pd.read_csv(e_file[0], sep='\t')
                    s = pd.read_csv(s_file, sep='\t').rename(columns=twelve2thirteen)

                    e = e[~e['mutation_type'].isna()]
                    del e['cryptic']

                    e = e.sort_values(['scaffold', 'position']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'position']).reset_index(drop=True)

                else:
                    assert False, name

                # Re-arange column order
                assert set(e.columns) == set(s.columns),\
                    [name,
                     set(e.columns) - set(s.columns),\
                     set(s.columns) - set(e.columns)]

                s = s[list(e.columns)]

                assert compare_dfs2(e, s, verbose=True), name

        # CHECK ATTRIBUTES
        sAdb = sol_IS._get_attributes_file()
        eAdb = exp_IS._get_attributes_file()

        # Handle name changes
        o2n = {'genes_SNP_density':'genes_SNP_count', 'read_report':'mapping_info'}

        for i, row in sAdb.iterrows():
            print("checking {0}".format(i))

            if i in ['location', 'version', 'bam_loc', 'genes_fileloc', 'window_table']:
                continue

            s = sol_IS.get(i)
            if i in o2n:    i = o2n[i]
            e = exp_IS.get(i)

            if i in ['scaffold_list']:
                if not compare_lists(e, s, verbose=True, check_order=False):
                    print("{0} is not the same".format(i))
                    print(e)
                    print(s)
                    assert e == s, i

            elif i in ['scaffold2length', 'scaffold2bin', 'bin2length']:
                assert compare_dicts(e, s, verbose=True), i

            elif i in ['window_table', 'raw_linkage_table', 'raw_snp_table', 'cumulative_scaffold_table',
                        'cumulative_snv_table', 'mapping_info', 'genes_table', 'genes_coverage',
                        'genes_clonality', 'genes_SNP_count', 'SNP_mutation_types']:

                # TRANSLATE THE OLD VERSION
                s = s.rename(columns=twelve2thirteen).rename(
                            columns={'pass_filter_cutoff':'pass_min_read_ani'})
                for r in del_thirteen:
                    if r in s.columns:
                        del s[r]
                for r in new_thirteen:
                    if r in e.columns:
                        del e[r]

                if i in ['window_table', 'mapping_info']:
                    e = e.sort_values(['scaffold']).reset_index(drop=True)
                    s = s.sort_values(['scaffold']).reset_index(drop=True)
                    e = e.rename(columns={'filtered_pairs':'filtered_read_pair_count'})

                if i in ['mapping_info']:
                    e = e[e['scaffold'] == 'all_scaffolds']
                    s = s[s['scaffold'] == 'all_scaffolds']

                    for c in list(s.columns):
                        s[c] = s[c].astype(e[c].dtype)

                    for c in ['median_insert']:
                        del e[c]
                        del s[c]

                if i in ['raw_linkage_table']:
                    continue

                    e = e.sort_values(['scaffold', 'position_A', 'position_B']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'position_A', 'position_B']).reset_index(drop=True)

                    # Delete random ones
                    rand = ['r2_normalized', 'd_prime_normalized']
                    for r in rand:
                        del e[r]
                        del s[r]

                if i in ['raw_snp_table', 'cumulative_snv_table', 'SNP_mutation_types']:
                    e = e.sort_values(['scaffold', 'position']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'position']).reset_index(drop=True)

                if i in ['genes_table']:
                    rand = ['partial']
                    for r in rand:
                        if r in e.columns:
                            del e[r]
                        if r in s.columns:
                            del s[r]

                if i in ['cumulative_scaffold_table']:
                    e = e.sort_values(['scaffold', 'mm']).reset_index(drop=True)
                    s = s.sort_values(['scaffold', 'mm']).reset_index(drop=True)

                    # Delete random ones
                    rand = ['nucl_diversity_rarefied', 'nucl_diversity_rarefied_median']
                    for r in rand:
                        if r in e.columns:
                            del e[r]
                        if r in s.columns:
                            del s[r]

                if i in ['genes_coverage', 'genes_clonality']:
                    e = e.sort_values(['gene', 'mm']).reset_index(drop=True)
                    s = s.sort_values(['gene', 'mm']).reset_index(drop=True)

                if i in ['genes_SNP_count']:
                    e = e.sort_values(['gene', 'mm']).reset_index(drop=True)
                    s = s.sort_values(['gene', 'mm']).reset_index(drop=True)

                    e = e[list(s.columns)]

                # Re-arange column order
                assert set(e.columns) == set(s.columns),\
                    [i,
                     set(e.columns) - set(s.columns),\
                     set(s.columns) - set(e.columns)]

                s = s[list(e.columns)]

                assert compare_dfs2(e, s, verbose=True), i

            elif i in ['scaffold_2_mm_2_read_2_snvs', 'clonTR']:
                pass

            elif i in ['covT', 'clonT']:
                assert compare_covTs(e, s), "{0} is not the same".format(i)

            else:
                print("YOUR NOT CHECKING {0}".format(i))
                print(s)

    def test17(self):
        '''
        Test scaffold failure
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program and make the split crash
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=3000 -d".format(self.script, self.failure_bam, \
            self.failure_fasta, base)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get a result
        Odb = Sprofile.get('cumulative_scaffold_table')
        assert len(Odb['scaffold'].unique()) == 1, Odb['scaffold'].unique()

        # Make sure the missing scaffold is reported
        rr = [f for f in glob.glob(base + '/log/*') if 'runtime' in f][0]
        got = 0
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'FailureScaffoldHeaderTesting' in line:
                    got += 1
        assert got == 3, got
        os.remove(rr)

        # Make it not crash on that scaffold
        importlib.reload(logging)
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=20000 -d".format(self.script, self.failure_bam, \
            self.failure_fasta, base)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get a result
        Odb = Sprofile.get('cumulative_scaffold_table')
        assert len(Odb['scaffold'].unique()) == 2, Odb['scaffold'].unique()

        # Make sure the missing scaffold is reported
        rr = [f for f in glob.glob(base + '/log/*') if 'runtime' in f][0]
        got = 0
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'FailureScaffoldHeaderTesting' in line:
                    got += 1
        assert got == 3, got
        os.remove(rr)

        # Make it crash on the gene profile
        importlib.reload(logging)
        cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=20000 -d -g {4}".format(self.script, self.failure_bam, \
            self.failure_fasta, base, self.failure_genes)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        Sprofile = inStrain.SNVprofile.SNVprofile(base)

        # Make sure you get a result
        Odb = Sprofile.get('cumulative_scaffold_table')
        assert len(Odb['scaffold'].unique()) == 2, Odb['scaffold'].unique()

        # Make sure the missing scaffold is reported
        rr = [f for f in glob.glob(base + '/log/*') if 'runtime' in f][0]
        got = 0
        with open(rr, 'r') as o:
            for line in o.readlines():
                line = line.strip()
                if 'FailureScaffoldHeaderTesting' in line:
                    got += 1
        assert got == 4, got
        os.remove(rr)

    def test18(self):
        '''
        Test providing a .fasta file with a really small scaffold
        '''
        base = self.test_dir + 'testR'

        cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation --debug".format(self.small_bam, self.small_fasta, \
            base)
        print(cmd)
        call(cmd, shell=True)

        IS = inStrain.SNVprofile.SNVprofile(base)

        db = IS.get('genome_level_info')
        assert db is not None
        assert len(glob.glob(base + '/log/*')) == 3


def _internal_verify_Sdb(Sdb):
    for scaff, d in Sdb.groupby('scaffold'):
        d = d.sort_values('mm')
        for thing in ['breadth_minCov', 'coverage', 'coverage_median']:
            assert d[thing].tolist() == sorted(d[thing].tolist()), d

        # clons = [1-c for c in d['mean_clonality']]
        # micros = [c for c in d['mean_microdiversity']]
        # for c, m in zip(clons, micros):
        #     if c == c:
        #         assert abs(c - m) < 0.0001, [c, m, scaff, d]

    sdb = Sdb[Sdb['breadth_minCov'] > 0]
    for col in sdb.columns:
        if col in ['nucl_diversity_rarefied', 'nucl_diversity_rarefied_median']:
            continue
        if len(sdb[sdb[col].isna()]) > 0:
            print(col)
            print(sdb[sdb[col].isna()])

    for i, row in Sdb.iterrows():
        if row['consensus_divergent_sites'] > 0:
            assert row['conANI_reference'] != 0

    assert Sdb['conANI_reference'].max() <= 1
    assert Sdb['breadth_minCov'].max() <= 1
    assert Sdb['breadth'].max() <= 1

    if 'breadth_rarefied' not in sdb.columns:
        assert len(sdb) == len(sdb.dropna())
    else:
        sdb = sdb[sdb['breadth_rarefied'] > 0]
        assert len(sdb) == len(sdb.dropna())

def _internal_verify_OdbSdb(Odb, Sdb):
    '''
    Odb = cumulative scaffold table
    Sdb = SNP table
    '''
    # Ensure internal consistancy between Sdb and Cdb at the lowest mm
    low_mm = Sdb['mm'].min()
    for scaff, db in Sdb[Sdb['mm'] == low_mm].groupby('scaffold'):
        snps = Odb['divergent_site_count'][(Odb['scaffold'] == scaff) & (Odb['mm'] \
                == low_mm)].fillna(0).tolist()[0]
        assert snps == len(db), [snps, len(db), scaff, low_mm]

    # Ensure internal consistancy between Sdb and Cdb at the highset mm
    odb = Odb.sort_values('mm').drop_duplicates(subset='scaffold', keep='last')
    for scaff, db in Sdb.sort_values('mm').drop_duplicates(subset=['scaffold'\
                    ,'position'], keep='last').groupby('scaffold'):
        snps = odb['divergent_site_count'][(odb['scaffold'] == scaff)].fillna(0).tolist()[0]
        assert snps == len(db), [snps, len(db), scaff]

def compare_dfs(db1, db2, round=4, verbose=False):
    '''
    Return True if dataframes are equal (order of dataframes doesn't matter)
    '''

    db1 = db1.fillna(0).round(round)
    db2 = db2.fillna(0).round(round)

    df = pd.concat([db1, db2], sort=True)
    df = df.reset_index(drop=True)
    df_gpby = df.groupby(list(df.columns))
    idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]

    identicle = (len(idx) == 0)
    if ((not identicle) and verbose):
        print("index: ", idx)
        print("db1: ",db1)
        print("db2: ",db2)
        print("df_gpby: ", str(df_gpby))

    return identicle

def compare_dfs2(db1, db2, round=4, verbose=False):
    '''
    Return True if dataframes are equal (order of dataframes doesn't matter)
    '''
    try:
        pd._testing.assert_frame_equal(db1, db2)
        return True
    except AssertionError as e:
        if verbose:
            print(e)
        return False

def compare_dfs3(db1, db2, round=4, verbose=False):
    '''
    Return True if dataframes are equal (order of dataframes doesn't matter)
    '''
    if len(db1) != len(db2):
        print("lengths are different")
        return False

    if set(db1.columns) != set(db2.columns):
        print("columns are different")
        return False

def compare_lists(l1, l2, verbose=False, check_order=True):
    if len(l1) != len(l2):
        print("lengths are different")
        return False

    if set(l1) != set(l2):
        print('elements are different')
        return False

    if check_order:
        for i, (o, t) in enumerate(zip(l1, l2)):
            if o != t:
                print('order is different')
                print('at position {0}, list 1 has {1} and 2 has {2}'.format(i, o, t))
                return False

    return True

def compare_dicts(d1, d2, verbose=False):
    if len(d1.keys()) != len(d2.keys()):
        print("lengths are different")
        return False

    if set(d1.keys()) != set(d2.keys()):
        print('elements are different')
        return False

    for o, oE in d1.items():
        if d2[o] != oE:
            print("{0} is {1} and {2}".format(o, oE, d2[o]))
            return False

    return True

def compare_covTs(c1, c2):
    old_scaffs = [x for x, i in c1.items() if len(i) > 0]
    new_scaffs = [x for x, i in c2.items() if len(i) > 0]
    if set(old_scaffs) != set(new_scaffs):
        print('Scaffolds are different')
        return False

    for scaffold in set(old_scaffs):
        s1 = c1[scaffold]
        s2 = c2[scaffold]

        if set(s1.keys()) != set(s2.keys()):
            print("Different mms for scaffold {0}".format(scaffold))
            return False

        for mm, s1_vals in s1.items():
            s2_vals = s2[mm]
            if not s1_vals.equals(s2_vals):
                print("Different values for scaffold {0}".format(scaffold))
                return False

    return True

class test_special():
    def setUp(self):
        self.test_dir = load_random_test_dir()
        self.tearDown()

        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_963_Ns.fasta.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_963_Ns.fasta'
        self.bbsam = load_data_loc() + \
            'bbmap_N5_271_010G1_scaffold_963.fasta.sorted.bam'
        self.bbfasta = load_data_loc() + \
            'N5_271_010G1_scaffold_963.fasta'
        self.IS = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
        self.small_fasta = load_data_loc() + \
            'SmallScaffold.fa'
        self.small_bam = load_data_loc() + \
            'SmallScaffold.fa.sorted.bam'
        self.genes = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.stb = load_data_loc() + \
            'GenomeCoverages.stb'


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
        '''
        Yell if the minor version is not correct
        '''
        assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.3')

        assert not inStrain.SNVprofile.same_versions('1.0.0', '0.8.3')

        assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.4')

        assert not inStrain.SNVprofile.same_versions('0.9.0', '0.8.4')

        assert not inStrain.SNVprofile.same_versions('1.1.0', '0.1.0')

    def test2(self):
        '''
        Make sure it works with Ns in the input
        '''
        # Run base
        base = self.test_dir + 'testNs'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()

    def test3(self):
        '''
        Try other mappers
        '''
        # Run bbmap and let it crash
        base = self.test_dir + 'bbMap'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(self.bbsam, \
            self.bbfasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()
        assert len(MCLdb) == 0

        # Run bbmap and make it work
        base = self.test_dir + 'bbMap'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation --use_full_fasta_header".format(self.bbsam, \
            self.bbfasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()
        assert len(MCLdb) > 0

    def test4(self):
        '''
        Test other --run_statistics
        '''
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


if __name__ == '__main__':
    test_strains().run()
    test_filter_reads().run()
    test_SNVprofile().run()
    test_gene_statistics().run()
    test_quickProfile().run()
    test_genome_wide().run()
    test_plot().run()
    test_special().run()
    test_readcomparer().run()

    print('everything is working swimmingly!')
