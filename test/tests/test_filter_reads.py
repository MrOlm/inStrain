"""
Run tests on inStrain filter_reads
"""

import glob
import importlib
import logging
import os
import shutil
from subprocess import call
import numpy as np

import pandas as pd
from Bio import SeqIO

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import tests.test_utils as test_utils


# noinspection PyTypeChecker
class TestFilterReads:
    def set_up(self):
        self.script = test_utils.get_script_loc('filter_reads')
        self.test_dir = test_utils.load_random_test_dir()
        self.sorted_bam = test_utils.load_data_loc() + \
                          'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa'
        self.readloc = test_utils.load_data_loc() + \
                       'N5_271_010G1.R1.fastq.gz'
        self.readloc2 = test_utils.load_data_loc() + \
                        'N5_271_010G1.R1.reads'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

        importlib.reload(logging)

    def tear_down(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        importlib.reload(logging)

    def run(self, min=0, max=4, tests='all'):
        if tests == 'all':
            tests = np.arange(min, max)

        for test_num in tests:
            self.set_up()
            print("\n*** Running test {0} ***\n".format(test_num))
            eval('self.test{0}()'.format(test_num))
            self.tear_down()

    def test0(self):
        """
        Test basic functionality
        """
        # Run program
        base = self.test_dir + 'test'

        cmd = "inStrain filter_reads {0} {1} -o {2}".format(self.sorted_bam,
                                                            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        files = glob.glob(base + '/*')
        assert len(files) == 1
        Rdb = pd.read_csv(files[0], sep='\t', header=1)
        assert len(Rdb) == 179, len(Rdb)

        # Make it make the detailed read report
        cmd = "inStrain filter_reads {0} {1} -o {2} --detailed_mapping_info".format(self.sorted_bam,
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
        """
        Compare version 1 and 2 of getting paired reads
        """
        # Load
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))  # set up .fasta file
        scaffolds = list(scaff2sequence.keys())

        # Run version 2 to get scaffold to pairs to info
        s2pair2info2 = inStrain.filter_reads.get_paired_reads_multi(self.sorted_bam, scaffolds)
        db2 = inStrain.filter_reads.make_detailed_mapping_info(s2pair2info2)

        # Run version 1 to get scaffold to pairs to info
        s2pair2info, scaff2total = inStrain.deprecated.deprecated_filter_reads.get_paired_reads_multi(self.sorted_bam,
                                                                                                      scaffolds,
                                                                                                      ret_total=True)
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
        kwargs = {'pairing_filter': 'non_discordant'}
        pair2info2N = inStrain.filter_reads.paired_read_filter(s2pair2info2, **kwargs)
        # filtered_2N = set(p for p, i in pair2info2N.items())
        filtered_2N = set([p for s, p2i in pair2info2N.items() for p, i in p2i.items()])
        assert filtered_1 != filtered_2N, [len(filtered_1), len(filtered_2N)]

        kwargs = {'pairing_filter': 'all_reads'}
        pair2info2N = inStrain.filter_reads.paired_read_filter(s2pair2info2, **kwargs)
        # filtered_2A = set(p for p, i in pair2info2N.items())
        filtered_2A = set([p for s, p2i in pair2info2N.items() for p, i in p2i.items()])
        assert filtered_2N != filtered_2A, [len(filtered_2N), len(filtered_2A)]
        assert filtered_1 != filtered_2A, [len(filtered_1), len(filtered_2A)]

        # Make sure that not priority reads changes this
        priority_reads_loc = self.readloc
        priority_reads = inStrain.filter_reads.load_priority_reads(priority_reads_loc)
        kwargs = {}
        pair2info2P = inStrain.filter_reads.paired_read_filter(s2pair2info2, priority_reads_set=priority_reads,
                                                               **kwargs)
        filtered_2P = set(p for p, i in pair2info2P.items())
        assert filtered_2 != filtered_2P, [len(filtered_2), len(filtered_2P)]

        # Make the old-style read report
        # Merge into single
        pair2info = {}
        for scaff, p2i in s2pair2info.items():
            for p, i in p2i.items():
                pair2info[p] = i
        RRo = inStrain.deprecated.deprecated_filter_reads.makeFilterReport(s2pair2info, scaff2total=scaff2total)

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
        assert RRo['unfiltered_pairs'].tolist()[0] == RR['unfiltered_pairs'].tolist()[0] == len(filtered_2) == len(
            filtered_1), \
            [RRo['unfiltered_pairs'].tolist()[0], RR['unfiltered_pairs'].tolist()[0], len(filtered_2), len(filtered_1)]

        o2n = {'pass_min_read_ani': 'pass_filter_cutoff'}
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
        RRs = inStrain.filter_reads.makeFilterReport2(s2pair2info2, scaffold_level_mapping_info=True,
                                                      pairTOinfo=pair2info2, priority_reads_set=priority_reads,
                                                      **kwargs)
        for col in ['singletons']:
            assert RRs['unfiltered_' + col].tolist()[0] > 0
        assert RRs['unfiltered_priority_reads'].tolist()[0] > 0

    def test2(self):
        """
        Compare method 1 and 2 of getting filtered reads
        """
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))  # set up .fasta file
        scaffolds = list(scaff2sequence.keys())

        # Try new method
        pair2infoF, RR = inStrain.filter_reads.load_paired_reads(self.sorted_bam, scaffolds)
        assert len(set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()])) \
               == int(RR['filtered_pairs'].tolist()[0])

        # Try old method

        # Load paired reads
        scaff2pair2info, scaff2total = inStrain.deprecated.deprecated_filter_reads.get_paired_reads_multi(
            self.sorted_bam, scaffolds, ret_total=True)
        # Merge into single
        pair2info = {}
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                pair2info[p] = i
        # Make read report
        logging.info('Making read report')
        RRo = inStrain.deprecated.deprecated_filter_reads.makeFilterReport(scaff2pair2info, scaff2total,
                                                                           pair2info=pair2info)
        # Filter the dictionary
        logging.info('Filtering reads')
        pair2infoFo = inStrain.deprecated.deprecated_filter_reads.filter_paired_reads_dict(pair2info)
        assert len(pair2infoFo.keys()) == int(RRo['filtered_pairs'].tolist()[0])

        # Compare
        assert set(pair2infoFo.keys()) == \
               set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()])
        for s, p2i in pair2infoF.items():
            for pair, info in p2i.items():
                assert info == pair2infoFo[pair], pair

        # Try new method with priority_reads
        kwargs = {"priority_reads": self.readloc2}
        pair2infoF, RR = inStrain.filter_reads.load_paired_reads(self.sorted_bam, scaffolds, **kwargs)

        pair2infoF_pairs = set([p for s, p2i in pair2infoF.items() for p, i in p2i.items()])

        PReads = inStrain.filter_reads.load_priority_reads(self.readloc2)
        assert set(pair2infoFo.keys()) != pair2infoF_pairs
        assert len(pair2infoF_pairs - set(pair2infoFo.keys())) > 0
        assert len(pair2infoF_pairs - set(pair2infoFo.keys()) - set(PReads)) == 0, \
            len(set(pair2infoF.keys()) - set(pair2infoFo.keys()) - set(PReads))

    def test3(self):
        """
        Test read filtering options
        """
        # Run on default
        base = self.test_dir + 'test'
        cmd = "inStrain filter_reads {0} {1} -o {2}".format(self.sorted_bam,
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
        cmd = "inStrain filter_reads {0} {1} -o {2} --pairing_filter all_reads".format(self.sorted_bam,
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
        cmd = "inStrain filter_reads {0} {1} -o {2} --pairing_filter non_discordant".format(self.sorted_bam,
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
        cmd = "inStrain filter_reads {0} {1} -o {2} --priority_reads {3}".format(self.sorted_bam,
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
        """
        Test that chaning the read filtering options actually changes the results
        """
        # Run base
        base = self.test_dir + 'testMatt'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(
            self.sorted_bam,
            self.fasta, base)
        call(cmd, shell=True)

        # Load the object
        Matt_object = inStrain.SNVprofile.SNVprofile(base)
        MCLdb = Matt_object.get_nonredundant_scaffold_table()

        # Run with priority_reads
        base = self.test_dir + 'testMatt'
        cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation --priority_reads {3}".format(
            self.sorted_bam,
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
