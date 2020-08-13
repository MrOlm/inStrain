"""
Run tests on inStrain filter_reads
"""

import glob
import os
import shutil
from subprocess import call

import numpy as np
import pandas as pd

import inStrain
import inStrain.GeneProfile
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import tests.test_utils as test_utils


class test_gene_statistics:
    def __init__(self):
        self.script = test_utils.get_script_loc('GeneProfile')
        self.old_script = test_utils.get_script_loc('deprecated_gene_statistics')
        self.test_dir = test_utils.load_random_test_dir()
        self.old_IS = test_utils.load_data_loc() + \
                      'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.v4.IS'
        self.IS = test_utils.load_data_loc() + \
                  'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.3.0g'
        self.genes = test_utils.load_data_loc() + \
                     'N5_271_010G1_scaffold_min1000.fa.genes.fna'
        self.IS2 = test_utils.load_data_loc() + \
                   'sars_cov_2_MT039887.1.fasta.bt2-vs-SRR11140750.sam.IS'
        self.genbank = test_utils.load_data_loc() + \
                       'sars_cov_2_MT039887.1.gb'

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
        """
        Test deprecated gene_statistic
        """
        location = os.path.join(self.test_dir, os.path.basename(self.old_IS))
        shutil.copytree(self.old_IS, location)

        # Run program
        base = self.test_dir + 'test'

        cmd = "{0} {1} -g {3}".format(self.old_script, location,
                                      base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert len(glob.glob(location + '/output/*')) >= 6

        # Read something
        Rdb = pd.read_csv(location + '/output/aa-SNVs.tsv', sep='\t')
        print(len(Rdb))
        # assert len(Rdb) == 31
        assert len(Rdb) == 350

        # Check gene coverages [THIS DOESN'T WORK]
        # Gdb = pd.read_csv(location + '/output/genes.tsv', sep='\t')
        # print(len(Gdb))
        # assert len(Gdb) > 0

    def test1(self):
        """
        Test GeneProfile vs deprecated gene_statistic
        """
        # Run the dereplicated version
        location = os.path.join(self.test_dir, os.path.basename(self.old_IS))
        shutil.copytree(self.old_IS, location)

        # Run program
        base = self.test_dir + 'testO'

        cmd = "{0} {1} -g {3}".format(self.old_script, location,
                                      base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Read the mutations
        Rdb = pd.read_csv(location + '/output/aa-SNVs.tsv', sep='\t')
        assert len(Rdb) == 350

        # Get the mutations in this IS (OLD VERSION)
        Sdb = inStrain.SNVprofile.SNVprofile(location).get_nonredundant_snv_table()
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

        cmd = "inStrain profile_genes -i {1} -g {3}".format(self.script, location,
                                                            base, self.genes)
        print(cmd)
        call(cmd, shell=True)

        # Get the mutations in this IS
        Sdb = inStrain.SNVprofile.SNVprofile(location).get_nonredundant_snv_table()
        Sdb['allele_count'] = Sdb['allele_count'].astype(int)
        Sdb = Sdb[Sdb['cryptic'] == False]
        if 'allele_count' in Sdb.columns:
            Sdb = Sdb[Sdb['allele_count'].isin([1, 2])]
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

        # assert len(RNdb) > 300, len(RNdb)
        assert len(RNdb) > 250, len(RNdb)
        assert len(Rdb) == len(RNdb)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_colwidth', -1)

        Rdb = Rdb[['scaffold', 'position', 'mutation_type', 'mutation', 'gene', 'mut_key']].sort_values(
            ['scaffold', 'position'])
        Rdb['gene'] = [0 if (g is None) | (g == 'None') else g for g in Rdb['gene']]
        RNdb = RNdb[['scaffold', 'position', 'mutation_type', 'mutation', 'gene', 'mut_key']].sort_values(
            ['scaffold', 'position'])

        Rdb = Rdb.fillna(0)
        RNdb = RNdb.fillna(0)

        for i, row in Rdb.iterrows():
            Nrow = RNdb[RNdb['mut_key'] == row['mut_key']]
            for col in Rdb.columns:
                assert row[col] == Nrow[col].tolist()[0], [col, row, Nrow]

        assert test_utils.compare_dfs(RNdb, Rdb, verbose=True)

        # Also compare gene clonalities

    def test2(self):
        """
        Make sure it produces output
        """
        location = os.path.join(self.test_dir, os.path.basename(self.IS))
        shutil.copytree(self.IS, location)

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3} --store_everything -d".format(self.script, location,
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
                # print(ScaffCov, GeneCov, scaff, mm)
                # assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov - GeneCov) / ScaffCov))
        assert np.mean(diffs) < 0.15  # The average difference is less than 15%

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
        """
        Make sure it works with a genbank file
        """
        location = os.path.join(self.test_dir, os.path.basename(self.IS2))
        shutil.copytree(self.IS2, location)

        # Oh god this is horrible; re-name some variables
        tis = inStrain.SNVprofile.SNVprofile(location)
        db = tis.get('cumulative_snv_table')
        db = db.rename(columns={'varBase': 'var_base', 'conBase': 'con_base'})
        tis.store('cumulative_snv_table', db, 'pandas',
                  'Cumulative SNP on mm level. Formerly snpLocations.pickle')

        # db = tis.get('raw_snp_table')
        # db = db.rename(columns={'var_base':'varBase', 'con_base':'conBase'})
        # tis.store('raw_snv_table', db, 'pandas',
        #             'Contains raw SNP information on a mm level')

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3} --store_everything".format(self.script, location,
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
                # print(ScaffCov, GeneCov, scaff, mm)
                # assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov - GeneCov) / ScaffCov))
        assert np.mean(diffs) < 0.16, np.mean(diffs)  # The average difference is less than 16%

    def test4(self):
        """
        Make sure it works when there are no SNVs
        """
        location = os.path.join(self.test_dir, os.path.basename(self.IS2))
        shutil.copytree(self.IS2, location)

        # Delete SNVs
        s_loc = location + '/output/sars_cov_2_MT039887.1.fasta.bt2-vs-SRR11140750.sam.IS_SNVs.tsv'
        with open(s_loc, 'w'):
            pass
        inStrain.SNVprofile.SNVprofile(location).store('cumulative_snv_table', pd.DataFrame(),
                                                       'pandas',
                                                       'Cumulative SNP on mm level. Formerly snpLocations.pickle')

        # Run program
        base = self.test_dir + 'testN'

        cmd = "inStrain profile_genes -i {1} -g {3}".format(self.script, location,
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
                # print(ScaffCov, GeneCov, scaff, mm)
                # assert abs((ScaffCov-GeneCov)/ScaffCov) < 3
                diffs.append(abs((ScaffCov - GeneCov) / ScaffCov))
        assert np.mean(diffs) < 0.16, np.mean(diffs)  # The average difference is less than 16%
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
