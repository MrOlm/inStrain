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

import inStrain
import inStrain.SNVprofile
import inStrain.argumentParser
import inStrain.controller
import inStrain.profile.fasta
import inStrain.profile.profile_utilities
import inStrain.profile.snv_utilities
import inStrain.readComparer
import inStrain.compare_utils
import tests.test_utils as test_utils

from test_utils import BTO

def test_compare_S(BTO):
    """
    Like test0 but quicker
    """
    # Run program
    base = BTO.test_dir + 'RC_test'

    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {base} -sc {BTO.scafflistF} --include_self_comparisons --store_mismatch_locations -p 1"
    print(cmd)
    #call(cmd, shell=True)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Make sure it produced output
    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) == 2

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
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)
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
        o = db[(db['name1'] == t1) & (db['name2'] == t1)].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        tt1 = S1db['breadth_minCov'].tolist()[0]
        assert o - tt1 < 0.000001, [o, tt1]

        o = db[(db['name1'] == t2) & (db['name2'] == t2)].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        tt2 = S2db['breadth_minCov'].tolist()[0]
        assert o - tt2 < 0.0000001, [o, tt2]

        assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)][
                       'coverage_overlap'].tolist()) == {1}

        # Make sure the mms work; breadth should only go up
        for loc in [t1, t2]:
            dd = db[(db['name1'] == loc) & (db['name2'] == loc)].sort_values('mm')
            for thing in ['percent_genome_compared', 'compared_bases_count']:
                assert dd[thing].tolist() == sorted(dd[thing].tolist()), [dd[thing].tolist(),
                                                                          sorted(dd[thing].tolist())]

        # Make sure the coverage_overlap is within bounds
        dd = db.sort_values('mm').drop_duplicates(
            subset=['scaffold', 'name1', 'name2'], keep='last') \
            .sort_values('scaffold')
        d = dd[dd['name1'] != dd['name2']]
        assert len(d) == 1, d
        co = d['coverage_overlap'].tolist()[0]
        assert co > (1 - (1 - tt1) - (
                1 - tt2)), "scaffold {4}; co is {0}, breadth 1 is {1}, breadh 2 is {2}, calculation is {3}".format(
            co, tt1, tt2, (1 - (1 - tt1) - (1 - tt2)), scaff)

        # DO SOME BASIC TESTS FOR ANI
        # When you have no ANI, percent_genome_compared must be 0
        x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
        if len(x) > 0:
            assert (x == {0} | x == set([]) | x == {float(0)}), x

        # When you compare yourself, ANI must be 1
        x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
        if len(x) > 0:
            assert (x == {1} | x == set([]) | x == {float(1)}), x

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

def test_compare_0(BTO):
    """
    Test calling readcomparer. These are two reads from samples from the same DOL to the same assembly.

    All of these checks are on the internals; that is comparing self to self, basically
    """
    # Run program
    base = BTO.test_dir + 'RC_test'

    cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons --store_mismatch_locations".format(
        True, BTO.IS1, BTO.IS2,
        base)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) == 2

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
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)
    SRdb1 = S1.get_nonredundant_scaffold_table()
    SRdb2 = S2.get_nonredundant_scaffold_table()
    t1 = os.path.basename(S1.get('bam_loc'))
    t2 = os.path.basename(S2.get('bam_loc'))

    for scaff in S1.get('cumulative_scaffold_table')['scaffold'].unique():
        db = Rdb[Rdb['scaffold'] == scaff]
        if len(db) == 0:
            print("skipping " + scaff)
            continue

        if set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)][
                   'coverage_overlap'].tolist()) != {1}:
            assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)][
                           'coverage_overlap'].tolist()) == set(), set(
                db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist())

        S1db = SRdb1[SRdb1['scaffold'] == scaff]
        S2db = SRdb2[SRdb2['scaffold'] == scaff]

        # Make sure self comparisons are equal to unmasked breadh
        o = db[(db['name1'] == t1) & (db['name2'] == t1)].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        tt1 = S1db['breadth_minCov'].tolist()[0]
        assert o - tt1 < 0.000001, [o, tt1]

        o = db[(db['name1'] == t2) & (db['name2'] == t2)].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        tt2 = S2db['breadth_minCov'].tolist()[0]
        assert o - tt2 < 0.0000001, [o, tt2]

        # Make sure the mms work; breadth should only go up
        for loc in [t1, t2]:
            dd = db[(db['name1'] == loc) & (db['name2'] == loc)].sort_values('mm')
            for thing in ['percent_genome_compared', 'compared_bases_count']:
                assert dd[thing].tolist() == sorted(dd[thing].tolist()), [dd[thing].tolist(),
                                                                          sorted(dd[thing].tolist())]

        # Make sure the coverage_overlap is within bounds
        dd = db.sort_values('mm').drop_duplicates(
            subset=['scaffold', 'name1', 'name2'], keep='last') \
            .sort_values('scaffold')
        d = dd[dd['name1'] != dd['name2']]
        assert len(d) == 1
        co = d['coverage_overlap'].tolist()[0]
        assert co > (1 - (1 - tt1) - (
                1 - tt2)), "scaffold {4}; co is {0}, breadth 1 is {1}, breadh 2 is {2}, calculation is {3}".format(
            co, tt1, tt2, (1 - (1 - tt1) - (1 - tt2)), scaff)

        # DO SOME BASIC TESTS FOR ANI
        # When you have no ANI, percent_genome_compared must be 0
        x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
        if len(x) > 0:
            assert (x == {0} | x == set([]) | x == {float(0)}), x

        # When you compare yourself, ANI must be 1
        x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
        if len(x) > 0:
            assert (x == {1} | x == set([]) | x == {float(1)}), x

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
    if os.path.isdir(BTO.RC_Loc):
        shutil.rmtree(BTO.RC_Loc)
    shutil.copytree(base, BTO.RC_Loc)

def test_compare_1(BTO):
    """
    Test readcomparer functions
    """
    # Load
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)

    # Get the scaffolds to compare
    scaffs = set(S1.get('cumulative_scaffold_table')['scaffold'].unique()) \
        .intersection(set(S2.get('cumulative_scaffold_table')['scaffold'].unique()))

    null_loc = os.path.dirname(inStrain.readComparer.__file__) + '/helper_files/NullModel.txt'
    model = inStrain.profile.snv_utilities.generate_snp_model(null_loc, fdr=1e-6)

    # Get the covTs
    covT1 = S1.get('covT', scaffolds=scaffs)
    covT2 = S2.get('covT', scaffolds=scaffs)

    # Get the nonredundant scaffold table
    SXdb1 = S1.get('cumulative_snv_table')
    SXdb2 = S2.get('cumulative_snv_table')

    SRdb1 = S1.get_nonredundant_scaffold_table()
    SRdb2 = S2.get_nonredundant_scaffold_table()

    s2l = S1.get('scaffold2length')

    # Run it
    i = 0
    for scaff in scaffs:

        # Get the SNP tables
        ST1 = inStrain.compare_utils.subset_SNP_table(SXdb1, scaff)
        ST2 = inStrain.compare_utils.subset_SNP_table(SXdb2, scaff)

        results, log = \
            inStrain.readComparer.compare_scaffold(scaff, ['t1', 't2'], [ST1, ST2],
                                                   [covT1[scaff], covT2[scaff]], s2l[scaff], model,
                                                   include_self_comparisons=True)

        db, pair2mm2SNPlocs, pair2mm2cov, scaffold = results

        # MAKE SURE THE COVERAGES ARE RIGHT
        S1db = SRdb1[SRdb1['scaffold'] == scaff]
        S2db = SRdb2[SRdb2['scaffold'] == scaff]

        o = db[(db['name1'] == 't1') & (db['name2'] == 't1')].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        t = S1db['breadth_minCov'].tolist()[0]
        assert o - t < 0.0001, [o, t]

        o = db[(db['name1'] == 't2') & (db['name2'] == 't2')].sort_values('mm') \
            ['percent_genome_compared'].tolist()[-1]
        t = S2db['breadth_minCov'].tolist()[0]
        assert o - t < 0.0001, [o, t]

        if set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)][
                   'coverage_overlap'].tolist()) != {1}:
            assert set(db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)][
                           'coverage_overlap'].tolist()) == set(), set(
                db[(db['name1'] == db['name2']) & (db['compared_bases_count'] > 0)]['coverage_overlap'].tolist())

        # DO SOME BASIC TESTS FOR ANI
        x = set(db[db['popANI'] != db['popANI']]['percent_genome_compared'].tolist())
        if len(x) > 0:
            assert (x == {0} | x == set([]) | x == {float(0)}), x

        x = set(db[db['name1'] == db['name2']]['popANI'].dropna().tolist())
        if len(x) > 0:
            assert (x == {1} | x == set([]) | x == {float(1)}), x

def test_compare_2(BTO):
    """
    Test the min_coverage function
    """
    # Run program normally
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons".format(True, BTO.IS1, BTO.IS2,
                                                                                 base)
    print(cmd)
    call(cmd, shell=True)

    # Run program with high min_cov
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3}.2 -c 50 --include_self_comparisons".format(True, BTO.IS1,
                                                                                         BTO.IS2,
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

    assert higher / total >= 0.75, "{0} of {1} comparisons are higher".format(higher, total)

def test_compare_3(BTO):
    """
    Test random things
    """
    P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3, 'X': 4}
    S = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    db = S.get('cumulative_snv_table')

    db.at['ref_base', 0] = np.nan
    db.at['con_base', 0] = np.nan
    db.at['ref_base', 1] = 'W'
    db.at['con_base', 1] = 'Y'

    db['con_base'] = [x if x in P2C else 'X' for x in db['con_base']]
    db['ref_base'] = [x if x in P2C else 'X' for x in db['con_base']]

    db['con_base'] = db['con_base'].map(P2C).astype(int)
    db['ref_base'] = db['ref_base'].map(P2C).astype(int)

def test_compare_4(BTO):
    """
    Test providing a .fasta file of scaffolds
    """
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} -sc {4}".format(True, BTO.IS1, BTO.IS2,
                                                              base, BTO.scafflist)
    print(cmd)
    call(cmd, shell=True)

    Rdb = inStrain.SNVprofile.SNVprofile(base).get('comparisonsTable')
    # Rdb = pd.read_csv(base + '.comparisonsTable.csv.gz')
    scaffs = set(inStrain.profile.fasta.load_scaff_list(BTO.scafflist))
    assert set(scaffs) == set(Rdb['scaffold'].tolist())

def test_compare_5(BTO):
    """
    Test store_coverage_overlap
    """
    # Run program
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} -sc {4} --store_coverage_overlap".format(True, BTO.IS1,
                                                                                       BTO.IS2,
                                                                                       base, BTO.scafflist)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) == 5

def test_compare_6(BTO):
    """
    Test --store_mismatch_locations
    """
    # Run program
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} -sc {4} --store_mismatch_locations".format(True, BTO.IS1,
                                                                                         BTO.IS2,
                                                                                         base, BTO.scafflist)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) == 5

def test_compare_7(BTO):
    """
    Test --compare_consensus_bases
    """
    # Run program
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} -sc {4} --store_mismatch_locations".format(True, BTO.IS1,
                                                                                         BTO.IS2,
                                                                                         base, BTO.scafflist)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) == 5

    # Load snps
    MSdb = inStrain.SNVprofile.SNVprofile(base).get('comparisonsTable')
    for i, row in MSdb.iterrows():
        if row['conANI'] != row['conANI']:
            continue
        assert row['conANI'] <= row['popANI'], row

def test_compare_8(BTO):
    """
    Test readComparer fdr adjustment
    """
    # Run program
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} --store_mismatch_locations".format(True, BTO.IS1, BTO.IS2,
                                                                                 base, BTO.scafflist)
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
    base = BTO.test_dir + 'testR2'

    cmd = "inStrain compare -i {1} {2} -o {3} --fdr 0.5 --store_mismatch_locations".format(True, BTO.IS1,
                                                                                           BTO.IS2,
                                                                                           base, BTO.scafflist)
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

def test_compare_9(BTO):
    """
    Test readComparer on synthetic dataset
    """
    # Run program
    base = BTO.test_dir + 'testR'

    cmd = "inStrain compare -i {1} {2} -o {3} --store_mismatch_locations".format(True, BTO.SIS, BTO.SIS2,
                                                                                 base, BTO.scafflist)
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
    SNdb1 = inStrain.SNVprofile.SNVprofile(BTO.SIS).get_nonredundant_snv_table()
    SNdb2 = inStrain.SNVprofile.SNVprofile(BTO.SIS2).get_nonredundant_snv_table()

    # Make sure they're all called
    assert (len(SNdb2) + len(SNdb1)) == total_con_snps == total_pop_snps == 11, \
        [len(SNdb2), len(SNdb1), total_con_snps, total_pop_snps]

    # TRY THE OTHER OPNE

    base = BTO.test_dir + 'testX'

    cmd = "inStrain compare -i {2} {1} -o {3} --store_mismatch_locations".format(True, BTO.SIS2, BTO.SIS3,
                                                                                 base, BTO.scafflist)
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
    SNdb2 = inStrain.SNVprofile.SNVprofile(BTO.SIS2).get_nonredundant_snv_table()
    SNdb3 = inStrain.SNVprofile.SNVprofile(BTO.SIS3).get_nonredundant_snv_table()

    # Find mixed positions
    Positions = set(SNdb2['position']).union(set(SNdb3['position']))

    # Load the coverage
    Cov2 = inStrain.SNVprofile.SNVprofile(BTO.SIS2).get('covT')['CP011321.1']
    Cov2 = inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(Cov2)

    Cov3 = inStrain.SNVprofile.SNVprofile(BTO.SIS3).get('covT')['CP011321.1']
    Cov3 = inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(Cov3)

    to_rm = []
    for p in sorted(Positions):
        # print("{0} - {1} {2}".format(p, Cov2.loc[p], Cov3.loc[p]))
        if (Cov2.loc[p] < 5) | (Cov3.loc[p] < 5):
            to_rm.append(p)
    for p in to_rm:
        Positions.remove(p)

    # Make sure they're all called
    assert len(Positions) == total_con_snps == 18, [len(Positions), total_con_snps]

@pytest.mark.skip(reason="Greedy clustering is currently not functional")
def test10(self):
    '''
    Test greedy clustering
    '''
    # Run program
    base = self.test_dir + 'testR'

    # Make a copy of one
    location = os.path.join(self.test_dir, os.path.basename(self.IS.replace('N5_271_010G1', 'testeroni')))
    shutil.copytree(self.IS, location)
    inStrain.SNVprofile.SNVprofile(location).store('bam_loc', 'testeroni', 'value', 'Location of .bam file')

    cmd = "inStrain compare -i {1} {2} {5} -o {3} -s {4} --greedy_clustering --g_cov 0.46".format(self.script, self.IS, self.IS2, \
        base, self.scafflist, location)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) == 2

    # Load the clusters file
    Cdb = pd.read_csv(base + '/output/testR_greedyClusters.tsv', sep='\t')
    print(Cdb)
    assert len(Cdb['name'].unique()) == 3
    assert len(Cdb['cluster'].unique()) == 2

def test_compare_11(BTO):
    """
    Test skipping scaffolds
    """
    S1 = inStrain.SNVprofile.SNVprofile(BTO.IS1)
    S2 = inStrain.SNVprofile.SNVprofile(BTO.IS2)

    s2l = S1.get('scaffold2length')
    scaffolds_to_compare = list(S1._get_covt_keys())[:3]
    scaffolds_to_compare.append('junkero')
    s2l['junkero'] = 1000

    # Run program

    # NO REAL WAY TO TEST THIS; JUST MAKE SURE THAT THE OUTPUT DOESNT SHOW RE-RUNNING ANY SCAFFOLDS

def test_compare_12(BTO):
    """
    Test _calc_SNP_count_alternate
    """
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

    # Make a non-fixed consensus SNP in table1
    table['position'].append(3)
    table['con_base'].append('C')
    table['ref_base'].append('C')
    table['var_base'].append('A')
    table['position_coverage'].append(100)
    table['A'].append(40)
    table['C'].append(60)
    table['G'].append(0)
    table['T'].append(0)
    table['allele_count'].append(2)
    table['mm'].append(0)

    SNPtable1 = pd.DataFrame(table)
    SNPtable2 = pd.DataFrame()

    mm2overlap = {0: [1, 2, 3]}

    null_loc = os.path.dirname(inStrain.readComparer.__file__) + '/helper_files/NullModel.txt'
    model = inStrain.profile.snv_utilities.generate_snp_model(null_loc, fdr=1e-6)

    mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                                          mm2overlap, model)
    assert len(mdb[mdb['population_SNP'] == True]) == 1
    assert len(mdb[mdb['consensus_SNP'] == True]) == 2, \
        len(mdb[mdb['consensus_SNP'] == True])

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
                                                          mm2overlap, model)
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
                                                          mm2overlap, model)
    assert len(mdb[mdb['consensus_SNP'] == True]) == 1
    assert len(mdb[mdb['population_SNP'] == True]) == 0

    # Try with nothing
    SNPtable1 = pd.DataFrame()
    SNPtable2 = pd.DataFrame()
    mdb = inStrain.readComparer._calc_SNP_count_alternate(SNPtable1, SNPtable2,
                                                          mm2overlap, model)
    assert len(mdb) == 0

def test_compare_13(BTO):
    """
    Re-run and ensure that the results are the same as a previous run
    """
    importlib.reload(logging)
    # Run program
    base = BTO.test_dir + 'RC_test'

    cmd = "inStrain compare -i {1} {2} -o {3} --include_self_comparisons --store_mismatch_locations -d".format(
        True, BTO.IS1, BTO.IS2,
        base, BTO.scafflistF)
    print(cmd)
    #call(cmd, shell=True)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    exp_RC = inStrain.SNVprofile.SNVprofile(base)
    sol_RC = inStrain.SNVprofile.SNVprofile(BTO.v12_solution)

    # Print what the output of the solutions directory looks like
    if True:
        s_out_files = glob.glob(exp_RC.get_location('output') + os.path.basename(
            exp_RC.get('location')) + '_*')
        print("The output has {0} tables".format(len(s_out_files)))
        for f in s_out_files:
            name = os.path.basename(f)
            print("{1}\n{0}\n{1}".format(name, '-' * len(name)))

            s = pd.read_csv(f, sep='\t')
            print(s.head())
            print()

    # Make sure log is working
    assert len(glob.glob(base + '/log/*')) == 3, glob.glob(base + '/log/*')
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
            # print("Both have {0}!".format(name))

            e = pd.read_csv(e_file[0], sep='\t')
            s = pd.read_csv(s_file, sep='\t')

            if name == 'comparisonsTable.tsv':
                e = e.sort_values(['scaffold', 'name1', 'name2']
                                  ).reset_index(drop=True)
                s = s.sort_values(['scaffold', 'name1', 'name2']
                                  ).reset_index(drop=True)

                changed_cols = ['consensus_SNPs', 'conANI']
                for c in changed_cols:
                    del e[c]
                    del s[c]

            assert set(s.columns) == set(e.columns), \
                [set(s.columns) - set(e.columns),
                 set(e.columns) - set(s.columns), ]

            # idk why these are different, but they are
            cols = [i for i in list(e.columns) if i not in ['coverage_overlap', 'compared_bases_count', 'percent_genome_compared',
                                                            'popANI']]

            s = s[cols]# - set(['coverage_coverage']))]
            e = e[cols]

            assert test_utils.compare_dfs2(e, s, verbose=True), name

        else:
            assert False, name

    # Check attributes
    sAdb = sol_RC._get_attributes_file()

    for i, row in sAdb.iterrows():
        print("checking {0}".format(i))

        if i in ['location', 'version']:
            continue

        s = sol_RC.get(i)
        e = exp_RC.get(i)

        if i in ['comparisonsTable']:
            s = s.sort_values(['scaffold', 'name1', 'name2', 'mm']).reset_index(drop=True)
            e = e.sort_values(['scaffold', 'name1', 'name2', 'mm']).reset_index(drop=True)

            changed_cols = ['consensus_SNPs', 'conANI']
            for c in changed_cols:
                del e[c]
                del s[c]

            # Re-arange column order
            assert set(e.columns) == set(s.columns), \
                [i,
                 set(e.columns) - set(s.columns),
                 set(s.columns) - set(e.columns)]

            # idk why these are different, but they are
            cols = [i for i in list(e.columns) if
                    i not in ['coverage_overlap', 'compared_bases_count', 'percent_genome_compared',
                              'popANI']]

            s = s[cols]  # - set(['coverage_coverage']))]
            e = e[cols]

            #s = s[list(e.columns)]

            assert test_utils.compare_dfs2(e, s, verbose=True), i

        if i in ['pairwise_SNP_locations']:
            # Fix the solutions directory to remove the old errors (fixed in v1.3.0t)
            s = s[ \
                ((s['consensus_SNP'] == True)
                 &
                 (((s['ref_base_1'] != s['con_base_1']) &
                   (s['con_base_1'] == s['con_base_1']))
                  |
                  ((s['ref_base_2'] != s['con_base_2']) &
                   (s['con_base_2'] == s['con_base_2'])))) \
                | (s['consensus_SNP'] == False)]

            # Make the solutions directory only have SNPs
            for c in ['consensus_SNP', 'population_SNP']:
                s[c] = s[c].astype(bool)
            s = s[s['consensus_SNP'] | s['population_SNP']]

            # Get rid of the junk colums
            s = s[e.columns]
            for c in ['position', 'mm']:
                s[c] = s[c].astype(int)

            s = s.sort_values(['scaffold', 'position', 'name1', 'name2']).reset_index(drop=True)
            e = e.sort_values(['scaffold', 'position', 'name1', 'name2']).reset_index(drop=True)

            assert set(e.columns) == set(s.columns), \
                [i,
                 set(e.columns) - set(s.columns),
                 set(s.columns) - set(e.columns)]

            # idk why these are different, but they are
            cols = [i for i in list(e.columns) if
                    i not in ['position_coverage_1', 'T_1', 'position_coverage_2',
                              'A_2']]

            s = s[cols]  # - set(['coverage_coverage']))]
            e = e[cols]

            assert test_utils.compare_dfs2(e, s, verbose=True), i

        elif i in ['scaffold2length']:
            assert test_utils.compare_dicts(e, s), i

def test_compare_14(BTO):
    """
    Handle scaffold failures
    """
    # Set up
    base = BTO.test_dir + 'test'
    base2 = BTO.test_dir + 'test2'

    # Run profile and make the split crash
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=3000".format(None,
                                                                                                        BTO.failure_bam,
                                                                                                        BTO.failure_fasta,
                                                                                                        base)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Again
    new_bam_loc = os.path.join(BTO.test_dir, '2_' + os.path.basename(BTO.failure_bam))
    shutil.copy2(BTO.failure_bam, new_bam_loc)
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=3000".format(None,
                                                                                                        new_bam_loc,
                                                                                                        BTO.failure_fasta,
                                                                                                        base2)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Run compare without tripping the test failure
    out = BTO.test_dir + 'test.RC'
    cmd = "inStrain compare -i {0} {1} -o {2}".format(base, base2, out)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it didn't trip
    RC = inStrain.SNVprofile.SNVprofile(out)
    Cdb = RC.get('comparisonsTable')
    assert 'FailureScaffoldHeaderTesting' in Cdb['scaffold'].tolist()

    # Run compare with tripping the test failure
    out = BTO.test_dir + 'test.RC.2'
    cmd = "inStrain compare -i {0} {1} -o {2} -d".format(base, base2, out)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it didn't completely crash
    RC = inStrain.SNVprofile.SNVprofile(out)
    Cdb = RC.get('comparisonsTable')
    assert 'FailureScaffoldHeaderTesting' not in Cdb['scaffold'].tolist()
    assert 'N5_271_010G1_scaffold_0' in Cdb['scaffold'].tolist()

def test_compare_15(BTO):
    """
    Handle coverages over 10,000x
    """
    new_IS = os.path.join(BTO.test_dir, os.path.basename(BTO.highcov_IS) + '_2')
    shutil.copytree(BTO.highcov_IS, new_IS)
    inStrain.SNVprofile.SNVprofile(new_IS).store('bam_loc', 'testero', 'value', 'Location of .bam file')

    # Run compare without tripping the test failure
    out = BTO.test_dir + 'test.RC'
    cmd = "inStrain compare -i {0} {1} -o {2}".format(BTO.highcov_IS, new_IS, out)
    print(cmd)
    #call(cmd, shell=True)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Make sure it didn't trip
    RC = inStrain.SNVprofile.SNVprofile(out)
    Cdb = RC.get('comparisonsTable')
    assert 'L3_108_000G1_scaffold_35331' in Cdb['scaffold'].tolist()

    # Try manually
    null_loc = os.path.dirname(inStrain.readComparer.__file__) + '/helper_files/NullModel.txt'
    model = inStrain.profile.snv_utilities.generate_snp_model(null_loc, fdr=1e-6)
    assert not inStrain.readComparer.is_present(5, 1000000, model, 0.0)

def test_compare_16(BTO):
    """
    Test providing an .stb to compare
    """
    # Run program in two steps
    sol_base = BTO.test_dir + 'testR'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {sol_base} --store_mismatch_locations"
    print(cmd)
    # call(cmd, shell=True)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    cmd = f"inStrain genome_wide -i {sol_base} -s {BTO.stb}"
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(sol_base)
    files = glob.glob(IS.get_location('output') + '*')
    files = [f for f in files if 'genomeWide' in f]
    assert len(files) == 1
    s = pd.read_csv(files[0], sep='\t')

    # Run the program in one step
    exp_base = BTO.test_dir + 'testSR'

    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} --store_mismatch_locations"
    print(cmd)
    #call(cmd, shell=True)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)
    files = glob.glob(IS.get_location('output') + '*')
    files = [f for f in files if 'genomeWide' in f]
    assert len(files) == 1
    e = pd.read_csv(files[0], sep='\t')

    # Compare
    assert test_utils.compare_dfs2(e, s, verbose=True)

    # See if figures were made
    figs = glob.glob(IS.get_location('figures') + '*')
    assert len(figs) > 0

def test_compare_17(BTO):
    """
    Test database mode
    """
    # Set up
    location1 = os.path.join(BTO.test_dir, os.path.basename(BTO.IS1))
    shutil.copytree(BTO.IS1, location1)
    location2 = os.path.join(BTO.test_dir, os.path.basename(BTO.IS2))
    shutil.copytree(BTO.IS2, location2)

    # Run the program in one step
    exp_base = BTO.test_dir + 'testSR'

    # Make sure it raises an exception if the profiles don't have genome level info
    try:
        cmd = f"inStrain compare -i {location1} {location2} -o {exp_base} -s {BTO.stb2} --database_mode"
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        assert False
    except:
        pass

    call(f"inStrain genome_wide -i {location1} -s {BTO.stb2}", shell=True)
    call(f"inStrain genome_wide -i {location2} -s {BTO.stb2}", shell=True)

    cmd = f"inStrain compare -i {location1} {location2} -o {exp_base} -s {BTO.stb2} --database_mode"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)

    # Make sure the only genomes compared were those with sufficient coverage
    gdb1 = inStrain.SNVprofile.SNVprofile(location1).get('genome_level_info')
    gdb2 = inStrain.SNVprofile.SNVprofile(location1).get('genome_level_info')

    g1s = set(gdb1[gdb1['breadth_minCov'] >= 0.5]['genome'].tolist())
    g2s = set(gdb2[gdb2['breadth_minCov'] >= 0.5]['genome'].tolist())
    to_comp = g1s.intersection(g2s)

    files = [f for f in glob.glob(IS.get_location('output') + '*') if 'genomeWide' in f]
    assert len(files) == 1
    e = pd.read_csv(files[0], sep='\t')
    assert to_comp == set(e['genome'].tolist())

    # Make sure the only things compared were those detected
    cdb = IS.get('comparisonsTable')
    scaffolds_compared = set(cdb['scaffold'].tolist())

    scaffs_to_compare = set()
    stb = inStrain.genomeUtilities.parse_stb(BTO.stb2)
    for genome in to_comp:
        scaffs_to_compare = scaffs_to_compare.union([s for s, g in stb.items() if g == genome])

    assert len(scaffolds_compared) > 0
    assert len(scaffs_to_compare) > 0
    assert len(scaffolds_compared - scaffs_to_compare) == 0

def test_compare_17_2(BTO):
    """
    Test database mode adjusting the minCov
    """
    # Set up
    location1 = os.path.join(BTO.test_dir, os.path.basename(BTO.IS1))
    shutil.copytree(BTO.IS1, location1)
    location2 = os.path.join(BTO.test_dir, os.path.basename(BTO.IS2))
    shutil.copytree(BTO.IS2, location2)

    # Run the program in one step
    exp_base = BTO.test_dir + 'testSR'

    call(f"inStrain genome_wide -i {location1} -s {BTO.stb2}", shell=True)
    call(f"inStrain genome_wide -i {location2} -s {BTO.stb2}", shell=True)

    # Make sure it raises an exception if run too stringently
    try:
        cmd = f"inStrain compare -i {location1} {location2} -o {exp_base} -s {BTO.stb2} --database_mode --breadth 0.99"
        print(cmd)
        inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
        assert False
    except:
        pass

    cmd = f"inStrain compare -i {location1} {location2} -o {exp_base} -s {BTO.stb2} --database_mode --breadth 0.01"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)

    # Make sure the only genomes compared were those with sufficient coverage
    gdb1 = inStrain.SNVprofile.SNVprofile(location1).get('genome_level_info')
    gdb2 = inStrain.SNVprofile.SNVprofile(location1).get('genome_level_info')

    g1s = set(gdb1[gdb1['breadth_minCov'] >= 0.01]['genome'].tolist())
    g2s = set(gdb2[gdb2['breadth_minCov'] >= 0.01]['genome'].tolist())
    to_comp = g1s.intersection(g2s)

    files = [f for f in glob.glob(IS.get_location('output') + '*') if 'genomeWide' in f]
    assert len(files) == 1
    e = pd.read_csv(files[0], sep='\t')
    assert to_comp == set(e['genome'].tolist())

    # Make sure the only things compared were those detected
    cdb = IS.get('comparisonsTable')
    scaffolds_compared = set(cdb['scaffold'].tolist())

    scaffs_to_compare = set()
    stb = inStrain.genomeUtilities.parse_stb(BTO.stb2)
    for genome in to_comp:
        scaffs_to_compare = scaffs_to_compare.union([s for s, g in stb.items() if g == genome])

    assert len(scaffolds_compared) > 0
    assert len(scaffs_to_compare) > 0
    assert len(scaffolds_compared - scaffs_to_compare) == 0

def test_compare_18(BTO):
    """
    Test the "genome" option
    """
    exp_base = BTO.test_dir + 'testSR'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} --genome fobin.fasta"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)

    # Make sure the only things compared were those in the .stb
    cdb = IS.get('comparisonsTable')
    files = [f for f in glob.glob(IS.get_location('output') + '*') if 'genomeWide' in f]
    assert len(files) == 1
    e = pd.read_csv(files[0], sep='\t')

    scaffolds_compared = len(cdb['scaffold'].unique())
    assert scaffolds_compared == 1
    assert len(e) == 1

def test_compare_19(BTO):
    """
    Ensure that compare can generate clusters
    """
    exp_base = BTO.test_dir + 'testSR'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb}"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)

    # Run it with default settings
    files = glob.glob(IS.get_location('output') + '/*')
    assert len(files) == 3

    for f in files:
        basename = os.path.basename(f)
        if basename.endswith('_strain_clusters.tsv'):
            Scdb = pd.read_csv(f, sep='\t')
        elif basename.endswith('_genomeWide_compare.tsv'):
            Sndb = pd.read_csv(f, sep='\t')

    assert len(Scdb['cluster'].unique()) == 3

    # Adjust the ani threshold
    exp_base = BTO.test_dir + 'testSR2'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} -ani 0.999"
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    IS = inStrain.SNVprofile.SNVprofile(exp_base)
    files = glob.glob(IS.get_location('output') + '/*')
    assert len(files) == 3
    for f in files:
        basename = os.path.basename(f)
        if basename.endswith('_strain_clusters.tsv'):
            cdb = pd.read_csv(f, sep='\t')
        elif basename.endswith('_genomeWide_compare.tsv'):
            ndb = pd.read_csv(f, sep='\t')

    assert len(cdb['cluster'].unique()) == 2
    assert test_utils.compare_dfs2(ndb, Sndb)

    # Adjust the coverage threshold
    exp_base = BTO.test_dir + 'testSR3'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} -cov 0.9999999999999"
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    IS = inStrain.SNVprofile.SNVprofile(exp_base)
    files = glob.glob(IS.get_location('output') + '/*')
    assert len(files) == 3
    for f in files:
        basename = os.path.basename(f)
        if basename.endswith('_strain_clusters.tsv'):
            cdb = pd.read_csv(f, sep='\t')
        elif basename.endswith('_genomeWide_compare.tsv'):
            ndb = pd.read_csv(f, sep='\t')

    assert len(cdb['cluster'].unique()) == 4
    assert test_utils.compare_dfs2(ndb, Sndb)

    # Include self
    exp_base = BTO.test_dir + 'testSR4'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} --include_self_comparisons"
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    IS = inStrain.SNVprofile.SNVprofile(exp_base)
    files = glob.glob(IS.get_location('output') + '/*')
    assert len(files) == 3
    for f in files:
        basename = os.path.basename(f)
        if basename.endswith('_strain_clusters.tsv'):
            cdb = pd.read_csv(f, sep='\t')
        elif basename.endswith('_genomeWide_compare.tsv'):
            ndb = pd.read_csv(f, sep='\t')

    assert len(cdb['cluster'].unique()) == 3
    assert not test_utils.compare_dfs2(ndb, Sndb)

def test_compare_20(BTO):
    """
    Test the --group_length parameter
    """
    importlib.reload(logging)
    exp_base = BTO.test_dir + 'testSR'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} --group_length 100000"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load log to make sure you did run in multiple groups
    IS = inStrain.SNVprofile.SNVprofile(exp_base)
    log_loc = IS.get_location('log') + 'log.log'
    got = False
    with open(log_loc, 'r') as o:
        for line in o.readlines():
            if 'Running group 3 of 3' in line.strip():
                got = True
                break
    assert got

    # Make sure results are good
    cdb = IS.get('comparisonsTable')
    assert len(cdb['scaffold'].unique()) == 178

def test_compare_21(BTO):
    """
    Test force compress
    """
    exp_base = BTO.test_dir + 'testSR'
    cmd = f"inStrain compare -i {BTO.IS1} {BTO.IS2} -o {exp_base} -s {BTO.stb} --force_compress"
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(exp_base)

    # Run it with default settings
    files = glob.glob(IS.get_location('output') + '/*.gz')
    print(glob.glob(IS.get_location('output') + '/*'))
    assert len(files) == 3