"""
Run tests on inStrain profile
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
import inStrain.argumentParser
import inStrain.controller
import inStrain.deprecated
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.logUtils
import inStrain.profile
import inStrain.profile.fasta
import tests.test_utils as test_utils
import inStrain.polymorpher

from test_utils import BTO

# class test_profile:
#     def __init__(BTO):
#         BTO.script = test_utils.get_script_loc('inStrain')
#         BTO.test_dir = test_utils.load_random_test_dir()
#         BTO.bam1 = test_utils.load_data_loc() + \
#                           'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
#         BTO.fasta = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_min1000.fa'
#         BTO.failure_bam = test_utils.load_data_loc() + \
#                            'N5_271_010G1_scaffold_failureScaffold.sorted.bam'
#         BTO.single_scaff = test_utils.load_data_loc() + \
#                             'N5_271_010G1_scaffold_101.fasta'
#         BTO.fasta_extra = test_utils.load_data_loc() + \
#                            'N5_271_010G1_scaffold_min1000_extra.fa'
#         BTO.small_fasta = test_utils.load_data_loc() + \
#                            'SmallScaffold.fa'
#         BTO.small_bam = test_utils.load_data_loc() + \
#                          'SmallScaffold.fa.sorted.bam'
#         BTO.extra_single_scaff = test_utils.load_data_loc() + \
#                                   'N5_271_010G1_scaffold_101_extra.fasta'
#         BTO.failure_fasta = test_utils.load_data_loc() + \
#                              'N5_271_010G1_scaffold_failureScaffold.fa'
#         BTO.failure_genes = test_utils.load_data_loc() + \
#                              'N5_271_010G1_scaffold_failureScaffold.fa.genes.fna.fa'
#         BTO.cc_solution = test_utils.load_data_loc() + \
#                            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam.CB'
#         BTO.pp_snp_solution = test_utils.load_data_loc() + \
#                                'strainProfiler_v0.3_results/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted' \
#                                '.bam_SP_snpLocations.pickle '
#         BTO.cc_snp_solution = test_utils.load_data_loc() + \
#                                'v0.4_results/test_0.98.freq'
#         BTO.v12_solution = test_utils.load_data_loc() + \
#                             'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.2.14'
#         BTO.sam = test_utils.load_data_loc() + \
#                    'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam'
#         BTO.IS = test_utils.load_data_loc() + \
#                   'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
#         BTO.scafflist = test_utils.load_data_loc() + \
#                          'scaffList.txt'
#         BTO.genes = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_min1000.fa.genes.fna'
#         BTO.stb = test_utils.load_data_loc() + \
#                    'GenomeCoverages.stb'
#         BTO.genes = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_min1000.fa.genes.fna'
# 
#         twelve2thirteen, del_thirteen, new_thirteen = test_utils.get_twelve2thriteen()
#         BTO.twelve2thirteen = twelve2thirteen
#         BTO.del_thirteen = del_thirteen
#         BTO.new_thirteen = new_thirteen
# 
#     def setUp(BTO, destroy=True):
#         if destroy:
#             if os.path.isdir(BTO.test_dir):
#                 shutil.rmtree(BTO.test_dir)
#             os.mkdir(BTO.test_dir)
# 
#         importlib.reload(logging)
# 
#     def tearDown(BTO):
#         logging.shutdown()
#         if os.path.isdir(BTO.test_dir):
#             shutil.rmtree(BTO.test_dir)
# 
#     def run(BTO, min=0, max=19, tests='all', cleanUp=True):
#         # YOU HAVE TO RUN THIS ONE ON ITS OWN, BECUASE IT MESSES UP FUTURE RUNS
#         # BTO.setUp()
#         # BTO.test0()
#         # BTO.tearDown()
# 
#         if tests == 'all':
#             tests = np.arange(min, max+1)
# 
#         for test_num in tests:
#             BTO.setUp()
#             print("\n*** Running {1} test {0} ***\n".format(test_num, BTO.__class__))
#             eval('BTO.test{0}()'.format(test_num))
#             if cleanUp:
#                 BTO.tearDown()
# 
#         # BTO.setUp(destroy=True)
#         # BTO.test16()
#         # BTO.tearDown()

def test_profile_0(BTO):
    """
    Compare Matts to CCs methodology
    """
    # Run Matts program
    base = BTO.test_dir + 'testMatt'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --min_mapq 2 --skip_genome_wide".format(
        True, BTO.bam1,
        BTO.fasta, base)
    print(cmd)
    sys_args = cmd.split(' ')
    args = inStrain.argumentParser.parse_args(sys_args[1:])
    inStrain.controller.Controller().main(args)

    # Load the object
    Matt_object = inStrain.SNVprofile.SNVprofile(base)

    # Run CCs program
    base = BTO.test_dir + 'testCC'
    cmd = "{0} {1} {2} -o {3} -l 0.95".format(True, BTO.bam1,
                                              BTO.fasta, base)
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
    CPdb = CPdb.drop_duplicates(subset=['scaff', 'loc'])
    CPdb = CPdb.sort_values(['scaff', 'loc'])

    # Load Matts beautiful object
    MPdb = Matt_object.get_nonredundant_snv_table()

    # Allowing for cryptic SNPs, make sure Matt calls everything CC does
    MS = set(["{0}-{1}".format(x, y) for x, y in zip(MPdb['scaffold'], MPdb['position'])])
    CS = set(["{0}-{1}".format(x, y) for x, y in zip(CPdb['scaff'], CPdb['loc'])])
    assert len(MS) > 0
    assert len(CS - MS) == 0, CS - MS

    # Not allowing for cyptic SNPs, make sure CC calls everything Matt does
    MPdb = MPdb[MPdb['cryptic'] == False]
    MPdb = MPdb[MPdb['allele_count'] >= 2]
    MS = set(["{0}-{1}".format(x, y) for x, y in zip(MPdb['scaffold'], MPdb['position'])])
    CS = set(["{0}-{1}".format(x, y) for x, y in zip(CPdb['scaff'], CPdb['loc'])])
    assert len(MS - CS) == 0, MS - CS

    ## COMPARE CLONALITY

    # Parse CCs dumb clonality table
    CLdb = CC_object.clonality_table
    p2c = CLdb.set_index('position')['clonality'].to_dict()

    # Load Matt's beautiful table
    MCLdb = Matt_object.get_clonality_table()

    # print(set(p2c.keys()) - set(["{0}_{1}".format(s, p) for s, p in zip(MCLdb['scaffold'], MCLdb['position'])]))
    assert len(MCLdb) == len(CLdb), (len(MCLdb), len(CLdb))
    for i, row in MCLdb.dropna().iterrows():
        assert (p2c["{0}_{1}".format(row['scaffold'], row['position'])]
                - row['clonality']) < .001, (row, p2c["{0}_{1}".format(row['scaffold'], row['position'])])

    ## COMPARE LINKAGE
    CLdb = CC_object.r2linkage_table
    CLdb['position_A'] = [eval(str(x))[0].split(':')[0].split('_')[-1] for x in CLdb['total_A']]
    CLdb['position_B'] = [eval(str(x))[0].split(':')[0].split('_')[-1] for x in CLdb['total_B']]
    CLdb['scaffold'] = [x.split(':')[0] for x in CLdb['Window']]

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
    MS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in
              zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])])
    CS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in
              zip(CLdb['scaffold'], CLdb['position_A'], CLdb['position_B'])])
    assert len(CS - MS) <= 1, [CS - MS]
    # At scaffold N5_271_010G1_scaffold_110 from position 525 to 546 you end up in an edge case
    # where you skip it because you have absolutely no minor alleles to counterbalance it. It's fine,
    # CC just reports an r2 of np.nan, and this seems like the easiest way to handle it

    # Not allowing for cyptic SNPs, make sure CC calls everything Matt does
    MLdb = MLdb[MLdb['cryptic'] == False]
    MS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in
              zip(MLdb['scaffold'], MLdb['position_A'], MLdb['position_B'])])
    CS = set(["{0}-{1}-{2}".format(x, y, z) for x, y, z in
              zip(CLdb['scaffold'], CLdb['position_A'], CLdb['position_B'])])
    assert len(MS - CS) == 0, [len(MS), len(CS), len(MS - CS), MS - CS]

def test_profile_1(BTO):
    """
    Basic test- Make sure whole version doesn't crash when run from the command line
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 -d".format(True, BTO.bam1,
                                                           BTO.fasta_extra, base)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    assert os.path.isdir(base)
    assert len(glob.glob(base + '/output/*')) == 5, glob.glob(base + '/output/*')

    # Make sure the output makes sense
    S1 = inStrain.SNVprofile.SNVprofile(base)
    db = S1.get('cumulative_scaffold_table')
    assert len(db) > 0
    test_utils._internal_verify_Sdb(db)

    # Make sure it doesn't mess up at the lower-case bases (I put a lower-case in scaffold N5_271_010G1_scaffold_0 at a c; make sure it's not there)
    sdb = S1.get('raw_snp_table')
    assert len(sdb) > 0
    assert len(sdb[sdb['ref_base'] == 'c']) == 0

def test_profile_2(BTO):
    """
    Test filter reads; make sure CCs and Matt's agree
    """
    # Set up
    positions, total_length = inStrain.deprecated.deprecated_filter_reads.get_fasta(BTO.fasta)
    min_read_ani = 0.98

    # Run initial filter_reads
    subset_reads, Rdb = inStrain.deprecated.deprecated_filter_reads.filter_reads(BTO.bam1, positions,
                                                                                 total_length,
                                                                                 filter_cutoff=min_read_ani,
                                                                                 max_insert_relative=3,
                                                                                 min_insert=50, min_mapq=2)

    # Run Matts filter_reads
    scaff2sequence = SeqIO.to_dict(SeqIO.parse(BTO.fasta, "fasta"))  # set up .fasta file
    s2l = {s: len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}  # Get scaffold2length
    scaffolds = list(s2l.keys())
    subset_reads2 = inStrain.deprecated.deprecated_filter_reads.filter_paired_reads(BTO.bam1,
                                                                                    scaffolds,
                                                                                    filter_cutoff=min_read_ani,
                                                                                    max_insert_relative=3,
                                                                                    min_insert=50, min_mapq=2)

    # Run Matts filter_reads in a different way
    pair2info = inStrain.deprecated.deprecated_filter_reads.get_paired_reads(BTO.bam1, scaffolds)
    pair2infoF = inStrain.deprecated.deprecated_filter_reads.filter_paired_reads_dict(pair2info,
                                                                                      filter_cutoff=min_read_ani,
                                                                                      max_insert_relative=3,
                                                                                      min_insert=50, min_mapq=2)
    subset_readsF = list(pair2infoF.keys())

    # Run Matts filter_reads in a different way still
    scaff2pair2infoM, Rdb = inStrain.filter_reads.load_paired_reads(
        BTO.bam1, scaffolds, min_read_ani=min_read_ani,
        max_insert_relative=3, min_insert=50, min_mapq=2)
    # pair2infoMF = inStrain.filter_reads.paired_read_filter(scaff2pair2infoM)
    # pair2infoMF = inStrain.filter_reads.filter_scaff2pair2info(pair2infoMF,
    #                 min_read_ani=min_read_ani, max_insert_relative=3,
    #                 min_insert=50, min_mapq=2)

    subset_readsMF = set()
    for scaff, pair2infoC in scaff2pair2infoM.items():
        subset_readsMF = subset_readsMF.union(pair2infoC.keys())
    # subset_readsMF = list(pair2infoMF.keys())

    assert (set(subset_reads2) == set(subset_reads) == set(subset_readsF) == set(subset_readsMF)), \
        [len(subset_reads2), len(subset_reads), len(subset_readsF), len(subset_readsMF)]

    # Make sure the filter report is accurate
    # Rdb = inStrain.filter_reads.makeFilterReport2(scaff2pair2infoM, pairTOinfo=pair2infoMF, min_read_ani=min_read_ani, max_insert_relative=3,
    #                                     min_insert=50, min_mapq=2)
    assert int(Rdb[Rdb['scaffold'] == 'all_scaffolds'] \
                   ['unfiltered_pairs'].tolist()[0]) \
           == len(list(pair2info.keys()))
    assert int(Rdb[Rdb['scaffold'] == 'all_scaffolds']['filtered_pairs'].tolist()[0]) \
           == len(subset_reads)

    # Try another cutuff
    positions, total_length = inStrain.deprecated.deprecated_filter_reads.get_fasta(BTO.fasta)
    min_read_ani = 0.90

    # Run initial filter_reads
    subset_reads, Rdb = inStrain.deprecated.deprecated_filter_reads.filter_reads(BTO.bam1, positions,
                                                                                 total_length,
                                                                                 filter_cutoff=min_read_ani,
                                                                                 max_insert_relative=3,
                                                                                 min_insert=50, min_mapq=2)

    # Run Matts filter_reads
    scaff2sequence = SeqIO.to_dict(SeqIO.parse(BTO.fasta, "fasta"))  # set up .fasta file
    s2l = {s: len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}  # Get scaffold2length
    scaffolds = list(s2l.keys())

    scaff2pair2infoM, Rdb = inStrain.filter_reads.load_paired_reads(
        BTO.bam1, scaffolds, min_read_ani=min_read_ani,
        max_insert_relative=3, min_insert=50, min_mapq=2)
    pair2infoMF_keys = set()
    for scaff, pair2infoC in scaff2pair2infoM.items():
        pair2infoMF_keys = pair2infoMF_keys.union(pair2infoC.keys())
    # Scaff2pair2infoM = inStrain.filter_reads.get_paired_reads_multi(BTO.bam1, scaffolds)
    # pair2infoMF = inStrain.filter_reads.paired_read_filter(scaff2pair2infoM)
    # pair2infoMF = inStrain.filter_reads.filter_paired_reads_dict2(pair2infoMF,
    #                 min_read_ani=min_read_ani, max_insert_relative=3,
    #                 min_insert=50, min_mapq=2)

    subset_reads2 = pair2infoMF_keys

    assert (set(subset_reads2) == set(subset_reads))

def test_profile_3(BTO):
    """
    Testing scaffold table (breadth and coverage) vs. calculate_breadth
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --skip_genome_wide".format(True, BTO.bam1,
                                                                              BTO.fasta, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    Sprofile = inStrain.SNVprofile.SNVprofile(base)
    Odb = Sprofile.get('cumulative_scaffold_table')

    # Verify coverage table
    test_utils._internal_verify_Sdb(Odb)

    # Compare to calculate_coverage
    Cdb = pd.read_csv(BTO.cc_solution)
    s2c = Cdb.set_index('scaffold')['coverage'].to_dict()
    s2b = Cdb.set_index('scaffold')['breadth'].to_dict()
    for scaff, db in Odb.groupby('scaffold'):
        db = db.sort_values('mm', ascending=False)
        assert (db['coverage'].tolist()[0] - s2c[scaff]) < .1, \
            [db['coverage'].tolist()[0], s2c[scaff]]
        assert (db['breadth'].tolist()[0] - s2b[scaff]) < .01, \
            [db['breadth'].tolist()[0], s2b[scaff]]

    # Verify SNP calls
    Sdb = Sprofile.get('cumulative_snv_table')
    test_utils._internal_verify_OdbSdb(Odb, Sdb)

def test_profile_4(BTO):
    """
    Test store_everything and database mode
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --skip_plot_generation -s {4}".format(
        True, BTO.bam1,
        BTO.fasta, base, BTO.stb)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    Sprofile = inStrain.SNVprofile.SNVprofile(base)

    # Make sure you stored a heavy object
    assert Sprofile.get('testeroni') is None
    assert Sprofile.get('covT') is not None
    assert Sprofile.get('mm_to_position_graph') is not None

    # Run database mode
    base2 = BTO.test_dir + 'test2'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --database_mode --skip_plot_generation -s {4}".format(
        True, BTO.bam1,
        BTO.fasta, base2, BTO.stb)
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
    assert not test_utils.compare_dfs2(mdb1, mdb2, verbose=True)

    # Make sure you have more skip mm level
    mdb1 = Sprofile.get('genome_level_info')
    mdb2 = Sprofile2.get('genome_level_info')
    assert 'mm' in mdb1
    assert 'mm' not in mdb2

    # Make sure you skip junk genomes
    # assert len(set(mdb1['genome']) - set(mdb2['genome'])) > 0

def test_profile_5(BTO):
    """
    Test one thread
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 1 --skip_genome_wide".format(True, BTO.bam1,
                                                                                   BTO.fasta, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    Sprofile = inStrain.SNVprofile.SNVprofile(base)

    # Make sure you get a result
    Odb = Sprofile.get('cumulative_scaffold_table')
    assert len(Odb['scaffold'].unique()) == 178, Odb

def test_profile_6(BTO):
    """
    Test the case where only one scaffold is  preset at all in the .bam file AND it has no SNPs
    """
    # Run program
    base = BTO.test_dir + 'test'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.99 -p 1 --skip_genome_wide".format(True, BTO.bam1,
                                                                                   BTO.single_scaff, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    Sprofile = inStrain.SNVprofile.SNVprofile(base)

    # Load output
    Odb = Sprofile.get('cumulative_scaffold_table')
    print(Odb)
    test_utils._internal_verify_Sdb(Odb)

def test_profile_7(BTO):
    """
    Test the case where a scaffold is not preset at all in the .bam file

    Also test being able to adjust read filtering parameters
    """
    # Run program
    base = BTO.test_dir + 'test'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.80 -p 6 --store_everything --skip_genome_wide --skip_plot_generation".format(
        True, BTO.bam1,
        BTO.extra_single_scaff, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    Sprofile = inStrain.SNVprofile.SNVprofile(base)

    # Make sure scaffold table is OK
    Odb = Sprofile.get('cumulative_scaffold_table')
    test_utils._internal_verify_Sdb(Odb)

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

        cmd = "inStrain profile {1} {2} -o {3} --{4} {5} -p 6 --store_everything --skip_plot_generation".format(
            True, BTO.bam1,
            BTO.extra_single_scaff, base, thing, val)
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
        # assert rdb[rdb['scaffold'] == 'all_scaffolds']["pass_{0}".format(thing)].tolist()[0] == 0

        reads = len(Sprofile.get('Rdic').keys())
        print(Sprofile.get('Rdic'))
        assert reads < ORI_READS

def test_profile_8(BTO):
    """
    Test the ability to make and sort .bam files from .sam
    """
    importlib.reload(logging)
    # Copy sam to test dir
    new_sam = os.path.join(BTO.test_dir, os.path.basename(BTO.sam))
    shutil.copyfile(BTO.sam, new_sam)

    # Run program
    base = BTO.test_dir + 'test'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.80 -p 6 --store_everything --skip_genome_wide -d".format(
        True, new_sam,
        BTO.extra_single_scaff, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

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

def test_profile_9(BTO):
    """
    Test the debug option

    v0.5.1 - Actually this should happen all the time now...
    v1.2.0 - This test is obsolete now
    """
    pass
    # # Run Matts program
    # base = BTO.test_dir + 'testMatt'
    # # cmd = "{0} {1} {2} -o {3} -l 0.95 --store_everything --debug".format(True, BTO.bam1, \
    # #     BTO.fasta, base)
    # cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --store_everything --skip_genome_wide".format(True, BTO.bam1, \
    #     BTO.fasta, base)
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

def test_profile_10(BTO):
    """
    Test min number of reads filtered and min number of genome coverage
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_scaffold_reads 10 --skip_genome_wide".format(True,
                                                                                                      BTO.bam1,
                                                                                                      BTO.fasta,
                                                                                                      base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    # Make sure you actually filtered out the scaffolds
    Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
    Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
    print("{0} of {1} scaffolds have >10 reads".format(len(Rdb[Rdb['filtered_pairs'] >= 10]),
                                                       len(Rdb)))
    assert len(Sdb['scaffold'].unique()) == len(Rdb[Rdb['filtered_pairs'] >= 10]['scaffold'].unique()) - 1

    # Try min_genome_coverage

    # Make sure it fails with no .stb
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 5 --skip_genome_wide".format(True,
                                                                                                      BTO.bam1,
                                                                                                      BTO.fasta,
                                                                                                      base)
    exit_code = call(cmd, shell=True)
    assert exit_code == 1

    # Run it with an .stb
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 -s {4} --min_genome_coverage 5 --skip_genome_wide".format(True,
                                                                                                      BTO.bam1,
                                                                                                      BTO.fasta,
                                                                                                      base,
                                                                                                      BTO.stb2)
    exit_code = call(cmd, shell=True)
    assert exit_code == 0

    # Make sure you actually filtered out the scaffolds
    Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
    Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
    assert len(Rdb) == 179
    print(Sdb)
    assert len(Sdb) == 42

    # Make sure empty scaffolds don't mess it up
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 5  -s {4} --skip_genome_wide".format(
        True, BTO.bam1, BTO.fasta_extra, base, BTO.stb2)
    exit_code = call(cmd, shell=True)
    assert exit_code == 0, exit_code

    # Make sure you actually filtered out the scaffolds
    Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
    Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
    assert len(Rdb) == 180, len(Rdb)
    assert len(Sdb) == 42, len(Sdb)

    # Run it with an .stb and coverage that cannot be hit
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 100  -s {4} --skip_genome_wide".format(
        True, BTO.bam1, BTO.fasta, base, BTO.stb2)
    exit_code = call(cmd, shell=True)
    assert exit_code == 1

    # Run it with an .stb and coverage that is low
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 --min_genome_coverage 1.1  -s {4} --skip_genome_wide".format(
        True, BTO.bam1, BTO.fasta, base, BTO.stb2)
    exit_code = call(cmd, shell=True)
    Sprofile = inStrain.SNVprofile.SNVprofile(base)

    # Make sure you actually filtered out the scaffolds
    Sdb = pd.read_csv(glob.glob(base + '/output/*scaffold_info.tsv')[0], sep='\t')
    Rdb = pd.read_csv(glob.glob(base + '/output/*mapping_info.tsv')[0], sep='\t', header=1)
    assert len(Rdb) == 179
    assert len(Sdb) > 42

def test_profile_11(BTO):
    """
    Test skip mm profiling

    # NOTE! THIS TEST FAILS WITH NEW PYSAM (ONLY WORKS WITH PYSAM=0.15.4). NEED TO FIX
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} --skip_mm_profiling --skip_genome_wide".format(True,
                                                                                          BTO.bam1,
                                                                                          BTO.fasta, base)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    IS = inStrain.SNVprofile.SNVprofile(base)

    # Make sure you get the same results
    scaffdb = IS.get_nonredundant_scaffold_table().reset_index(drop=True)
    correct_scaffdb = inStrain.SNVprofile.SNVprofile(BTO.IS).get_nonredundant_scaffold_table().reset_index(
        drop=True)

    cols = ['scaffold', 'length', 'breadth', 'coverage']
    # assert test_utils.compare_dfs2(scaffdb[cols].sort_values(['scaffold']).reset_index(drop=True),
    #                               correct_scaffdb[cols].sort_values(['scaffold']).reset_index(drop=True))
    assert test_utils.compare_dfs(scaffdb[cols].sort_values(['scaffold']).reset_index(drop=True),
                                  correct_scaffdb[cols].sort_values(['scaffold']).reset_index(drop=True),
                                  verbose=True, round=4)

    # Make sure you dont have the raw mm
    sdb = IS.get('cumulative_scaffold_table')
    print(sdb.head())
    assert set(sdb['mm'].tolist()) == {0}, set(sdb['mm'].tolist())

def test_profile_12(BTO):
    """
    Test scaffolds_to_profile
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --min_mapq 2 --scaffolds_to_profile {4} --skip_genome_wide".format(
        True, BTO.bam1,
        BTO.fasta, base, BTO.scafflist)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    IS = inStrain.SNVprofile.SNVprofile(base)

    # Make sure you get the same results
    scaffdb = IS.get_nonredundant_scaffold_table().reset_index(drop=True)
    assert set(scaffdb['scaffold'].unique()) == set(inStrain.profile.fasta.load_scaff_list(BTO.scafflist))

def test_profile_13(BTO):
    """
    Make sure that covT, clonT, and the SNP table agree on coverage
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 --min_mapq 2 --scaffolds_to_profile {4} -p 1 --skip_genome_wide".format(
        True, BTO.bam1,
        BTO.fasta, base, BTO.scafflist)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))
    IS = inStrain.SNVprofile.SNVprofile(base)

    SRdb = IS.get('cumulative_snv_table')
    CovT = IS.get('covT')

    for i, row in SRdb.iterrows():
        cov = 0
        for mm, covs in CovT[row['scaffold']].items():
            if mm <= row['mm']:
                if row['position'] in covs:
                    cov += covs[row['position']]
        assert row['position_coverage'] == cov, [cov, row['position_coverage'], row]

def test_profile_14(BTO):
    """
    Basic test- Make sure genes and genome_wide can be run within the profile option

    Make sure logging produces the run statistics
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98".format(True, BTO.bam1,
                                                                  BTO.fasta, base, BTO.genes)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    assert os.path.isdir(base)
    assert len(glob.glob(base + '/output/*')) == 6, len(glob.glob(base + '/output/*'))
    assert len(glob.glob(base + '/log/*')) == 3

    # Make sure the output makes sense
    S1 = inStrain.SNVprofile.SNVprofile(base)
    db = S1.get('cumulative_scaffold_table')
    test_utils._internal_verify_Sdb(db)

    # Read the log
    for l in glob.glob(base + '/log/*'):
        if 'runtime_summary.txt' in l:
            with open(l, 'r') as o:
                for line in o.readlines():
                    print(line.strip())

def test_profile_15(BTO):
    """
    Basic test- Make sure genes and genome_wide can be run within the profile option
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98 --rarefied_coverage 10".format(True,
                                                                                         BTO.bam1,
                                                                                         BTO.fasta, base,
                                                                                         BTO.genes)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    assert os.path.isdir(base)
    assert len(glob.glob(base + '/output/*')) == 6

    # Make sure the output makes sense
    S1 = inStrain.SNVprofile.SNVprofile(base)
    db = S1.get('cumulative_scaffold_table')
    test_utils._internal_verify_Sdb(db)
    clontR = S1.get('clonTR')
    counts0 = sum([len(x[2]) if 2 in x else 0 for s, x in clontR.items()])

    # Make sure its in the genome_wide table
    gdb = pd.read_csv(glob.glob(base + '/output/*genome_info*.tsv')[0], sep='\t')
    assert 'nucl_diversity' in gdb.columns, gdb.head()

    # Run again with different rarefied coverage
    base = BTO.test_dir + 'test2'
    cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98 --rarefied_coverage 50".format(True,
                                                                                         BTO.bam1,
                                                                                         BTO.fasta, base,
                                                                                         BTO.genes)
    print(cmd)
    call(cmd, shell=True)

    S1 = inStrain.SNVprofile.SNVprofile(base)
    db = S1.get('cumulative_scaffold_table')
    test_utils._internal_verify_Sdb(db)
    clontR = S1.get('clonTR')
    counts2 = sum([len(x[2]) if 2 in x else 0 for s, x in clontR.items()])

    assert counts0 > counts2, [counts0, counts2]

def test_profile_16(BTO):
    """
    Make sure the results exactly match a run done with inStrain verion 1.2.14
    """
    importlib.reload(logging)
    # Run program
    base = BTO.test_dir + 'test'
    cmd = "inStrain profile {1} {2} -o {3} -g {4} --skip_plot_generation -p 6 -d".format(True,
                                                                                         BTO.bam1,
                                                                                         BTO.fasta, base,
                                                                                         BTO.genes)
    print(cmd)
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split(' ')[1:]))

    exp_IS = inStrain.SNVprofile.SNVprofile(base)
    sol_IS = inStrain.SNVprofile.SNVprofile(BTO.v12_solution)

    # Print what the output of the solutions directory looks like
    if True:
        s_out_files = glob.glob(exp_IS.get_location('output') + os.path.basename(
            exp_IS.get('location')) + '_*')
        print("The output has {0} tables".format(len(s_out_files)))
        for f in s_out_files:
            name = os.path.basename(f)
            print("{1}\n{0}\n{1}".format(name, '-' * len(name)))

            if 'mapping_info.tsv' in name:
                s = pd.read_csv(f, sep='\t', header=1)
            else:
                s = pd.read_csv(f, sep='\t')
            print(s.head())
            print()

    # MAKE SURE LOG IS WORKING
    assert len(glob.glob(base + '/log/*')) == 3, base
    Ldb = exp_IS.get_parsed_log()
    rdb = inStrain.logUtils.load_multiprocessing_log(Ldb)

    LOGGED_SCAFFOLDS = set(rdb[(rdb['command'] == 'MergeProfile') & (rdb['multi_log_type'] == 'WorkerLog')]['unit'].tolist())
    TRUE_SCAFFOLDS = \
        set(exp_IS.get_nonredundant_scaffold_table()['scaffold'].tolist())
    assert (LOGGED_SCAFFOLDS == TRUE_SCAFFOLDS)

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
            # print("Both have {0}!".format(name))

            e = pd.read_csv(e_file[0], sep='\t')
            s = pd.read_csv(s_file, sep='\t').rename(columns=BTO.twelve2thirteen)

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

                for col in ['mutation_type', 'mutation', 'gene', 'class']:
                    del e[col]

            if name in ['scaffold_info.tsv']:
                e = e.sort_values(['scaffold']).reset_index(drop=True)
                s = s.sort_values(['scaffold']).reset_index(drop=True)

                # TRANSLATE THE OLD VERSION
                s = s.rename(columns=BTO.twelve2thirteen)
                for r in BTO.del_thirteen:
                    if r in s.columns:
                        del s[r]
                for r in BTO.new_thirteen:
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
                s = s.rename(columns=BTO.twelve2thirteen)
                for r in BTO.del_thirteen:
                    if r in s.columns:
                        del s[r]
                for r in BTO.new_thirteen:
                    if r in e.columns:
                        del e[r]

                rand = ['SNS_count', 'divergent_site_count', 'partial']
                for r in rand:
                    if r in e.columns:
                        del e[r]
                    if r in s.columns:
                        del s[r]

            # Re-arange column order
            assert set(e.columns) == set(s.columns), \
                [name,
                 set(e.columns) - set(s.columns),
                 set(s.columns) - set(e.columns)]
            s = s[list(e.columns)]
            e = e[list(e.columns)]

            assert test_utils.compare_dfs2(e, s, verbose=True), name

        else:
            # print("Both dont have {0}!".format(name))
            if name in ['read_report.tsv']:
                e_file = [e for e in e_out_files if 'mapping_info.tsv' in os.path.basename(e)]

                e = pd.read_csv(e_file[0], sep='\t', header=1)
                s = pd.read_csv(s_file, sep='\t', header=1).rename(
                    columns={'pass_filter_cutoff': 'pass_min_read_ani'})

                e = e.sort_values(['scaffold']).reset_index(drop=True)
                s = s.sort_values(['scaffold']).reset_index(drop=True)

                e = e[e['scaffold'] == 'all_scaffolds']
                s = s[s['scaffold'] == 'all_scaffolds']

                for c in list(s.columns):
                    s[c] = s[c].astype(e[c].dtype)

                for c in ['median_insert']:  # calculated in a different way
                    del e[c]
                    del s[c]

            elif name in ['genomeWide_scaffold_info.tsv']:
                e_file = [e for e in e_out_files if 'genome_info.tsv' in os.path.basename(e)]

                e = pd.read_csv(e_file[0], sep='\t')
                s = pd.read_csv(s_file, sep='\t').rename(columns=BTO.twelve2thirteen)

                for r in BTO.del_thirteen:
                    if r in s.columns:
                        del s[r]
                s = s.rename(columns=BTO.twelve2thirteen)

                NEW_HERE = {'coverage_median', 'SNV_count', 'SNS_count',
                            'nucl_diversity', 'filtered_read_pair_count'}
                for c in BTO.new_thirteen.union(NEW_HERE):
                    if c in e.columns:
                        del e[c]

                # Remove the ones that are gained by the new read filtering
                for c in e.columns:
                    if c.startswith('reads_'):
                        del e[c]

                e = e.sort_values(['genome']).reset_index(drop=True)
                s = s.sort_values(['genome']).reset_index(drop=True)

                # Re-order columns
                assert set(e.columns) == set(s.columns), \
                    [set(e.columns) - set(s.columns),
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
                s = pd.read_csv(s_file, sep='\t').rename(columns=BTO.twelve2thirteen) \
                    .rename(columns={'pass_filter_cutoff': 'pass_min_read_ani'})

                # TRANSLATE THE OLD VERSION
                s = s.rename(columns=BTO.twelve2thirteen)
                for r in BTO.del_thirteen:
                    if r in s.columns:
                        del s[r]
                for r in BTO.new_thirteen:
                    if r in e.columns:
                        del e[r]

                new_cols = list(set(e.columns) - set(s.columns))
                for c in new_cols:
                    del e[c]

                removed_cols = ['unfiltered_reads', 'pass_pairing_filter', 'pass_min_mapq', 'unfiltered_singletons',
                                'filtered_priority_reads', 'mean_insert_distance', 'mean_pair_length',
                                'pass_min_insert', 'pass_max_insert', 'pass_min_read_ani', 'mean_PID',
                                'mean_mistmaches', 'median_insert', 'mean_mapq_score', 'unfiltered_pairs',
                                'filtered_singletons', 'unfiltered_priority_reads', 'filtered_pairs', 'scaffolds']
                for r in removed_cols:
                    if r in s.columns:
                        del s[r]

                e = e.sort_values(['genome']).reset_index(drop=True)
                s = s.sort_values(['genome']).reset_index(drop=True)
                assert set(e.columns) == set(s.columns), [set(s.columns) - set(e.columns)]

            elif name in ['SNP_mutation_types.tsv']:
                e_file = [e for e in e_out_files if 'SNVs.tsv' in os.path.basename(e)]

                e = pd.read_csv(e_file[0], sep='\t')
                s = pd.read_csv(s_file, sep='\t').rename(columns=BTO.twelve2thirteen)

                e = e[~e['mutation_type'].isna()]
                del e['cryptic']
                del e['class']

                e = e.sort_values(['scaffold', 'position']).reset_index(drop=True)
                s = s.sort_values(['scaffold', 'position']).reset_index(drop=True)

            else:
                assert False, name

            # Re-arange column order
            assert set(e.columns) == set(s.columns), \
                [name,
                 set(e.columns) - set(s.columns),
                 set(s.columns) - set(e.columns)]

            s = s[list(e.columns)]

            assert test_utils.compare_dfs2(e, s, verbose=True), name

    # CHECK ATTRIBUTES
    sAdb = sol_IS._get_attributes_file()

    # Handle name changes
    o2n = {'genes_SNP_density': 'genes_SNP_count', 'read_report': 'mapping_info'}

    for i, row in sAdb.iterrows():
        print("checking {0}".format(i))

        if i in ['location', 'version', 'bam_loc', 'genes_fileloc', 'window_table']:
            continue

        s = sol_IS.get(i)
        if i in o2n:    i = o2n[i]
        e = exp_IS.get(i)

        if i in ['scaffold_list']:
            if not test_utils.compare_lists(e, s, check_order=False):
                print("{0} is not the same".format(i))
                print(e)
                print(s)
                assert e == s, i

        elif i in ['scaffold2length', 'scaffold2bin', 'bin2length']:
            assert test_utils.compare_dicts(e, s), i

        elif i in ['window_table', 'raw_linkage_table', 'raw_snp_table', 'cumulative_scaffold_table',
                   'cumulative_snv_table', 'mapping_info', 'genes_table', 'genes_coverage',
                   'genes_clonality', 'genes_SNP_count', 'SNP_mutation_types']:

            # TRANSLATE THE OLD VERSION
            s = s.rename(columns=BTO.twelve2thirteen).rename(
                columns={'pass_filter_cutoff': 'pass_min_read_ani'})
            for r in BTO.del_thirteen:
                if r in s.columns:
                    del s[r]
            for r in BTO.new_thirteen:
                if r in e.columns:
                    del e[r]

            if i in ['window_table', 'mapping_info']:
                e = e.sort_values(['scaffold']).reset_index(drop=True)
                s = s.sort_values(['scaffold']).reset_index(drop=True)
                e = e.rename(columns={'filtered_pairs': 'filtered_read_pair_count'})

            if i in ['mapping_info']:
                e = e[e['scaffold'] == 'all_scaffolds']
                s = s[s['scaffold'] == 'all_scaffolds']

                for c in list(s.columns):
                    s[c] = s[c].astype(e[c].dtype)

                for c in ['median_insert']:
                    del e[c]
                    del s[c]

            if i in ['raw_linkage_table']:
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
                e = e.sort_values(['scaffold', 'gene']).reset_index(drop=True)
                s = s.sort_values(['scaffold', 'gene']).reset_index(drop=True)

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
            assert set(e.columns) == set(s.columns), \
                [i,
                 set(e.columns) - set(s.columns),
                 set(s.columns) - set(e.columns)]

            s = s[list(e.columns)]

            assert test_utils.compare_dfs2(e, s, verbose=True), i

        elif i in ['scaffold_2_mm_2_read_2_snvs', 'clonTR']:
            pass

        elif i in ['covT', 'clonT']:
            assert test_utils.compare_covTs(e, s), "{0} is not the same".format(i)

        else:
            print("YOUR NOT CHECKING {0}".format(i))
            print(s)

def test_profile_17(BTO):
    """
    Test scaffold failure
    """
    importlib.reload(logging)

    # Set up
    base = BTO.test_dir + 'test'

    # Run program and make the split crash
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=3000 -d".format(
        True, BTO.failure_bam,
        BTO.failure_fasta, base)
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
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=20000 -d".format(
        True, BTO.failure_bam,
        BTO.failure_fasta, base)
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
    cmd = "inStrain profile {1} {2} -o {3} -l 0.95 -p 6 --skip_genome_wide --window_length=20000 -d -g {4}".format(
        True, BTO.failure_bam,
        BTO.failure_fasta, base, BTO.failure_genes)
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

def test_profile_18(BTO):
    """
    Test providing a .fasta file with a really small scaffold
    """
    base = BTO.test_dir + 'testR'

    cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation --debug".format(BTO.small_bam, BTO.small_fasta,
                                                                                  base)
    print(cmd)
    call(cmd, shell=True)

    IS = inStrain.SNVprofile.SNVprofile(base)

    db = IS.get('genome_level_info')
    assert db is not None
    assert len(glob.glob(base + '/log/*')) == 3

def test_profile_19(BTO):
    """
    Test the case where a scaffold is not in the .bam but is in the .stb
    """
    # Run program
    base = BTO.test_dir + 'test'
    cmd = "inStrain profile {1} {2} -o {3} -l 0.80 -p 6 -s {4} --store_everything --skip_plot_generation".format(
        True, BTO.bam1,
        BTO.extra_single_scaff, base,
        BTO.extra_single_scaff.replace('.fasta', '.stb'))
    print(cmd)
    sys_args = cmd.split(' ')
    args = inStrain.argumentParser.parse_args(sys_args[1:])
    inStrain.controller.Controller().main(args)

    # Make sure scaffold table is OK
    floc = os.path.join(base, 'output', 'test_genome_info.tsv')
    db = pd.read_csv(floc, sep='\t')
    assert db['true_scaffolds'].iloc[0] == 2
    assert db['length'].iloc[0] == 1348

def test_profile_20(BTO):
    """
    Test the "force compress" option
    """
    # Set up
    base = BTO.test_dir + 'test'

    # Run program
    cmd = "inStrain profile {1} {2} -o {3} -g {4} -l 0.98 --force_compress".format(True, BTO.bam1,
                                                                  BTO.fasta, base, BTO.genes)
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    assert os.path.isdir(base)
    print(glob.glob(base + '/output/*'))
    assert len(glob.glob(base + '/output/*.gz')) == 6, len(glob.glob(base + '/output/*.gz'))
    assert len(glob.glob(base + '/log/*')) == 3

    # Make sure the output makes sense
    S1 = inStrain.SNVprofile.SNVprofile(base)
    db = S1.get('cumulative_scaffold_table')
    test_utils._internal_verify_Sdb(db)

# This is from test_compare.py
def test_UPDATE_COMPARE_TEST_DATA(BTO):
    """
    Run inStrain on bam1 and bam2, and store the results where IS1 and IS2 are
    """

    # Run the first one
    base1 = BTO.test_dir + os.path.basename(BTO.IS1)
    cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(BTO.bam1, BTO.fasta,
                                                                          base1)
    print(cmd)
    code = call(cmd, shell=True)
    assert code == 0, code

    # Copy to new location
    if os.path.isdir(BTO.IS1):
        shutil.rmtree(BTO.IS1)
    shutil.copytree(base1, BTO.IS1)
    #
    # # Run the second one
    # base2 = BTO.test_dir + os.path.basename(BTO.IS2)
    # cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(BTO.bam2, BTO.fasta,
    #                                                                       base2)
    # print(cmd)
    # code = call(cmd, shell=True)
    # assert code == 0, code
    #
    # # Copy to new location
    # if os.path.isdir(BTO.IS2):
    #     shutil.rmtree(BTO.IS2)
    # shutil.copytree(base2, BTO.IS2)
    #
