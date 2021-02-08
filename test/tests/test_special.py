"""
Run tests on inStrain plotting
"""

import glob
import os
import shutil
from subprocess import call, run, PIPE
import numpy as np

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.irep_utilities
import inStrain.profile.fasta
import tests.test_utils as test_utils
import inStrain.logUtils

from test_utils import BTO

# class test_special:
#     def __init__(BTO):
#         BTO.test_dir = test_utils.load_random_test_dir()
#         BTO.bam1 = test_utils.load_data_loc() + \
#                           'N5_271_010G1_scaffold_963_Ns.fasta.sorted.bam'
#         BTO.fasta = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_963_Ns.fasta'
#         BTO.bbsam = test_utils.load_data_loc() + \
#                      'bbmap_N5_271_010G1_scaffold_963.fasta.sorted.bam'
#         BTO.bbfasta = test_utils.load_data_loc() + \
#                        'N5_271_010G1_scaffold_963.fasta'
#         BTO.IS = test_utils.load_data_loc() + \
#                   'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
#         BTO.small_fasta = test_utils.load_data_loc() + \
#                            'SmallScaffold.fa'
#         BTO.small_bam = test_utils.load_data_loc() + \
#                          'SmallScaffold.fa.sorted.bam'
#         BTO.genes = test_utils.load_data_loc() + \
#                      'N5_271_010G1_scaffold_min1000.fa.genes.fna'
#         BTO.stb = test_utils.load_data_loc() + \
#                    'GenomeCoverages.stb'
#         BTO.IS1 = test_utils.load_data_loc() + \
#                    'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.forRC.IS'
# 
#     def setUp(BTO):
#         BTO.tearDown()
# 
#     def tearDown(BTO):
#         if os.path.isdir(BTO.test_dir):
#             shutil.rmtree(BTO.test_dir)
# 
#     def run(BTO, min=1, max=5, tests='all'):
# 
#         if tests == 'all':
#             tests = np.arange(min, max+1)
# 
#         for test_num in tests:
#             BTO.setUp()
#             print("\n*** Running {1} test {0} ***\n".format(test_num, BTO.__class__))
#             eval('BTO.test{0}()'.format(test_num))
#             BTO.tearDown()

def test_special_1(BTO):
    """
    Yell if the minor version is not correct
    """
    assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.3')

    assert not inStrain.SNVprofile.same_versions('1.0.0', '0.8.3')

    assert inStrain.SNVprofile.same_versions('0.8.3', '0.8.4')

    assert not inStrain.SNVprofile.same_versions('0.9.0', '0.8.4')

    assert not inStrain.SNVprofile.same_versions('1.1.0', '0.1.0')

def test_special_2(BTO):
    """
    Make sure it works with Ns in the input

    Also make sure it logs dependencies
    """
    # Run base
    base = BTO.test_dir + 'testNs'
    cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(
        BTO.bam1,
        BTO.fasta, base)
    print(cmd)
    call(cmd, shell=True)

    # Load the log and make sure dependency report is there
    log_loc = os.path.join(base, 'log/log.log')
    got = False
    with open(log_loc, 'r') as o:
        for line in o.readlines():
            if "$ DEPENDENCY REPORT $" in line:
                got = True
    assert got

def test_special_3(BTO):
    """
    Try other mappers
    """
    # Run bbmap and let it crash
    base = BTO.test_dir + 'bbMap'
    cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation".format(BTO.bbsam,
                                                                                                     BTO.bbfasta,
                                                                                                     base)
    print(cmd)
    call(cmd, shell=True)

    # Load the object
    Matt_object = inStrain.SNVprofile.SNVprofile(base)
    MCLdb = Matt_object.get_nonredundant_scaffold_table()
    assert len(MCLdb) == 0

    # Run bbmap and make it work
    base = BTO.test_dir + 'bbMap'
    cmd = "inStrain profile {0} {1} -o {2} -l 0.95 --skip_genome_wide --skip_plot_generation --use_full_fasta_header".format(
        BTO.bbsam,
        BTO.bbfasta, base)
    print(cmd)
    call(cmd, shell=True)

    # Load the object
    Matt_object = inStrain.SNVprofile.SNVprofile(base)
    MCLdb = Matt_object.get_nonredundant_scaffold_table()
    assert len(MCLdb) > 0

def test_special_4(BTO):
    """
    Test other --run_statistics
    """
    # Super quick profile run
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS))
    cmd = "inStrain profile {0} {1} -o {2} -s {3} -g {4} --skip_plot_generation".format(
        BTO.small_bam, BTO.small_fasta, location, BTO.stb, BTO.genes)
    print(cmd)
    call(cmd, shell=True)

    cmd = "inStrain other --run_statistics {0} --debug".format(location)
    print(cmd)
    call(cmd, shell=True)
    assert len(glob.glob(location + '/log/*')) == 3, location

def test_special_5(BTO):
    """
    Test --version
    """
    # Super quick profile run
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS))
    cmd = "inStrain profile {0} {1} -o {2} -s {3} -g {4} --version".format(
        BTO.small_bam, BTO.small_fasta, location, BTO.stb, BTO.genes)
    print(cmd)
    result = run(cmd, shell=True, stdout=PIPE)

    # Make sure the code is right
    assert result.returncode == 0

    # Make sure the program wasn't run
    assert not os.path.isdir(location)

    # Make sure the version was printed to STDOUT
    m = result.stdout.decode("utf-8").strip()
    assert 'inStrain version' in m

def test_special_6(BTO):
    """
    Some specific tests for logging
    """
    # Copy over
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS1))
    shutil.copytree(BTO.IS1, location)

    cmd = "inStrain other --run_statistics {0} --debug".format(location)
    print(cmd)
    call(cmd, shell=True)
    logs = glob.glob(location + '/log/*')
    assert len(logs) == 4, logs


    # cmd = "inStrain profile {0} {1} -o {2} -s {3} -g {4}".format(
    #     BTO.small_bam, BTO.small_fasta, location, BTO.stb, BTO.genes)
    # print(cmd)
    # result = run(cmd, shell=True, stdout=PIPE)
    #
    # # Make sure the code is right
    # assert result.returncode == 0

    # Load a log file
    # loc = glob.glob(location + '/log/*runtime*')[0]
    # with open(loc) as l:
    #     for line in l.readlines():
    #         print(line.strip())
