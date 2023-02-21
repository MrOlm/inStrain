"""
Run tests on inStrain genome_wide
"""

import glob
import os
import pickle
import shutil
from subprocess import call

import numpy as np
import pandas as pd

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.irep_utilities
import inStrain.profile.fasta
import tests.test_utils as test_utils

from test_utils import BTO

def test_genome_level_0(BTO):
    """
    Just run the program
    """
    # Run the dereplicated version
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_nogenes))
    shutil.copytree(BTO.IS_nogenes, location)
    gloc = location + '/output/test_genome_info.tsv'
    print(gloc)
    os.remove(gloc)

    cmd = "inStrain genome_wide -i {0} -s {1} -d".format(location, BTO.stb)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    files = glob.glob(IS.get_location('output') + '*')
    print(files)
    files = [f for f in files if 'genome_info' in f]
    assert len(files) == 1, files

    # Make sure mm level is False in the output
    for f in files:
        db = pd.read_csv(f, sep='\t')
        assert 'mm' not in db.columns
        print(db.head())

    # Make sure mm level is True in the raw data
    db = IS.get('genome_level_info')
    assert 'mm' in db.columns
    # print(IS)

def test_genome_level_1(BTO):
    """
    Test the different ways of importing the .stb file
    """
    # Setup
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_nogenes))
    shutil.copytree(BTO.IS_nogenes, location)
    for f in glob.glob(location + '/output/*.tsv'):
        os.remove(f)

    # Run with a normal .stb
    cmd = "inStrain genome_wide -i {0} -s {1}".format(location, BTO.stb)
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
    cmd = "inStrain genome_wide -i {0} -s {1} --skip_mm_profiling".format(location, BTO.stb)
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
    cmd = "inStrain genome_wide -i {0}".format(location, BTO.stb)
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
    cmd = "inStrain genome_wide -i {0} -s {1}".format(location, BTO.single_scaff)
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
    cmd = "inStrain genome_wide -i {0} -s {1} {2} {3}".format(location, BTO.single_scaff,
                                                              BTO.fasta, BTO.extra_single_scaff)
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

def test_genome_level_2(BTO):
    """
    Test running on an RC object
    """
    # Run the dereplicated version
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.RC_Loc))
    shutil.copytree(BTO.RC_Loc, location)

    cmd = "inStrain genome_wide -i {0} -s {1}".format(location, BTO.stb)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    files = glob.glob(IS.get_location('output') + '*')
    files = [f for f in files if 'genomeWide' in f]
    assert len(files) == 2

def test_genome_level_3(BTO):
    """
    Test running on the mm level
    """
    # Run the dereplicated version
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_nogenes))
    shutil.copytree(BTO.IS_nogenes, location)
    for f in glob.glob(location + '/output/*.tsv'):
        os.remove(f)

    cmd = "inStrain genome_wide -i {0} -s {1} --mm_level".format(location, BTO.stb)
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

def test_genome_level_4(BTO):
    """
    Test running RC object on the mm level
    """
    # Run the dereplicated version
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.RC_Loc))
    shutil.copytree(BTO.RC_Loc, location)
    for f in glob.glob(location + '/output/*.tsv'):
        os.remove(f)

    cmd = "inStrain genome_wide -i {0} -s {1} --mm_level".format(location, BTO.stb)
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

def test_genome_level_5(BTO):
    """
    Test iRep methods
    """
    # Make sure it can run on a random thing
    iRep, iRep_junk = inStrain.irep_utilities.calculate_iRep_from_coverage_array(
        np.array([5] * 50000 + [6] * 50000 + [7] * 50000), 1)
    assert abs(iRep - 1.62) < 0.01

    # Load pickled test data
    Test_sets = pickle.load(open(BTO.iRep_test_set, "rb"))

    # Make sure they all work
    for t_set in Test_sets:
        order, length, gc, windows, OLTwidnows, LTwidnows, raw_iRep, iRep = t_set

        # Make windows into df
        Idb = pd.DataFrame({'index': windows[0], 'coverage': windows[1]})
        GCdb = pd.DataFrame({'index': gc[0], 'GC_content': gc[1]})
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

def test_genome_level_6(BTO):
    """
    Test running on database mode (there is a genome with no reads)
    """
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_nogenes))

    # Make sure empty scaffolds don't mess it up
    cmd = "inStrain profile {1} {2} -o {3} -l 0.98 -s {4} --skip_plot_generation".format(
        'junk', BTO.bam1, BTO.extra_single_scaff, location, BTO.stb2)
    print(cmd)
    exit_code = call(cmd, shell=True)
    assert exit_code == 0, exit_code

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
