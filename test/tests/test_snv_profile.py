"""
Run tests on inStrain SNV profile
"""

import importlib
import logging
import os
import shutil
from subprocess import call
import numpy as np

import inStrain
import inStrain.SNVprofile
import tests.test_utils as test_utils

from test_utils import BTO

def test_SNV_profile_0(BTO):
    """
    Test basic things
    """
    # Make an object
    location = BTO.test_dir + 'test.IS'
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

def test_SNV_profile_1(BTO):
    """
    Test changing from old to new version
    """
    # Make a copy of this file
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_v03))
    shutil.copyfile(BTO.IS_v03, location)

    # Run the command
    cmd = "inStrain other --old_IS {1}".format(True, location)
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
        if name == 'location':
            pass

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
            print("Cant compare {0}".format(name))
            # assert nIS.get(name) == getattr(oIS, name), name
