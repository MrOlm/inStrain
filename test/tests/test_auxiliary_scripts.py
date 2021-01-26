"""
Run tests on auxillary scripts
"""

import os
import numpy as np
import shutil
import glob
from subprocess import call

import pandas as pd

import inStrain
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.profile.fasta
import tests.test_utils as test_utils

from test_utils import BTO

def test_rarefaction_curve_0(BTO):
    """
    Just run the program
    """
    # Run program
    base = BTO.test_dir + 'test'
    sl = test_utils.get_aux_script_loc('rarefaction_curve')

    cmd = f"{sl} -b {BTO.bam1} -o {base} -f {BTO.fasta} -s {BTO.stb}"
    print(cmd)
    call(cmd, shell=True)

    # Load log
    log_loc = os.path.join(base, 'log/log.log')
    assert os.path.isfile(log_loc), log_loc

    # Load output
    outfiles = glob.glob(os.path.join(base, 'output/*'))
    Odb = pd.read_csv([o for o in outfiles if 'rarefaction_table' in o][0])
    Sdb = pd.read_csv([o for o in outfiles if 'subset_table' in o][0])
    assert set(Odb['samtools_subset'].tolist()) == set(Sdb['samtools_subset'].tolist())
