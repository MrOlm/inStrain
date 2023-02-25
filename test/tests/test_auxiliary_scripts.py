"""
Run tests on auxillary scripts
"""

import os
import numpy as np
import shutil
import pytest
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

@pytest.mark.skipif(not inStrain.utils.find_program('sambamba')[1], reason="requires sambamba")
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

def test_recluster_instrain(BTO):
    """
    Just run the program on two settings to make sure the -a arguments works
    """
    # Run program
    outfile = BTO.test_dir + 'test.tsv'
    sl = test_utils.get_aux_script_loc('recluster_instrain_compare')
    db_loc = os.path.join(BTO.RC_Loc, 'output/RC_test_genomeWide_compare.tsv')

    cmd = f"{sl} --version"
    print(cmd)
    call(cmd, shell=True)

    cmd = f"{sl} -d {db_loc} -o {outfile} -a 99"
    print(cmd)
    call(cmd, shell=True)

    # Load output
    outfiles = glob.glob(outfile + '*')
    assert len(outfiles) == 1
    Odb = pd.read_csv(outfiles[0], sep='\t')
    assert len(Odb['cluster'].unique()) == 2

    cmd = f"{sl} -d {db_loc} -o {outfile} -a 99.995"
    print(cmd)
    call(cmd, shell=True)

    # Load output
    outfiles = glob.glob(outfile + '*')
    assert len(outfiles) == 1
    Odb = pd.read_csv(outfiles[0], sep='\t')
    assert len(Odb['cluster'].unique()) == 3
