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

from tests.test_utils import BTO

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


def _load_refilter_script():
    import importlib.util
    loc = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../auxiliary_scripts/refilter_profile.py')
    spec = importlib.util.spec_from_file_location('refilter_profile', loc)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _profile_for_refilter(BTO, name, extra=''):
    import importlib
    import logging
    import inStrain.argumentParser
    import inStrain.controller
    importlib.reload(logging)
    out = os.path.join(BTO.test_dir, name)
    cmd = f"profile {BTO.bam1} {BTO.fasta} -o {out} -s {BTO.stb} -g {BTO.genes} -p 6 --skip_plot_generation {extra}"
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split()))
    return out


def test_refilter_profile(BTO):
    """
    Re-filtering a default profile to stricter thresholds should give the same results as profiling with those
    thresholds in the first place (except linkage, which needs the reads)
    """
    import importlib
    import logging
    refilter = _load_refilter_script()

    default = _profile_for_refilter(BTO, 'default.IS')
    strict = _profile_for_refilter(BTO, 'strict.IS', '-c 10 -f 0.1')

    importlib.reload(logging)
    refiltered = os.path.join(BTO.test_dir, 'refiltered.IS')
    refilter.refilter_profile(refilter.parse_args(
        ['-i', default, '-o', refiltered, '-c', '10', '-f', '0.1', '-g', BTO.genes]))

    # The re-filtered profile must actually differ from the default one
    n_default = len(inStrain.SNVprofile.SNVprofile(default).get_nonredundant_snv_table())
    n_strict = len(inStrain.SNVprofile.SNVprofile(strict).get_nonredundant_snv_table())
    assert n_strict < n_default, (n_strict, n_default)

    # Linkage needs the reads, so it isn't re-calculated. Rarefied nucleotide diversity is a random subsample that
    # differs between any two runs (and doesn't depend on min_cov / min_freq), so it can't be compared either
    linkage_columns = ['linked_SNV_count', 'SNV_distance_mean', 'r2_mean', 'd_prime_mean',
                       'nucl_diversity_rarefied', 'nucl_diversity_rarefied_median']
    for table, key in [('SNVs', ['scaffold', 'position']), ('scaffold_info', ['scaffold']),
                       ('genome_info', ['genome']), ('gene_info', ['gene'])]:
        e = inStrain.SNVprofile.SNVprofile(strict).generate(table, store=False, return_table=True)
        g = inStrain.SNVprofile.SNVprofile(refiltered).generate(table, store=False, return_table=True)
        e = e.drop(columns=[c for c in linkage_columns if c in e.columns])
        g = g.drop(columns=[c for c in linkage_columns if c in g.columns])
        assert sorted(e.columns) == sorted(g.columns), (table, set(e.columns) ^ set(g.columns))
        e = e.sort_values(key).reset_index(drop=True)
        g = g.sort_values(key).reset_index(drop=True)[list(e.columns)]
        pd.testing.assert_frame_equal(e, g, check_dtype=False, check_categorical=False, obj=table)

    # Thresholds can't be lowered
    with pytest.raises(ValueError):
        refilter.refilter_profile(refilter.parse_args(
            ['-i', default, '-o', os.path.join(BTO.test_dir, 'lower.IS'), '-c', '3', '-f', '0.1']))
