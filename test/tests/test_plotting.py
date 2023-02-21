"""
Run tests on inStrain plotting
"""

import glob
import importlib
import os
import shutil
from subprocess import call
import numpy as np

import inStrain
import inStrain.argumentParser
import inStrain.SNVprofile
import inStrain.deprecated.deprecated_filter_reads
import inStrain.filter_reads
import inStrain.irep_utilities
import inStrain.profile.fasta
import tests.test_utils as test_utils

from test_utils import BTO

def test_plotting_1(BTO, view=False):
    """
    Make sure all plots are made when called from command line
    """
    # Run the IS plots
    FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
            'readANI_distribution.pdf', 'MajorAllele_frequency_plot.pdf',
            'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
            'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
            'GeneHistogram_plot.pdf']

    # !! Maybe when stable you can do this; for now you need to re-run it !!
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_plotting))
    shutil.copytree(BTO.IS_plotting, location)
    for f in glob.glob(location + '/log/*'):
        os.remove(f)

    # Run program
    # location = BTO.test_dir + 'test'
    # cmd = "inStrain profile {1} {2} -o {3} -g {4} -s {5} --skip_plot_generation -p 6 -d".format(
    #         'junk', BTO.bam1, BTO.fasta, location, BTO.genes, BTO.stb)
    # print(cmd)
    # call(cmd, shell=True)
    # os.remove(location + '/log/log.log')

    cmd = "inStrain plot -i {0} -d".format(location)
    print(cmd)

    call(cmd, shell=True)

    # sys_args = cmd.split(' ')
    # args = inStrain.argumentParser.parse_args(sys_args[1:])
    # inStrain.controller.Controller().main(args)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    # Make sure all figures are made
    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000, fig

    # Make sure logging works
    rr = [f for f in glob.glob(location + '/log/*') if 'runtime' in f][0]
    got = False
    with open(rr, 'r') as o:
        for line in o.readlines():
            line = line.strip()
            if 'Plot' in line:
                got += 1
            # print(line)
    assert got == 10, got  # Its in there twice for random reasons

    if view:
        assert False, "Get on the figures here: " + IS.get_location('figures')

def test_plotting_2(BTO):
    """
    Make sure all RC plots are made
    """
    # Run the IS plots
    FIGS = ['inStrainCompare_dendrograms.pdf']

    location = os.path.join(BTO.test_dir, os.path.basename(BTO.RC_Loc))
    shutil.copytree(BTO.RC_Loc, location)

    cmd = "inStrain genome_wide -i {0} -s {1}".format(location, BTO.stb)
    print(cmd)
    call(cmd, shell=True)

    cmd = "inStrain plot -i {0} -d".format(location)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')
    #print(figs)

    # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
    for F in FIGS:
        assert len([f for f in figs if F in f]) >= 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000

def test_plotting_3(BTO):
    """
    Test the breadth cutoff
    """
    # Run the IS plots
    FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
            'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
            'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
            'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
            'GeneHistogram_plot.pdf']

    # FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf']#, 'genomeWide_microdiveristy_metrics.pdf',
    # #         'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
    # #         'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
    # #         'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
    # #         'GeneHistogram_plot.pdf']

    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_plotting))
    shutil.copytree(BTO.IS_plotting, location)

    cmd = "inStrain plot -i {0} -mb 0.97".format(location)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000

def test_plotting_4(BTO):
    """
    Test the genome name
    """
    # Run the IS plots
    FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf', 'genomeWide_microdiveristy_metrics.pdf',
            'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
            'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
            'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
            'GeneHistogram_plot.pdf']

    # FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf']#, 'genomeWide_microdiveristy_metrics.pdf',
    # #         'MajorAllele_frequency_plot.pdf', 'readANI_distribution.pdf',
    # #         'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
    # #         'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
    # #         'GeneHistogram_plot.pdf']

    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_plotting))
    shutil.copytree(BTO.IS_plotting, location)

    cmd = "inStrain plot -i {0} -g maxbin2.maxbin.001.fasta".format(location)
    # cmd = "inStrain plot -i {0} -g poop".format(location)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000

# THIS WORKS ON MO's DEV MACHINE ONLY
def test_plotting_5(BTO):
    """
    Test the genome name
    """
    # Run the IS plots
    FIGS = ['genomeWide_microdiveristy_metrics.pdf',
            'MajorAllele_frequency_plot.pdf',
            'LinkageDecay_plot.pdf', 'ScaffoldInspection_plot.pdf',
            'ReadFiltering_plot.pdf', 'LinkageDecay_types_plot.pdf',
            'GeneHistogram_plot.pdf']

    if not os.path.isdir(BTO.LOCALONLY_IS_plotting):
        print("THIS TEST ONLY WORKS ON MO'S DEVELOPMENT MACHINE")
        return

    location = BTO.LOCALONLY_IS_plotting
    for f in glob.glob(location + '/log/*'):
        os.remove(f)
    for f in glob.glob(location + '/figures/*'):
        os.remove(f)

    # Pick a range of genomes
    g = "GUT_GENOME000001.fna GUT_GENOME057980.fna GUT_GENOME143760.fna GUT_GENOME265791.fna"

    cmd = f"inStrain plot -d -i {location} -g {g}"
    # cmd = "inStrain plot -i {0} -g poop".format(location)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    # assert sorted([os.path.basename(f) for f in figs]) == sorted(FIGS), set(FIGS) - set([os.path.basename(f) for f in figs])
    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000

def test_plot_10_squareform(BTO):
    """
    Test the ability of plot 10 to handle empty comparisons
    """
    # Run the IS plots
    FIGS = ['inStrainCompare_dendrograms.pdf']

    location = os.path.join(BTO.test_dir, os.path.basename(BTO.plot10_tester))
    shutil.copytree(BTO.plot10_tester, location)

    cmd = "inStrain plot -i {0} -d -pl 10".format(location)
    print(cmd)
    call(cmd, shell=True)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 2000


def test_BreadthCurve_plot_1(BTO, view=False):
    """
    Make sure plot1 is made from command line
    """
    # Run the IS plots
    FIGS = ['CoverageAndBreadth_vs_readMismatch.pdf']

    # !! Maybe when stable you can do this; for now you need to re-run it !!
    location = os.path.join(BTO.test_dir, os.path.basename(BTO.IS_plotting))
    shutil.copytree(BTO.IS_plotting, location)
    for f in glob.glob(location + '/log/*'):
        os.remove(f)

    cmd = "inStrain plot -i {0} -pl 1 2 -d".format(location)
    print(cmd)
    sys_args = cmd.split(' ')
    args = inStrain.argumentParser.parse_args(sys_args[1:])
    inStrain.controller.Controller().main(args)

    # Load output
    IS = inStrain.SNVprofile.SNVprofile(location)
    figs = glob.glob(IS.get_location('figures') + '*')

    # Make sure all figures are made
    for F in FIGS:
        assert len([f for f in figs if F in f]) == 1, F
    for fig in figs:
        assert os.path.getsize(fig) > 1000, fig

    view = False
    if view:
        assert False, "Get on the figures here: " + IS.get_location('figures')