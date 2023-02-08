import glob

from test_utils import BTO
from subprocess import call

import inStrain
import inStrain.argumentParser
import inStrain.parse_annotations

def test_PA_1(BTO):
    """
    Test running parse_annotations on default
    """
    # Run program
    base = BTO.test_dir + 'testA'

    cmd = f"inStrain parse_annotations -i {BTO.IS1} {BTO.IS2} -o {base} -a {BTO.anno_loc}"
    print(cmd)
    inStrain.parse_annotations.PAController(inStrain.argumentParser.parse_args(cmd.split(' ')[1:])).main()

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) >= 3

    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) >= 8

def test_PA_2(BTO):
    """
    Test running parse_annotations from the command line
    """
    # Run program
    base = BTO.test_dir + 'testA'

    cmd = f"/Users/mattolm/Programs/inStrain/bin/inStrain parse_annotations -i {BTO.IS1} {BTO.IS2} -o {base} -a {BTO.anno_loc}"
    print(cmd)
    call(cmd, shell=True)

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    print(rawfiles)

def test_PA_3(BTO):
    """
    Test running parse_annotations without a minimum genomes breadth
    """
    # Run program
    base = BTO.test_dir + 'testA'

    cmd = f"inStrain parse_annotations -i {BTO.IS1} {BTO.IS2} -o {base} -a {BTO.anno_loc} -b 0"
    print(cmd)
    inStrain.parse_annotations.PAController(inStrain.argumentParser.parse_args(cmd.split(' ')[1:])).main()

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) >= 3

    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) < 8



