import glob
import pandas as pd

from test_utils import BTO
import test_utils
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
    assert len(rawfiles) >= 2

    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) >= 8

    odb = pd.read_csv([f for f in outfiles if 'LongFormData.csv' in f][0])
    assert len(odb) > 0

    # Make sure not all 0s
    for loc in glob.glob(base + '/output/ParsedGeneAnno_*'):
        Ddb = pd.read_csv(loc)

        got = False
        for c, i in Ddb.items():
            if c == 'sample':
                continue

            if sum(i) > 0:
                got = True

        assert got, loc

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
    assert len(rawfiles) >= 2

    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) >= 8

    odb = pd.read_csv([f for f in outfiles if 'LongFormData.csv' in f][0])
    assert len(odb) > 0


def test_PA_3(BTO):
    """
    Test running parse_annotations without a minimum genomes breadth
    """
    # Run program
    base = BTO.test_dir + 'testA'

    cmd = f"inStrain parse_annotations -i {BTO.IS1} {BTO.IS2} -o {base} -a {BTO.anno_loc} -b 0 --store_rawdata"
    print(cmd)
    inStrain.parse_annotations.PAController(inStrain.argumentParser.parse_args(cmd.split(' ')[1:])).main()

    # Make sure it produced output
    rawfiles = glob.glob(base + '/raw_data/*')
    assert len(rawfiles) >= 3

    outfiles = glob.glob(base + '/output/*')
    assert len(outfiles) < 8

    odb = pd.read_csv([f for f in outfiles if 'LongFormData.csv' in f][0])
    assert len(odb) > 0

def test_PA_4(BTO):
    """
    Test loading and re-running the long-form data
    """
    # Run program
    base = BTO.test_dir + 'testA'

    cmd = f"inStrain parse_annotations -i {BTO.IS1} {BTO.IS2} -o {base} -a {BTO.anno_loc}"
    print(cmd)
    inStrain.parse_annotations.PAController(inStrain.argumentParser.parse_args(cmd.split(' ')[1:])).main()

    # Load the long from data using the API
    outfiles = glob.glob(base + '/output/*')
    lloc = [f for f in outfiles if 'LongFormData.csv' in f][0]
    sloc = [f for f in outfiles if 'SampleAnnotationTotals.csv' in f][0]

    sdb = pd.read_csv(sloc)
    s2a2g2vals = inStrain.parse_annotations.load_s2a2g2vals(lloc)

    # Re-create the output files
    metric2table = inStrain.parse_annotations.create_annotation_tables2(sdb, s2a2g2vals)

    # Make sure they're the same
    outloc = base + '/output/'
    m2n = {'long_data': 'LongFormData.csv'}
    for metric, table in metric2table.items():
        if metric in m2n:
            name = m2n[metric]
        else:
            name = 'ParsedGeneAnno_' + metric + '.csv'

        ori_table = pd.read_csv(outloc + name)
        assert test_utils.compare_dfs2(ori_table, table)

