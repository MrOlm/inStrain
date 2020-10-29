#!/usr/bin/env python
"""
Run tests
"""

import warnings

import pandas as pd

import tests.test_compare
import tests.test_filter_reads
import tests.test_genes
import tests.test_genome_level
import tests.test_plotting
import tests.test_profile
import tests.test_quick_profile
import tests.test_snv_profile
import tests.test_special

warnings.filterwarnings("ignore")

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)

# This is from test_compare.py
def UPDATE_COMPARE_TEST_DATA(self):
    """
    Run inStrain on bam1 and bam2, and store the results where IS1 and IS2 are
    """
    # Run the first one
    base1 = self.test_dir + os.path.basename(self.IS1)
    cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(self.bam1, self.fasta,
                                                                          base1)
    print(cmd)
    code = call(cmd, shell=True)
    assert code == 0, code

    # Copy to new location
    if os.path.isdir(self.IS1):
        shutil.rmtree(self.IS1)
    shutil.copytree(base1, self.IS1)

    # Run the second one
    base2 = self.test_dir + os.path.basename(self.IS2)
    cmd = "inStrain profile {0} {1} -o {2} --skip_plot_generation".format(self.bam2, self.fasta,
                                                                          base2)
    print(cmd)
    code = call(cmd, shell=True)
    assert code == 0, code

    # Copy to new location
    if os.path.isdir(self.IS2):
        shutil.rmtree(self.IS2)
    shutil.copytree(base2, self.IS2)

if __name__ == '__main__':
    tests.test_compare.test_readcomparer().run(tests=['18'], cleanUp=True)

    # tests.test_profile.test_profile().run(tests=[16], cleanUp=False)
    # tests.test_compare.test_readcomparer().run(tests=[13], cleanUp=False)

    # tests.test_profile.test_profile().run()
    # tests.test_filter_reads.TestFilterReads().run()
    # tests.test_snv_profile.test_SNVprofile().run()
    # tests.test_genes.test_gene_statistics().run()
    # tests.test_quick_profile.test_quickProfile().run()
    # tests.test_genome_level.test_genome_wide().run()
    # tests.test_plotting.test_plot().run()
    # tests.test_compare.test_readcomparer().run()
    # tests.tests_special.test_special().run()

    print('everything is working swimmingly!')
