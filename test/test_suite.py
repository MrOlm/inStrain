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
import tests.tests_special

warnings.filterwarnings("ignore")

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)


if __name__ == '__main__':
    # tests.test_profile.test_profile().run(tests=[16])
    # tests.test_compare.test_readcomparer().run(tests=[13])

    tests.test_profile.test_profile().run()
    tests.test_filter_reads.TestFilterReads().run()
    tests.test_snv_profile.test_SNVprofile().run()
    tests.test_genes.test_gene_statistics().run()
    tests.test_quick_profile.test_quickProfile().run()
    tests.test_genome_level.test_genome_wide().run()
    tests.test_plotting.test_plot().run()
    tests.test_compare.test_readcomparer().run()
    tests.tests_special.test_special().run()

    print('everything is working swimmingly!')
