"""
Methods to help with testing
"""

import os

import pandas as pd


def get_script_loc(script):
    """
    Relies on being run from test_suite.py (not ideal)
    """
    if script == 'inStrain':
        return os.path.join(str(os.getcwd()),
                            '../bin/inStrain')
    if script == 'filter_reads':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/filter_reads.py')
    if script == 'deprecated_gene_statistics':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/deprecated/deprecated_gene_statistics.py')
    if script == 'SNVprofile':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/SNVprofile.py')
    if script == 'readcomparer':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/readComparer.py')
    if script == 'quickProfile':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/quickProfile.py')
    if script == 'GeneProfile':
        return os.path.join(str(os.getcwd()),
                            '../inStrain/GeneProfile.py')


def load_random_test_dir():
    """
    Relies on being run from test_suite.py (not ideal)
    """
    loc = os.path.join(str(os.getcwd()),
                       'test_backend/testdir/')
    return loc


def load_data_loc():
    """
    Relies on being run from test_suite.py (not ideal)
    """
    return os.path.join(str(os.getcwd()),
                        'test_data/')


def get_twelve2thriteen():
    """
    Return information that helps with the 12 -> 13 transition
    """
    # Name translations going from version 1.2 to 1.3
    twelve2thirteen = {
        "baseCoverage": "position_coverage",
        "refBase": "ref_base",
        "conBase": 'con_base',
        "varBase": 'var_base',
        "varFreq": "var_freq",
        "conFreq": "con_freq",
        'refFreq': 'ref_freq',
        'coverage_median': 'coverage_median',
        'std_cov': 'coverage_std',
        'microdiversity': 'nucl_diversity',
        'mean_microdiversity': 'nucl_diversity',
        'median_microdiversity': 'nucl_diversity_median',
        'masked_breadth': 'breadth_minCov',
        'unmaskedBreadth': 'breadth_minCov',
        'rarefied_breadth': 'breadth_rarefied',
        'expected_breadth': 'breadth_expected',
        'rarefied_mean_microdiversity': 'nucl_diversity_rarefied',
        'rarefied_median_microdiversity': 'nucl_diversity_rarefied_median',
        'median_cov': 'coverage_median',
        'Reference_SNPs': 'SNS_count',
        'conANI': 'conANI_reference',
        'popANI': 'popANI_reference',
        'SNPs': 'divergent_site_count',
        'true_length': 'length',
        'total_SNPs': 'divergent_site_count',
        'population_SNPs': 'population_divergent_sites',
        'consensus_SNPs': 'consensus_divergent_sites',
        'filtered_pairs': 'filtered_read_pair_count',
        'SNS_S_count': '',
    }

    # Removed from version 1.3
    del_thirteen = {'bases_w_0_coverage',
                    'mean_clonality',
                    'median_clonality',
                    'BiAllelic_SNPs',
                    'SNPs_per_bp',
                    'MultiAllelic_SNPs',
                    'clonality',
                    'min_ANI'
                    }

    new_thirteen = {'linked_SNV_count',
                    'SNV_distance_mean',
                    'coverage_SEM',
                    'r2_mean',
                    'iRep_GC_corrected',
                    'iRep',
                    'd_prime_mean',
                    'coverage_SEM',
                    'SNV_count',
                    'SNS_S_count',
                    'S_sites',
                    'SNV_N_count',
                    'SNS_N_count',
                    'N_sites',
                    'dNdS_substitutions',
                    'SNV_S_count',
                    'pNpS_variants',
                    'gene_length',
                    'class'
                    }

    return twelve2thirteen, del_thirteen, new_thirteen


def compare_dfs2(db1, db2, verbose=False):
    """
    Return True if dataframes are equal (order of dataframes doesn't matter)
    """
    try:
        # noinspection PyProtectedMember
        pd._testing.assert_frame_equal(db1, db2)
        return True
    except AssertionError as e:
        if verbose:
            print(e)
        return False


def compare_lists(l1, l2, check_order=True):
    if len(l1) != len(l2):
        print("lengths are different")
        return False

    if set(l1) != set(l2):
        print('elements are different')
        return False

    if check_order:
        for i, (o, t) in enumerate(zip(l1, l2)):
            if o != t:
                print('order is different')
                print('at position {0}, list 1 has {1} and 2 has {2}'.format(i, o, t))
                return False

    return True


def _internal_verify_Sdb(Sdb):
    for scaff, d in Sdb.groupby('scaffold'):
        d = d.sort_values('mm')
        for thing in ['breadth_minCov', 'coverage', 'coverage_median']:
            assert d[thing].tolist() == sorted(d[thing].tolist()), d

        # clons = [1-c for c in d['mean_clonality']]
        # micros = [c for c in d['mean_microdiversity']]
        # for c, m in zip(clons, micros):
        #     if c == c:
        #         assert abs(c - m) < 0.0001, [c, m, scaff, d]

    sdb = Sdb[Sdb['breadth_minCov'] > 0]
    for col in sdb.columns:
        if col in ['nucl_diversity_rarefied', 'nucl_diversity_rarefied_median']:
            continue
        if len(sdb[sdb[col].isna()]) > 0:
            print(col)
            print(sdb[sdb[col].isna()])

    for i, row in Sdb.iterrows():
        if row['consensus_divergent_sites'] > 0:
            assert row['conANI_reference'] != 0

    assert Sdb['conANI_reference'].max() <= 1
    assert Sdb['breadth_minCov'].max() <= 1
    assert Sdb['breadth'].max() <= 1

    if 'breadth_rarefied' not in sdb.columns:
        assert len(sdb) == len(sdb.dropna())
    else:
        sdb = sdb[sdb['breadth_rarefied'] > 0]
        assert len(sdb) == len(sdb.dropna())


def _internal_verify_OdbSdb(Odb, Sdb):
    """
    Odb = cumulative scaffold table
    Sdb = SNP table
    """
    # Ensure internal consistancy between Sdb and Cdb at the lowest mm
    low_mm = Sdb['mm'].min()
    for scaff, db in Sdb[Sdb['mm'] == low_mm].groupby('scaffold'):
        snps = Odb['divergent_site_count'][(Odb['scaffold'] == scaff) & (Odb['mm']
                                                                         == low_mm)].fillna(0).tolist()[0]
        assert snps == len(db), [snps, len(db), scaff, low_mm]

    # Ensure internal consistancy between Sdb and Cdb at the highset mm
    odb = Odb.sort_values('mm').drop_duplicates(subset='scaffold', keep='last')
    for scaff, db in Sdb.sort_values('mm').drop_duplicates(subset=['scaffold'
        , 'position'], keep='last').groupby('scaffold'):
        snps = odb['divergent_site_count'][(odb['scaffold'] == scaff)].fillna(0).tolist()[0]
        assert snps == len(db), [snps, len(db), scaff]


def compare_dfs(db1, db2, round=4, verbose=False):
    """
    Return True if dataframes are equal (order of dataframes doesn't matter)
    """

    db1 = db1.fillna(0).round(round)
    db2 = db2.fillna(0).round(round)

    df = pd.concat([db1, db2], sort=True)
    df = df.reset_index(drop=True)
    df_gpby = df.groupby(list(df.columns))
    idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]

    identicle = (len(idx) == 0)
    if (not identicle) and verbose:
        print("index: ", idx)
        print("db1: ", db1)
        print("db2: ", db2)
        print("df_gpby: ", str(df_gpby))

    return identicle


def compare_dicts(d1, d2):
    if len(d1.keys()) != len(d2.keys()):
        print("lengths are different")
        return False

    if set(d1.keys()) != set(d2.keys()):
        print('elements are different')
        return False

    for o, oE in d1.items():
        if d2[o] != oE:
            print("{0} is {1} and {2}".format(o, oE, d2[o]))
            return False

    return True


def compare_covTs(c1, c2):
    old_scaffs = [x for x, i in c1.items() if len(i) > 0]
    new_scaffs = [x for x, i in c2.items() if len(i) > 0]
    if set(old_scaffs) != set(new_scaffs):
        print('Scaffolds are different')
        return False

    for scaffold in set(old_scaffs):
        s1 = c1[scaffold]
        s2 = c2[scaffold]

        if set(s1.keys()) != set(s2.keys()):
            print("Different mms for scaffold {0}".format(scaffold))
            return False

        for mm, s1_vals in s1.items():
            s2_vals = s2[mm]
            if not s1_vals.equals(s2_vals):
                print("Different values for scaffold {0}".format(scaffold))
                return False

    return True
