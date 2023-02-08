"""
Methods to help with testing
"""

import os

import pandas as pd
import shutil
import logging
import importlib

import pytest

class TestingClass():
    def teardown(self):
        importlib.reload(logging)
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)
        importlib.reload(logging)

@pytest.fixture()
def BTO():
    """
    Basic test object

    This object makes no copies of anything; just has references and does setup / cleanup
    """
    # Set up
    self = TestingClass()
    dl = load_data_loc()

    self.test_dir = load_random_test_dir()

    self.fasta =  dl + 'N5_271_010G1_scaffold_min1000.fa'
    self.failure_fasta = dl + 'N5_271_010G1_scaffold_failureScaffold.fa'
    self.highcov_fasta = dl + 'L3_108_000G1_scaffold_35331.fasta'
    self.single_scaff = dl + 'N5_271_010G1_scaffold_101.fasta'
    self.extra_single_scaff = dl + 'N5_271_010G1_scaffold_101_extra.fasta'
    self.fasta_extra = dl + 'N5_271_010G1_scaffold_min1000_extra.fa'
    self.small_fasta = dl + 'SmallScaffold.fa'
    self.bbfasta = dl +  'N5_271_010G1_scaffold_963.fasta'

    self.bam1 = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
    self.bam2 = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam'
    self.failure_bam = dl + 'N5_271_010G1_scaffold_failureScaffold.sorted.bam'
    self.highcov_bam = dl + 'NIHL.delta.fasta-vs-L2_019_000G1.L3_108_000G1_scaffold_35331.sorted.bam'
    self.small_bam = dl + 'SmallScaffold.fa.sorted.bam'
    self.bbsam = dl + 'bbmap_N5_271_010G1_scaffold_963.fasta.sorted.bam'

    self.IS = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS'
    self.IS1 = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.forRC.IS'
    self.IS2 = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.forRC.IS'
    self.highcov_IS = dl + 'L3_108_000G1_scaffold_35331.IS'
    self.old_IS = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.v4.IS'
    self.IS_nogenes = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.IS.v1.3.0g'
    self.IS_nogenes2 = dl + 'sars_cov_2_MT039887.1.fasta.bt2-vs-SRR11140750.sam.IS'
    self.IS_plotting = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.ForPlotting.IS'
    self.IS_v03 = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.v3.IS.pickle'

    self.RC_Loc = dl + 'readComparer_vCurrent.RC'
    self.v12_solution = dl + 'readComparer_v1.3.0r.RC'
    self.plot10_tester = dl + 'plot10bugtester.RC'

    self.SIS = dl + 'Ecoli_ani.100.0.subset.sorted.bam.IS'
    self.SIS2 = dl + 'Ecoli_ani.99.9.subset.sorted.bam.IS'
    self.SIS3 = dl + 'Ecoli_ani.98.0.subset.sorted.bam.IS'

    self.scafflist = dl + 'scaffList.txt'
    self.scafflistF = dl + 'scaffList.fasta'

    self.stb = dl + 'N5_271_010G1.maxbin2.stb'
    self.stb2 = dl + 'GenomeCoverages.stb'

    self.readloc = dl + 'N5_271_010G1.R1.fastq.gz'
    self.readloc2 = dl + 'N5_271_010G1.R1.reads'

    self.old_genes_script = get_script_loc('deprecated_gene_statistics')
    self.genes = dl + 'N5_271_010G1_scaffold_min1000.fa.genes.fna'
    self.genbank = dl + 'sars_cov_2_MT039887.1.gb'
    self.failure_genes = dl + 'N5_271_010G1_scaffold_failureScaffold.fa.genes.fna.fa'
    self.anno_loc = dl + 'N5_271_010G1_scaffold_min1000.fa.genes.faa.KOfam_anno.txt'

    self.iRep_test_set = dl + 'test_iRep.p'

    self.cc_solution = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam.CB'
    self.pp_snp_solution = dl + 'strainProfiler_v0.3_results/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted' \
                                '.bam_SP_snpLocations.pickle'

    self.sam = dl + 'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam'
    self.special_bam = dl + 'SP_sam.bam'

    self.numbam = dl + 'liver.bam'
    self.numfasta = dl + 'transcriptome_chopped.fa'

    self.LOCALONLY_IS_plotting = '/Users/mattolm/Programs/testing_house/UHGG_reps.fasta-vs-N5_216_039G1.a.IS'

    self.teardown()
    yield self
    self.teardown()


def get_script_loc(script):
    """
    Relies on being run from test_docker.py (not ideal)
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

def get_aux_script_loc(script):
    """
    Relies on being run from within inStrain/test (not ideal)
    """
    if script == 'rarefaction_curve':
        return os.path.join(str(os.getcwd()),
                            '../auxiliary_scripts/rarefaction_curve.py')


def load_random_test_dir():
    """
    Relies on being run from test_docker.py (not ideal)
    """
    loc = os.path.join(str(os.getcwd()),
                       'test_backend/testdir/')
    return loc


def load_data_loc():
    """
    Relies on being run from test_docker.py (not ideal)
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

    for x, y in zip(db1.iterrows(), db2.iterrows()):
        i, row1 = x
        i, row2 = y

        for col in df.columns:
            if row1[col] != row2[col]:
                print(f"1 has {row1}, 2 has {row2}")

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
