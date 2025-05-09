#!/usr/bin/env python

import os
import sys
import json
import h5py
import time
import psutil
import pickle
import logging
import argparse
import traceback
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing
import concurrent.futures
from concurrent import futures
from collections import defaultdict

from Bio import SeqIO

if __name__ != '__main__':
    from ._version import __version__

import inStrain.profile.snv_utilities
import inStrain.profile.fasta

import inStrain.SNVprofile
import inStrain.controller
import inStrain.logUtils
import inStrain.polymorpher

C2P = {0:'A', 1:'C', 2:'T', 3:'G'}
def compare_scaffold(scaffold, cur_names, SNPtables, covTs, mLen, null_model, **kwargs):
    '''
    This is the money method thats going to be multithreaded eventually

    Arguments:
        scaffold: name of scaffold (needs to match the SNPtables)
        names: names of the samples; must be in the same order as covTs and SNPtables
        sProfiles: a list of SNVprofiles that correspond to the names list
        mLen: length of the scaffold

    Returns:
        [DataFrame of "scaffold, sample1, sample2, mm, coverage_overlap, ANI",
        DataFrame of SNP locs,
        pair2mm2covOverlap,
        name of scaffold], log
    '''
    log_message = inStrain.logUtils.get_worker_log('Compare', scaffold, 'start')

    # Load arguments
    min_cov = kwargs.get('min_cov', 5)
    min_freq = kwargs.get('min_freq', 5)
    debug = kwargs.get('debug', False)
    fdr = kwargs.get('fdr', 1e-6)
    store_coverage = kwargs.get('store_coverage_overlap', False)
    store_mm_locations = kwargs.get('store_mismatch_locations', False)
    include_self_comparisons = kwargs.get('include_self_comparisons', False)

    # Load pooling arguments
    run_pooling = kwargs.get('run_pooling', False)
    name2Rdic = kwargs.get('name2Rdic', {})
    name2bam_loc = kwargs.get('name2bamloc', {})

    # For testing purposes
    if ((scaffold == 'FailureScaffoldHeaderTesting') & (debug)):
        assert False

    if len(cur_names) < 2:
        results = [pd.DataFrame(), pd.DataFrame(), {}, 'skip_{0}'.format(scaffold)]
        log_message += '\n' + inStrain.logUtils.get_worker_log('Compare', scaffold, 'end')
        return (results, log_message)

    # Iterate through pairs
    i = 0
    snpLocs = []
    pair2mm2covOverlap = {}
    table = defaultdict(list)
    for covT1, SNPtable1_ori, name1 in zip(covTs, SNPtables, cur_names):
        i += 1
        j = 0
        for covT2, SNPtable2_ori, name2 in zip(covTs, SNPtables, cur_names):
            j += 1

            if i > j:
                continue

            if (not include_self_comparisons) & (i == j):
                continue

            logging.debug("{2} {0} vs {1} ({3} {4})".format(name1, name2, scaffold,
                        i, j))

            if debug:
                pid = os.getpid()
                process  = psutil.Process(os.getpid())
                bytes_used = process.memory_info().rss
                total_available_bytes = psutil.virtual_memory()
                nm = "\n{4} PID {0} end at {5} with {1} RAM. System has {2} of {3} available".format(
                        pid, bytes_used, total_available_bytes[1], total_available_bytes[0],
                        scaffold, time.time())
                logging.debug(nm)

            mm2overlap, mm2coverage = calc_mm2overlap_v3(covT1, covT2, min_cov=min_cov, verbose=False, debug=debug)
            Mdb = _calc_SNP_count_v3(SNPtable1_ori, SNPtable2_ori, mm2overlap, null_model, min_freq=min_freq, debug=debug)

            table = _update_overlap_table(table, scaffold, mm2overlap, mm2coverage, Mdb, name1, name2, mLen)

            if store_mm_locations:
                Mdb['name1'] = name1
                Mdb['name2'] = name2
                Mdb['scaffold'] = scaffold
                snpLocs.append(Mdb)

            if store_coverage:
                pair2mm2covOverlap['-vs-'.join(sorted([name1, name2]))] = mm2overlap

    # logging.debug("Returning {0} {1} {2}".format(scaffold, i, j))

    Cdb = pd.DataFrame(table)
    # if len(Cdb) > 0:
    #     Cdb['coverage_overlap'] = Cdb['coverage_overlap'].astype('int64')

    if len(snpLocs) > 0:
        Mdb = pd.concat(snpLocs, sort=False)
    else:
        Mdb = pd.DataFrame()

    # Run pooling analysis
    if run_pooling:
        PM = inStrain.polymorpher.PoolController(SNPtables, cur_names, name2Rdic, name2bam_loc, scaffold)
        PM.main()
        DSTdb = PM.DDST
        PMdb = PM.PST
        results = [Cdb, Mdb, pair2mm2covOverlap, scaffold, DSTdb, PMdb]

    else:
        results = [Cdb, Mdb, pair2mm2covOverlap, scaffold]

    log_message += '\n' + inStrain.logUtils.get_worker_log('Compare', scaffold, 'end')
    return (results, log_message)

def calc_mm2overlap_v3(covT1, covT2, min_cov=5, verbose=False, debug=False):
    '''
    Calculate mm2overlap for a pair of covTs

    Coverage is calculated as cov = len(coveredInBoth) / len(coveredInEither)
    This means that its the percentage of bases that are covered by both

    Returns:
        mm2overlap -> dictionary to array of "True" where there's overlap and "False" where both are compared but there's not overlap
        mm2coverage -> dictionary of mm -> the alignment coverage
    '''
    mm2overlap = {}
    mm2coverage = {}

    # if debug != False:
    #     scaffold, name1, name2 = debug

    mms = sorted(list(set(covT1.keys()).union(set(covT2.keys()))))
    cov1 = pd.Series(dtype='float64')
    cov2 = pd.Series(dtype='float64')
    for mm in mms:
        if mm in covT1:
            cov1 = cov1.add(covT1[mm], fill_value=0)
        if mm in covT2:
            cov2 = cov2.add(covT2[mm], fill_value=0)

        # Figure out where each has min coverage
        T1 = cov1[(cov1 >= min_cov)].index
        T2 = cov2[(cov2 >= min_cov)].index

        # Figure out the total possible overlap
        coveredInEither = np.union1d(T1, T2)#, assume_unique=True)
        coveredInBoth = np.intersect1d(T1, T2, assume_unique=True)

        # Calculate coverage
        if len(coveredInEither) > 0:
            cov = len(coveredInBoth) / len(coveredInEither)
        else:
            cov = 0

        # Save
        mm2overlap[mm] = coveredInBoth
        mm2coverage[mm] = cov

    return mm2overlap, mm2coverage

def calc_mm2overlap(covT1, covT2, min_cov=5, verbose=False, debug=False):
    '''
    Calculate mm2overlap for a pair of covTs

    Coverage is calculated as cov = len(coveredInBoth) / len(coveredInEither)
    This means that its the percentage of bases that are covered by both

    Returns:
        mm2overlap -> dictionary to array of "True" where there's overlap and "False" where both are compared but there's not overlap
        mm2coverage -> dictionary of mm -> the alignment coverage
    '''
    mm2overlap = {}
    mm2coverage = {}

    # if debug != False:
    #     scaffold, name1, name2 = debug

    mms = sorted(list(set(covT1.keys()).union(set(covT2.keys()))))
    cov1 = pd.Series(dtype='float64')
    cov2 = pd.Series(dtype='float64')
    for mm in mms:
        if mm in covT1:
            cov1 = cov1.add(covT1[mm], fill_value=0)
        if mm in covT2:
            cov2 = cov2.add(covT2[mm], fill_value=0)

        # Figure out where each has min coverage
        T1 = set(cov1[(cov1 >= min_cov)].index)
        T2 = set(cov2[(cov2 >= min_cov)].index)

        # Figure out the total possible overlap
        coveredInEither = T1.union(T2)

        # Figure out where there's overlap in both
        coveredInBoth = T1.intersection(T2)

        # Calculate coverage
        if len(coveredInEither) > 0:
            cov = len(coveredInBoth) / len(coveredInEither)
        else:
            cov = 0

        # Save
        mm2overlap[mm] = coveredInBoth
        mm2coverage[mm] = cov

    return mm2overlap, mm2coverage

def _gen_blank_Mdb(COLUMNS):
    '''
    COLUMNS = ['position', 'consensus_SNP', 'population_SNP', 'mm'
    'con_base_1', 'ref_base_1', 'var_base_1', 'position_coverage_1',
    'A_1', 'C_1', 'T_1', 'G_1',
    'con_base_2', 'ref_base_2', 'var_base_2', 'position_coverage_2',
    'A_2', 'C_2', 'T_2', 'G_2']
    '''

    return pd.DataFrame({c:[] for c in COLUMNS})

def _calc_SNP_count_v3(SNPtable1, SNPtable2, mm2overlap, null_model, min_freq=.05, debug=False):

    mm2ANI = {}
    mm2popANI = {}
    dbs = []

    # Get the null model for SNP calling
    model_to_use = null_model

    # Constant for the SNP dataframe
    SNP_COLUMNS = ['position', 'con_base', 'ref_base', 'var_base',
    'position_coverage', 'A', 'C', 'T', 'G']

    # Constant for the output dataframe
    OUT_COLUMNS = ['position', 'consensus_SNP', 'population_SNP', 'mm',
    'con_base_1', 'ref_base_1', 'var_base_1', 'position_coverage_1',
    'A_1', 'C_1', 'T_1', 'G_1',
    'con_base_2', 'ref_base_2', 'var_base_2', 'position_coverage_2',
    'A_2', 'C_2', 'T_2', 'G_2']

    RENAME_COLUMNS = ['con_base', 'ref_base', 'var_base', 'position_coverage',
                    'A', 'C', 'T', 'G']

    # Iterate mm levels
    for mm, cov_arr in mm2overlap.items():

        # Subset to bases that have coverage in both
        covs = set(cov_arr)
        #covs = cov_arr

        # These represent relevant counts at these posisions
        if len(SNPtable1) > 0:
            s1_all = SNPtable1[[(p in covs) for p in SNPtable1['position'].values]].drop_duplicates(
                        subset=['position'], keep='last')
            del s1_all['mm']
            if len(s1_all) == 0:
                s1_all = None
        else:
            #s1_all = _gen_blank_SNPdb(SNP_COLUMNS)
            s1_all = None

        if len(SNPtable2) > 0:
            s2_all = SNPtable2[[(p in covs) for p in SNPtable2['position'].values]].drop_duplicates(
                        subset=['position'], keep='last')
            del s2_all['mm']
            if len(s2_all) == 0:
                s2_all = None
        else:
            #s2_all = _gen_blank_SNPdb(SNP_COLUMNS)
            s2_all = None

        # Merge
        if (s1_all is None) & (s2_all is None):
            #Mdb = _gen_blank_Mdb(OUT_COLUMNS)
            Mdb = None

        elif s1_all is None:
            Mdb = s2_all.rename(columns={c:c + '_2' for c in RENAME_COLUMNS})
            for c in RENAME_COLUMNS:
                Mdb[c + '_1'] = np.nan

        elif s2_all is None:
            Mdb = s1_all.rename(columns={c:c + '_1' for c in RENAME_COLUMNS})
            for c in RENAME_COLUMNS:
                Mdb[c + '_2'] = np.nan

        else:
            Mdb = pd.merge(s1_all, s2_all, on='position', suffixes=('_1', '_2'), how='outer', copy=False)

        if Mdb is not None:
            # Vectorized conSNP calculation
            Mdb = calc_con_snps(Mdb)

            # Vectorized popSNP calculation
            Mdb = calc_pop_snps(Mdb, model_to_use, min_freq)

            

            Mdb['mm'] = mm
            Mdb = Mdb[OUT_COLUMNS]

            # Only keep SNPs
            Mdb = Mdb[Mdb['consensus_SNP'] | Mdb['population_SNP'] ]

            dbs.append(Mdb)

    if len(dbs) > 0:
        Mdb = pd.concat(dbs, sort=False)
    else:
        Mdb = _gen_blank_Mdb(OUT_COLUMNS)
    return Mdb

def call_pop_snps_vectorized_optimized(df, model_to_use, min_freq):
    '''
    Optimized vectorized function to call population SNPs
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing required columns for SNP identification
    model_to_use : dict
        Model dictionary for minimum count thresholds
    min_freq : float
        Minimum frequency threshold
        
    Returns:
    --------
    pandas.Series
        Boolean Series indicating population SNPs
    '''
    import pandas as pd
    import numpy as np
    
    # Initialize the result series with False
    result = pd.Series(False, index=df.index)
    
    # CASE 1: Are the consensus bases the same? If yes, not a SNP
    same_consensus = df['con_base_1'] == df['con_base_2']
    
    # CASE 2: Handle where one of the consensus bases is NaN
    con1_is_nan = df['con_base_1'].isna()
    con2_is_nan = df['con_base_2'].isna()
    one_consensus_nan = con1_is_nan | con2_is_nan
    
    # For rows where only con_base_1 is NaN, check if ref_base_2 is present in sample 2
    con1_nan_only = con1_is_nan & ~con2_is_nan
    if con1_nan_only.any():
        # Extract subset of data for this case
        subset = df.loc[con1_nan_only]
        
        # Vectorized approach for getting counts
        ref_bases = subset['ref_base_2'].values
        counts = pd.Series(0.0, index=subset.index)
        
        # Pre-compute all possible column names for faster lookups
        base_cols = [col for col in subset.columns if col.endswith('_2') and len(col) == 3]
        
        # Create a mapping of base to column for quick lookup
        base_to_col = {col[0]: col for col in base_cols}
        
        # Use numpy for faster assignment
        for base, col in base_to_col.items():
            mask = ref_bases == base
            if mask.any() and col in subset.columns:
                counts.loc[subset.index[mask]] = subset.loc[subset.index[mask], col].values
        
        totals = subset['position_coverage_2']
        # Check if the reference base is present
        is_present_mask = is_present_vectorized(counts, totals, model_to_use, min_freq)
        # If reference is NOT present, it's a population SNP
        result.loc[con1_nan_only] = ~is_present_mask
    
    # For rows where only con_base_2 is NaN, check if ref_base_1 is present in sample 1
    con2_nan_only = ~con1_is_nan & con2_is_nan
    if con2_nan_only.any():
        # Extract subset of data for this case
        subset = df.loc[con2_nan_only]
        
        # Vectorized approach for getting counts
        ref_bases = subset['ref_base_1'].values
        counts = pd.Series(0.0, index=subset.index)
        
        # Pre-compute all possible column names for faster lookups
        base_cols = [col for col in subset.columns if col.endswith('_1') and len(col) == 3]
        
        # Create a mapping of base to column for quick lookup
        base_to_col = {col[0]: col for col in base_cols}
        
        # Use numpy for faster assignment
        for base, col in base_to_col.items():
            mask = ref_bases == base
            if mask.any() and col in subset.columns:
                counts.loc[subset.index[mask]] = subset.loc[subset.index[mask], col].values
        
        totals = subset['position_coverage_1']
        # Check if the reference base is present
        is_present_mask = is_present_vectorized(counts, totals, model_to_use, min_freq)
        # If reference is NOT present, it's a population SNP
        result.loc[con2_nan_only] = ~is_present_mask
    
    # CASE 3: Handle where both consensus bases are present but different
    both_present_diff = ~same_consensus & ~one_consensus_nan
    if both_present_diff.any():
        subset = df.loc[both_present_diff]
        
        # Check if con_base_1 is present in sample 2
        con_bases_1 = subset['con_base_1'].values
        counts_1_in_2 = pd.Series(0.0, index=subset.index)
        
        # Pre-compute mappings for faster lookups
        base_cols_2 = [col for col in subset.columns if col.endswith('_2') and len(col) == 3]
        base_to_col_2 = {col[0]: col for col in base_cols_2}
        
        for base, col in base_to_col_2.items():
            mask = con_bases_1 == base
            if mask.any() and col in subset.columns:
                counts_1_in_2.loc[subset.index[mask]] = subset.loc[subset.index[mask], col].values
        
        is_1_in_2 = is_present_vectorized(counts_1_in_2, subset['position_coverage_2'], model_to_use, min_freq)
        
        # Check if con_base_2 is present in sample 1
        con_bases_2 = subset['con_base_2'].values
        counts_2_in_1 = pd.Series(0.0, index=subset.index)
        
        # Pre-compute mappings for faster lookups
        base_cols_1 = [col for col in subset.columns if col.endswith('_1') and len(col) == 3]
        base_to_col_1 = {col[0]: col for col in base_cols_1}
        
        for base, col in base_to_col_1.items():
            mask = con_bases_2 == base
            if mask.any() and col in subset.columns:
                counts_2_in_1.loc[subset.index[mask]] = subset.loc[subset.index[mask], col].values
        
        is_2_in_1 = is_present_vectorized(counts_2_in_1, subset['position_coverage_1'], model_to_use, min_freq)
        
        # Check for shared minor alleles - vectorized approach
        shared_minor = pd.Series(False, index=subset.index)
        
        # Check which column format is used ('morphia' or 'allele_count')
        if 'allele_count_1' in df.columns:
            multi_allele_1 = subset['allele_count_1'] > 1
            multi_allele_2 = subset['allele_count_2'] > 1
            both_multi = multi_allele_1 & multi_allele_2
            
            if both_multi.any():
                shared_minor.loc[both_multi] = subset.loc[both_multi, 'var_base_1'].values == subset.loc[both_multi, 'var_base_2'].values
                
        elif 'morphia_1' in df.columns:
            multi_allele_1 = subset['morphia_1'] > 1
            multi_allele_2 = subset['morphia_2'] > 1
            both_multi = multi_allele_1 & multi_allele_2
            
            if both_multi.any():
                shared_minor.loc[both_multi] = subset.loc[both_multi, 'var_base_1'].values == subset.loc[both_multi, 'var_base_2'].values
        
        # Determine if it's a population SNP: both consensus bases are different,
        # neither is present in the other sample, and they don't share a minor allele
        is_pop_snp = ~is_1_in_2 & ~is_2_in_1 & ~shared_minor
        result.loc[both_present_diff] = is_pop_snp
    
    return result


def calc_con_snps(df):
    '''
    Vectorized function to call SNPs based on consensus sequences
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing con_base_1, con_base_2, ref_base_1, ref_base_2 columns
        
    Returns:
    --------
    pandas.DataFrame
        The input DataFrame with 'consensus_SNP' column added
    '''
    # Extract all needed columns at once to avoid multiple DataFrame lookups
    con_base_1 = df['con_base_1']
    con_base_2 = df['con_base_2']
    ref_base_1 = df['ref_base_1']
    ref_base_2 = df['ref_base_2']
    
    # Create masks for NaN values - reuse them to avoid recalculation
    sample1_nan = con_base_1.isna()
    sample2_nan = con_base_2.isna()
    
    # Calculate the results directly without intermediate variables
    # Initialize with False array of same size as dataframe
    result = pd.Series(False, index=df.index)
    
    # Set result where both samples are present (most common case first)
    non_nan_mask = ~(sample1_nan | sample2_nan)
    if non_nan_mask.any():
        result.loc[non_nan_mask] = con_base_1.loc[non_nan_mask] != con_base_2.loc[non_nan_mask]
    
    # Set result where only sample2 is present (sample1 is NaN)
    if sample1_nan.any():
        result.loc[sample1_nan] = con_base_2.loc[sample1_nan] != ref_base_2.loc[sample1_nan]
    
    # Set result where only sample1 is present (sample2 is NaN)
    if sample2_nan.any():
        result.loc[sample2_nan] = con_base_1.loc[sample2_nan] != ref_base_1.loc[sample2_nan]
    
    # Assign the result directly without copying
    df['consensus_SNP'] = result
    
    return df


def is_present_vectorized(counts, totals, model_to_use, min_freq):
    """
    Optimized vectorized function to check if counts meet the threshold criteria
    
    Parameters:
    -----------
    counts : pandas.Series
        Count values
    totals : pandas.Series
        Total coverage values
    model_to_use : dict
        Model dictionary for minimum count thresholds
    min_freq : float
        Minimum frequency threshold
        
    Returns:
    --------
    pandas.Series
        Boolean Series indicating if counts meet threshold criteria
    """
    import numpy as np
    
    # Calculate frequencies directly (no intermediate Series needed)
    freqs = counts.values / totals.values
    
    # Convert model_to_use to sorted arrays for faster lookup
    coverage_levels = np.array(sorted(model_to_use.keys()), dtype=np.float64)
    threshold_values = np.array([model_to_use[k] for k in coverage_levels], dtype=np.float64)
    
    # Find appropriate thresholds in one vectorized operation
    # Use np.maximum to handle the case where totals are less than the smallest coverage level
    indices = np.maximum(0, np.searchsorted(coverage_levels, totals.values, side='right') - 1)
    thresholds = threshold_values[indices]
    
    # Combine conditions directly with numpy operations
    # Both operations at once: (counts >= thresholds) & (freqs >= min_freq)
    is_present = (counts.values >= thresholds) & (freqs >= min_freq)
    
    # Convert back to pandas Series with the original index
    return pd.Series(is_present, index=counts.index)

# Usage:
def calc_pop_snps(Mdb, model_to_use, min_freq):
    '''
    Apply vectorized population SNP calling to a DataFrame
    
    Parameters:
    -----------
    Mdb : pandas.DataFrame
        Input DataFrame
    model_to_use : dict
        Model for minimum count thresholds
    min_freq : float
        Minimum frequency threshold
        
    Returns:
    --------
    pandas.DataFrame
        Input DataFrame with 'population_SNP' column added
    '''
    Mdb.loc[:, 'population_SNP'] = call_pop_snps_vectorized_optimized(Mdb, model_to_use, min_freq)
    return Mdb

# def _calc_SNP_count_alternate(SNPtable1, SNPtable2, mm2overlap, null_model, min_freq=.05, debug=False):

#     mm2ANI = {}
#     mm2popANI = {}
#     dbs = []

#     # Get the null model for SNP calling
#     model_to_use = null_model

#     # Constant for the SNP dataframe
#     SNP_COLUMNS = ['position', 'con_base', 'ref_base', 'var_base',
#     'position_coverage', 'A', 'C', 'T', 'G']

#     # Constant for the output dataframe
#     OUT_COLUMNS = ['position', 'consensus_SNP', 'population_SNP', 'mm',
#     'con_base_1', 'ref_base_1', 'var_base_1', 'position_coverage_1',
#     'A_1', 'C_1', 'T_1', 'G_1',
#     'con_base_2', 'ref_base_2', 'var_base_2', 'position_coverage_2',
#     'A_2', 'C_2', 'T_2', 'G_2']

#     RENAME_COLUMNS = ['con_base', 'ref_base', 'var_base', 'position_coverage',
#                     'A', 'C', 'T', 'G']

#     # Iterate mm levels
#     for mm, cov_arr in mm2overlap.items():

#         # Subset to bases that have coverage in both
#         covs = set(cov_arr)

#         # These represent relevant counts at these posisions
#         if len(SNPtable1) > 0:
#             s1_all = SNPtable1[[(p in covs) for p in SNPtable1['position'].values]].drop_duplicates(
#                         subset=['position'], keep='last')
#             del s1_all['mm']
#             if len(s1_all) == 0:
#                 s1_all = None
#         else:
#             #s1_all = _gen_blank_SNPdb(SNP_COLUMNS)
#             s1_all = None

#         if len(SNPtable2) > 0:
#             s2_all = SNPtable2[[(p in covs) for p in SNPtable2['position'].values]].drop_duplicates(
#                         subset=['position'], keep='last')
#             del s2_all['mm']
#             if len(s2_all) == 0:
#                 s2_all = None
#         else:
#             #s2_all = _gen_blank_SNPdb(SNP_COLUMNS)
#             s2_all = None

#         # Merge
#         if (s1_all is None) & (s2_all is None):
#             #Mdb = _gen_blank_Mdb(OUT_COLUMNS)
#             Mdb = None

#         elif s1_all is None:
#             Mdb = s2_all.rename(columns={c:c + '_2' for c in RENAME_COLUMNS})
#             for c in RENAME_COLUMNS:
#                 Mdb[c + '_1'] = np.nan

#         elif s2_all is None:
#             Mdb = s1_all.rename(columns={c:c + '_1' for c in RENAME_COLUMNS})
#             for c in RENAME_COLUMNS:
#                 Mdb[c + '_2'] = np.nan

#         else:
#             Mdb = pd.merge(s1_all, s2_all, on='position', suffixes=('_1', '_2'), how='outer', copy=False)

#         if Mdb is not None:
#             Mdb.loc[:,'consensus_SNP'] = Mdb.apply(call_con_snps, axis=1)
#             Mdb.loc[:,'population_SNP'] = Mdb.apply(call_pop_snps, axis=1, args=(model_to_use, min_freq))

#             Mdb['mm'] = mm
#             Mdb = Mdb[OUT_COLUMNS]

#             # Only keep SNPs
#             Mdb = Mdb[Mdb['consensus_SNP'] | Mdb['population_SNP'] ]

#             dbs.append(Mdb)

#     if len(dbs) > 0:
#         Mdb = pd.concat(dbs, sort=False)
#     else:
#         Mdb = _gen_blank_Mdb(OUT_COLUMNS)
#     return Mdb

def call_con_snps(row):
    '''
    Call a SNP if the consensus sequnces aren't the same
    '''

    # This was only a SNP in the first sapmle
    if row['con_base_1'] != row['con_base_1']:
        return row['con_base_2'] != row['ref_base_2']

    # This was only a SNP in the second sapmle
    if row['con_base_2'] != row['con_base_2']:
        return row['con_base_1'] != row['ref_base_1']

    # This is a SNP in both samples
    return row['con_base_1'] != row['con_base_2']

def is_present(counts, total, model, min_freq):
    '''
    Return true if the base counts represented by "counts" are detected above background
    '''
    if total in model:
        min_bases = model[total]
    else:
        min_bases = model[-1]

    return (counts >= min_bases) and ((float(counts) / total) >= min_freq)

def call_pop_snps(row, model, min_freq):
    '''
    To be applied to a DataFrame

    Call a SNP if you can't find the consenus of 1 in 2 AND
    you can't find the consensus of 2 in 1 AND
    1 and 2 don't share a minor allele
    '''
    # Are the consensus bases the same?
    if row['con_base_1'] == row['con_base_2']:
        return False

    # Is it a SNP in only one? If so, see if the reference is still there
    if (row['con_base_1'] != row['con_base_1']) | (row['con_base_2'] != row['con_base_2']):

        # In this case, is consensus allele still detected?
        if (row['con_base_1'] != row['con_base_1']):
            count = row['{0}_2'.format(row['ref_base_2'])]
            total = row['position_coverage_2']
            if is_present(count, total, model, min_freq):
                return False

        elif (row['con_base_2'] != row['con_base_2']):
            count = row['{0}_1'.format(row['ref_base_1'])]
            total = row['position_coverage_1']
            if is_present(count, total, model, min_freq):
                return False

        return True

    ### OK, so it's a SNP in both ###

    # Look for con_base_1 in sample 2
    try:
        count = row["{0}_2".format(row['con_base_1'])]
        total = row['position_coverage_2']
        if is_present(count, total, model, min_freq):
            return False
    except:
        print(row)

    # Look for con_base_2 in sample 1
    count = row["{0}_1".format(row['con_base_2'])]
    total = row['position_coverage_1']
    if is_present(count, total, model, min_freq):
        return False

    # Look for minor in both samples
    if 'allele_count_1' in row:
        if (row['allele_count_1'] > 1) & (row['allele_count_2'] > 1):
            if row['var_base_1'] == row['var_base_2']:
                return False

    elif 'morphia_1' in row:
        if (row['morphia_1'] > 1) & (row['morphia_2'] > 1):
            if row['var_base_1'] == row['var_base_2']:
                return False

    return True

def _is_snp(db1, db2, position, p2c1, p2c2, min_freq, model_to_use, debug=False, compare_consensus_bases=False):
    '''
    Determine if the consensus base of db1 is a SNP in db2 and vice versa
    '''
    # Try and out quick
    if p2c1[position] == p2c2[position]:
        return False
    if compare_consensus_bases:
        return True

    # otime = time.time()
    # print("checkpoint {0} {1}".format(0, time.time() - otime))
    # These are sorted by mm above
    dd1 = db1[db1['position'] == position]#.sort_values('mm', ascending=False)#.drop_duplicates(subset='position', keep='last')
    dd2 = db2[db2['position'] == position]#.sort_values('mm', ascending=False)#.drop_duplicates(subset='position', keep='last')

    #print("checkpoint {0} {1}".format(1, time.time() - otime))

    assert len(dd1) > 0, [position, 'one']
    assert len(dd2) > 0, [position, 'two']

    #print("checkpoint {0} {1}".format(2, time.time() - otime))

    con1 = p2c1[position]
    con2 = p2c2[position]

    #print("{0} in {1} is {2}".format(position, p2c1, con1))

    #print("checkpoint {0} {1}".format(3, time.time() - otime))

    # Check if the consensus of db1 is a SNP in db2
    #counts2 = dd2.iloc[-1][C2P[con1]]
    counts2 = dd2.iloc[-1][con1]
    total = dd2.iloc[-1]['position_coverage']

    if total in model_to_use:
        min_bases = model_to_use[total]
    else:
        min_bases = model_to_use[-1]

    if (counts2 >= min_bases) and ((float(counts2) / total) >= min_freq):
        pass
    else:
        return True

    # Check the opposite
    #counts1 = dd1.iloc[-1][C2P[con2]]
    counts1 = dd1.iloc[-1][con2]
    total = dd1.iloc[-1]['position_coverage']

    if total in model_to_use:
        min_bases = model_to_use[total]
    else:
        min_bases = model_to_use[-1]

    if (counts1 >= min_bases) and ((float(counts1) / total) >= min_freq):
        pass
    else:
        return True

    # Looks like you're not a SNP, then!
    if debug:  print('exit2')

    #print("checkpoint {0} {1}".format(5, time.time() - otime))

    return False

def _update_overlap_table(table, scaffold, mm2overlap, mm2coverage, Mdb, name1, name2, mLen):
    '''
    covarage_overlap = the percentage of bases that are either covered or not covered in both
        - So if both scaffolds have 0 coverage, this will be 1
    percent_genome_compared = the percentage of bases in the scaffolds that are covered by both
        - So if both scaffolds have 0 coverave, this will be 0
    compared_bases_count = the number of considered bases
    '''
    rel_mms = set(list(mm2overlap.keys()))
    got_mms = set()

    for mm, mdb in Mdb.groupby('mm'):
        if mm not in rel_mms:
            continue
        got_mms.add(mm)

        overlap = mm2overlap[mm]
        bases = len(overlap)

        table['mm'].append(mm)
        table['scaffold'].append(scaffold)
        table['name1'].append(name1)
        table['name2'].append(name2)
        table['coverage_overlap'].append(mm2coverage[mm])
        table['compared_bases_count'].append(bases)
        table['percent_genome_compared'].append(bases/mLen)
        table['length'].append(mLen)

        snps = len(mdb[mdb['consensus_SNP'] == True])
        popsnps = len(mdb[mdb['population_SNP'] == True])

        table['consensus_SNPs'].append(snps)
        table['population_SNPs'].append(popsnps)

        if bases == 0:
            table['popANI'].append(np.nan)
            table['conANI'].append(np.nan)
        else:
            table['conANI'].append((bases - snps) / bases)
            table['popANI'].append((bases - popsnps) / bases)

    # Doing this to allow the group-by
    for mm in rel_mms - got_mms:
        overlap = mm2overlap[mm]
        bases = len(overlap)

        table['mm'].append(mm)
        table['scaffold'].append(scaffold)
        table['name1'].append(name1)
        table['name2'].append(name2)
        table['coverage_overlap'].append(mm2coverage[mm])
        table['compared_bases_count'].append(bases)
        table['percent_genome_compared'].append(bases / mLen)
        table['length'].append(mLen)

        snps = 0
        popsnps = 0

        table['consensus_SNPs'].append(snps)
        table['population_SNPs'].append(popsnps)

        if bases == 0:
            table['popANI'].append(np.nan)
            table['conANI'].append(np.nan)
        else:
            table['conANI'].append((bases - snps) / bases)
            table['popANI'].append((bases - popsnps) / bases)

    return table