'''
This contains a number of methods and functions to calculate the iRep metric

https://github.com/christophertbrown/iRep

https://www.nature.com/articles/nbt.3704
'''

import os
import sys
import glob
import scipy
import lmfit
import scipy.signal
import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from collections import defaultdict

def calculate_iRep_from_coverage_array(rcov, num_contigs, gc_windows=None):
    '''
    Argumemts:
        rcov: genome-wide array of coverage values
        num_contigs: number of contigs in genome
        gc_windows: GC content of windows to do GC coverage correction

    Returns:
        iRep: iRep value or np.nan if filtered out
        iRep_junk: Dictionary of other iRep outputs
    '''
    FilterChriteria = {'kept_windows':np.nan, # Calculated in coverage_windows
                       'avg_cov':np.nan, # Calculated on raw array
                       'r2':np.nan, # Calculated from iRep itself
                       'fragMbp':np.nan } # Calculated from iRep itself

    # Filter chriteria
    length = len(rcov)
    FilterChriteria['avg_cov'] = np.mean(rcov)
    FilterChriteria['fragMbp'] = num_contigs/(float(length)/1000000)

    # Calculate the windows
    oIdb = _iRep_windows(rcov)

    # Add GC if appropriate
    if gc_windows is not None:
        oIdb = pd.merge(oIdb, gc_windows, on='index')

    # Filter out junk windows
    Idb = _iRep_filter_windows(oIdb, on='coverage')
    FilterChriteria['kept_windows'] = len(Idb) / len(oIdb)

    # Get raw iRep values
    Idb.loc[:,'coverage_OLT'] = _iRep_log_transform(Idb['coverage'])
    iRep = _calc_iRep(Idb, length, on='coverage_OLT', FilterChriteria=FilterChriteria)
    FilterChriteria['unfiltered_raw_iRep'] = iRep

    # Get gc-corrected iRep values
    FilterChriteria['iRep_GC_corrected'] = False
    if gc_windows is not None:
        Idb = _iRep_gc_bias(Idb)
        Idb.loc[:,'coverage_LT'] = _iRep_log_transform(Idb['corrected_coverage'])
        iRep = _calc_iRep(Idb, length, on='coverage_LT')
        FilterChriteria['unfiltered_iRep'] = iRep
        FilterChriteria['iRep_GC_corrected'] = True

    # Get raw iRep values
    Idb.loc[:,'coverage_OLT'] = _iRep_log_transform(Idb['coverage'])
    iRep = _calc_iRep(Idb, length, on='coverage_OLT', FilterChriteria=FilterChriteria)
    FilterChriteria['unfiltered_raw_iRep'] = iRep

    # See if iRep passes
    if (FilterChriteria['kept_windows'] < 0.98) or \
        (FilterChriteria['avg_cov'] < 5) or \
        (FilterChriteria['r2'] < 0.9) or \
        (FilterChriteria['fragMbp'] > 175):

        iRep = np.nan

    return iRep, FilterChriteria

def generate_gc_windows(order, scaff2sequence, mask_edges=100):
    '''
    Calculate the GC content for windows across a .fasta file

    Argumnets:
        order = list of scaffolds in the order they should be in
        scaff2sequence = scaffold2sequence (scaff2sequence = SeqIO.to_dict(SeqIO.parse(fasta_loc, "fasta")))
        mask_edges = remove this many bp from start and end of contigs

    Modified for speed
    '''
    # Load the .fasta
    #scaff2sequence = SeqIO.to_dict(SeqIO.parse(file_loc, "fasta"))

    # Make the genome sequence into a list
    splits = []
    for scaff in order:
        seq = scaff2sequence[scaff]
        if mask_edges:
            seq = seq[100:len(seq)-100]
            splits.append(seq)
            #genome_seq.extend(seq)
    genome_seq = sum(splits, [])

    # Calculate GC content
    gcdb = _iRep_gc_content(genome_seq)

    return gcdb

def _iRep_gc_content(seq, window = 5000, slide = 100):
    """
    iRep gc_content message

    calculate gc content over sequence windows
    """
    # convert GC
    replacements = {'G':1, 'C':1, 'A':0, 'T':0, 'N':0}
    GC = [] # G - C
    for base in seq:
        try:
            GC.append(replacements[base.upper()])
        except:
            GC.append(0)
    # calculate gc content over sliding windows
    i = 0
    weights = np.ones(window)
    table = defaultdict(list)
    for gc in scipy.signal.fftconvolve(GC, weights, 'valid').tolist()[0::slide]:
        table['index'].append(i)
        table['GC_content'].append(gc/window)
        i += slide
    return pd.DataFrame(table)

def _iRep_windows(cov, window=5000, slide=100, FilterChriteria=None):
    '''
    Replicates *coverage_windows()*
    '''
    table = defaultdict(list)

    i = 0
    weights = np.ones(window)
    for c in scipy.signal.fftconvolve(cov, weights, 'valid').tolist()[0::slide]:
        table['index'].append(i)
        table['coverage'].append(c/window)
        i += slide

    return pd.DataFrame(table)

def _iRep_filter_windows(cov, on='coverage', mdif = float(8)):
    '''
    Replicates *filter_windows()*

    Remove windows with weird coverage

    Edited to improve speed in version 1.3.0l
    '''
    med = np.median(cov[on])

    return cov[[True if \
            ((y > 0) and (med > 0) and
            (abs(float(max([y, med])) / float(min([y, med]))) <= mdif))
            else False for y in cov[on]]]
    # keeper_index = []
    # for i, row in cov.iterrows():
    #     y = row[on]
    #     if y <= 0 or med <= 0:
    #         continue
    #     if abs(float(max([y, med])) / float(min([y, med]))) > mdif:
    #         continue
    #     keeper_index.append(row['index'])
    #
    # return cov[cov['index'].isin(keeper_index)]

def _iRep_log_transform(array):
    lt = []
    eps = 1e-50
    for i in array:
        if i < eps:
            lt.append(np.log2(eps))
        else:
            lt.append(np.log2(i))
    return lt

def _calc_iRep(db, length, on='coverage_OLT', FilterChriteria=None):
    '''
    Replicates iRep_from_windows
    '''
    Ys = sorted(list(db[on]))
    windows = len(Ys)

    if windows == 0:
        return np.nan

    dif = float(length)/float(windows)
    Xs = [int(i * dif) + 1 for i, value in enumerate(Ys, 0)]
    Xt, Yt = trim_data((Xs, Ys), xy = True)
    db = pd.DataFrame({'index':Xt, 'cov':Yt})
    m, b, fit, r2, info = fit_coverage((Xt, Yt, None, True))
    iRep = 2**(m * length)

    if FilterChriteria is not None:
        FilterChriteria['r2'] = r2

    return iRep

def trim_data(data, xy, p = 0.1):
    """
    RIGHT FROM iREP

    remove data from ends of sorted list
    """
    if xy is False:
        length = len(data)
        num = int(length * (p/2))
        return data[num:length - num]
    X, Y = data
    length = len(X)
    num = int(length * (p/2))
    return X[num:length - num], Y[num:length - num]

def fit_coverage(pars):
    """
    RIGHT FROM iREP

    fit line to sorted coverage values to get slope
    """
    x, y, info, return_fit = pars
    if len(x) <= 2: # make sure there is more data than parameters
        if return_fit is False:
            return (False, False, False, info)
        else:
            return (False, False, False, False, info)
    # Parameters
    Pars = lmfit.Parameters()
    ## y = mx + b
    Pars.add('m', value = 1, vary = True)
    Pars.add('b', value = 1)
    # fit data to model
    mi = lmfit.minimize(coverage_function, Pars, args = (x,), \
            kws = {'data':y, 'printPs':False}, method = 'leastsq')
    # calculate r-squared
    r2 = 1 - (mi.residual.var() / np.var(y))
    if return_fit is False:
        return (mi.params['m'].value, mi.params['b'].value, r2, info)
    # get fitted values
    fit = [x, coverage_function(mi.params, x)]
    return (mi.params['m'].value, mi.params['b'].value, fit, r2, info)

def coverage_function(pars, X, data = None, printPs = False):
    """
    RIGHT FROM iREP

    linear function for sorted coverage profile
    y = mx + b
    """
    m = pars['m'].value
    b = pars['b'].value
    if printPs is True:
        print('m: %s b: %s' % \
            ('{:,}'.format(int(m)), '{:,}'.format(int(b))))
    results = [float(m * x) + b for x in X]
    if data is None:
        return np.asarray(results)
    return np.asarray([y - data[i] for i, y in enumerate(results)]) # model - data

def _iRep_gc_bias(Idb, correction_threshold=0.0):
    '''
    iRep gc_bias method
    '''
    # fit line
    m, b, fit, r2, l = fit_coverage((Idb['GC_content'].tolist(), Idb['coverage'].tolist(), False, True))

    # remove outliers
    Idb.loc[:,'error'] = [abs(cov - (m * gc + b)) for gc, cov in zip(Idb['GC_content'], Idb['coverage'])]
    try:
        cutoff = sorted(Idb['error'].tolist(), reverse = True)[int(len(Idb['error'])*0.01)]
    except:
        cutoff = 0
    FIdb = Idb[~(Idb['error'] >= cutoff)]

    # re-fit with filtered data
    m, b, fit, r2, l = fit_coverage((FIdb['GC_content'].tolist(), FIdb['coverage'].tolist(), False, True))
    if r2 < correction_threshold:
        Idb['corrected_coverage'] = Idb['coverage']
        return Idb

    # correct coverge
    corrected = [[], []]
    av = np.average(Idb['coverage'])
    Idb['corrected_coverage'] = [cov + (av - (m * gc + b)) for cov, gc in zip(Idb['coverage'], Idb['GC_content'])]

    return Idb
