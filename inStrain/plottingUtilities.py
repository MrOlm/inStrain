#!/usr/bin/env python

import os
import sys
import h5py
import logging
import traceback
import warnings
import numpy as np
import scipy.cluster.hierarchy
import scipy.spatial.distance as ssd
from collections import defaultdict

import inStrain.SNVprofile
import inStrain.readComparer

import inStrain.profile.profile_utilities

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import seaborn as sns

import drep.d_cluster
import drep.d_analyze

matplotlib.rcParams['pdf.fonttype'] = 42

def mm_plot(db, left_val='breadth', right_val='coverage', title='',\
           minANI=0.9):
    '''
    The input db for this is "mm_genome_info" from "makeGenomeWide" in genomeUtilities.py
    '''
    db = db.sort_values('ANI_level')
    sns.set_style('white')

    # breadth
    fig, ax1 = plt.subplots()
    ax1.plot(db['ANI_level'], db[left_val], ls='-', color='blue')
    if left_val == 'breadth':
        ax1.plot(db['ANI_level'], estimate_breadth(db['coverage']), ls='--', color='lightblue')
    ax1.set_ylabel(left_val, color='blue')
    ax1.set_xlabel('Minimum read ANI level')
    ax1.set_ylim(0,1)

    # coverage
    ax2 = ax1.twinx()
    ax2.plot(db['ANI_level'], db[right_val], ls='-', color='red')
    ax2.set_ylabel(right_val, color='red')
    ax2.set_ylim(0,)

    # asthetics
    plt.xlim(1, max(minANI, db['ANI_level'].min()))

    plt.title(title)

def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage

    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1

# def genome_wide_plot(IS_locs, scaffolds, what='coverage', ANI_levels=[100, 98, 0], window_len=1000):
#     '''
#     Arguments:
#         IS_locs = list of IS objects
#         scaffolds = list of scaffolds to profile and plot (in order)

#     Keyword arguments:
#         ANI_levels = list of ANI levesl to plot
#         window_len = length of each window to profile
#     '''
#     if what == 'coverage':
#         item = 'covT'
#     elif what == 'clonality':
#         item = 'clonT'

#     # Load coverages for the scaffolds at each ANI level
#     dbs = []
#     for IS_loc in IS_locs:
#         IS = inStrain.SNVprofile.SNVprofile(IS_loc)
#         if what in ['coverage', 'clonality']:
#             wdb, breaks = load_windowed_coverage(IS, scaffolds, window_len=window_len, ANI_levels=ANI_levels, item=item)
#         elif what in ['linkage']:
#             wdb, breaks = load_windowed_linkage(IS, scaffolds, window_len=window_len, ANI_levels=ANI_levels)
#         elif what in ['snp_density']:
#             wdb, breaks = load_windowed_SNP_density(IS, scaffolds, window_len=window_len, ANI_levels=ANI_levels)

#         wdb['name'] = os.path.basename(IS_loc)
#         dbs.append(wdb)
#     Wdb = pd.concat(dbs, sort=True)

#     # Make the plot
#     multiple_coverage_plot(Wdb, breaks, thing=what)

#     return Wdb, breaks

def load_windowed_metrics(scaffolds, s2l, rLen, metrics=None, window_len=None, ANI_levels=[0, 100],
                        min_scaff_len=0, report_midpoints=False, covTs=False, clonTs=False,
                        raw_linkage_table=False, cumulative_snv_table=False):

    if metrics is None:
        metrics = ['coverage', 'nucl_diversity', 'linkage', 'snp_density']

    if type(metrics) != type([]):
        print("Metrics must be a list")
        return

    # Figure out the MMs needed
    #rLen = IS.get_read_length()
    mms = [_get_mm(None, ANI, rLen=rLen) for ANI in ANI_levels]

    # Sort the scaffolds
    #s2l = IS.get('scaffold2length')
    scaffolds = sorted(scaffolds, key=s2l.get, reverse=True)

    if min_scaff_len > 0:
        scaffolds = [s for s in scaffolds if s2l[s] >= min_scaff_len]

    # Figure out the window length
    if window_len == None:
        window_len = int(sum([s2l[s] for s in scaffolds]) / 100)
    else:
        window_len = int(window_len)

    # Calculate the breaks
    breaks = []
    midpoints = {}
    tally = 0
    for scaffold in scaffolds:
        midpoints[scaffold] = tally + int(s2l[scaffold] / 2)
        tally += s2l[scaffold]
        breaks.append(tally)


    dbs = []
    if 'coverage' in metrics:
        if covTs == False:
            logging.error("need covTs for coverage")
            raise Exception
        cdb = load_windowed_coverage_or_clonality('coverage', covTs, scaffolds, window_len, mms, ANI_levels, s2l)
        cdb['metric'] = 'coverage'
        dbs.append(cdb)
#     if 'clonality' in metrics:
#         cdb = load_windowed_coverage_or_clonality(IS, 'clonality', scaffolds, window_len, mms, ANI_levels, s2l)
#         cdb['metric'] = 'clonality'
#         dbs.append(cdb)
    if 'nucl_diversity' in metrics:
        if clonTs == False:
            logging.error("need clonTs for microdiversity")
            raise Exception
        cdb = load_windowed_coverage_or_clonality('nucl_diversity', clonTs, scaffolds, window_len, mms, ANI_levels, s2l)
        cdb['metric'] = 'nucl_diversity'
        dbs.append(cdb)
    if 'linkage' in metrics:
        if raw_linkage_table is False:
            logging.error("need raw_linkage_table for linkage")
            raise Exception
        cdb = load_windowed_linkage(raw_linkage_table, scaffolds, window_len, mms, ANI_levels, s2l)
        cdb['metric'] = 'linkage'
        dbs.append(cdb)
    if 'snp_density' in metrics:
        if cumulative_snv_table is False:
            logging.error("need cumulative_snv_table for snp_density")
            raise Exception
        if len(cumulative_snv_table) > 0:
            cdb = load_windowed_SNP_density(cumulative_snv_table, scaffolds, window_len, mms, ANI_levels, s2l)
            cdb['metric'] = 'snp_density'
            dbs.append(cdb)

    if len(dbs) > 0:
        Wdb = pd.concat(dbs, sort=True)
        Wdb = Wdb.rename(columns={'avg_cov':'value'})
    else:
        Wdb = pd.DataFrame()

    # Add blanks at the breaks
    table = defaultdict(list)
    for mm, ani in zip(mms, ANI_levels):
        for metric in Wdb['metric'].unique():
            for bre in breaks:
                table['scaffold'].append('break')
                table['mm'].append(mm)
                table['ANI'].append(ani)
                table['adjusted_start'].append(bre) # The minus one makes sure it doenst split things it shouldnt
                table['adjusted_end'].append(bre)
                table['value'].append(np.nan)
                table['metric'].append(metric)
    bdb = pd.DataFrame(table)
    Wdb = pd.concat([Wdb, bdb], sort=False)
    if len(Wdb) > 0:
        Wdb['midpoint'] = [np.mean([x, y]) for x, y in zip(Wdb['adjusted_start'], Wdb['adjusted_end'])]
        Wdb = Wdb.sort_values(['metric', 'mm', 'midpoint', 'scaffold'])

    if report_midpoints:
        return Wdb, breaks, midpoints
    else:
        return Wdb, breaks

def load_windowed_coverage_or_clonality(thing, covTs, scaffolds, window_len, mms, ANI_levels, s2l):
    '''
    Get the windowed coverage

    Pass in a clonTs for microdiversity and covTs for coverage
    '''
    if thing == 'coverage':
        item = 'covT'
    elif thing == 'nucl_diversity':
        item = 'clonT'
    else:
        print("idk what {0} is".format(thing))
        return

    # Get the covTs
    #covTs = IS.get(item, scaffolds=scaffolds)

    # Make the windows
    dbs = []
    tally = 0
    breaks = []
    for scaffold in scaffolds:
        if scaffold not in covTs:
            tally += s2l[scaffold]
            breaks.append(tally)
            continue
        else:
            covT = covTs[scaffold]

        for mm, ani in zip(mms, ANI_levels):
            if item == 'covT':
                cov = inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(covT, mm)
                if len(cov) == 0:
                    continue
                db = _gen_windowed_cov(cov, window_len, sLen=s2l[scaffold])

            elif item == 'clonT':
                cov = _get_basewise_clons3(covT, mm)
                if len(cov) == 0:
                    continue
                db = _gen_windowed_cov(cov, window_len, sLen=s2l[scaffold], full_len=False)
                db['avg_cov'] = [1 - x if x == x else x for x in db['avg_cov']]

            db['scaffold'] = scaffold
            db['mm'] = mm
            db['ANI'] = ani

            db['adjusted_start'] = db['start'] + tally
            db['adjusted_end'] = db['end'] + tally

            dbs.append(db)
        tally += s2l[scaffold]
        breaks.append(tally)

    if len(dbs) > 0:
        Wdb = pd.concat(dbs)
    else:
        Wdb = pd.DataFrame()
    return Wdb#, breaks

def load_windowed_linkage(Ldb, scaffolds, window_len, mms, ANI_levels, s2l, on='r2'):

    # Get the linkage table
    #Ldb = IS.get('raw_linkage_table')
    Ldb = Ldb[Ldb['scaffold'].isin(scaffolds)].sort_values('mm')
    got_scaffolds = set(Ldb['scaffold'].unique())

    # Make the windows
    dbs = []
    tally = 0
    breaks = []
    for scaffold in scaffolds:
        if scaffold not in got_scaffolds:
            tally += s2l[scaffold]
            breaks.append(tally)
            continue
        else:
            ldb = Ldb[Ldb['scaffold'] == scaffold]

        for mm, ani in zip(mms, ANI_levels):
            db = ldb[ldb['mm'] <= int(mm)].drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')
            cov = db.set_index('position_A')[on].sort_index()

            db = _gen_windowed_cov(cov, window_len, sLen=s2l[scaffold], full_len=False)
            db['scaffold'] = scaffold
            db['mm'] = mm
            db['ANI'] = ani

            db['adjusted_start'] = db['start'] + tally
            db['adjusted_end'] = db['end'] + tally

            dbs.append(db)
        tally += s2l[scaffold]
        breaks.append(tally)

    if len(dbs) > 0:
        Wdb = pd.concat(dbs)
    else:
        Wdb = pd.DataFrame()
    return Wdb

def load_windowed_SNP_density(Ldb, scaffolds, window_len, mms, ANI_levels, s2l):
    # Get the table
    #Ldb = IS.get('cumulative_snv_table')
    Ldb = Ldb[Ldb['scaffold'].isin(scaffolds)].sort_values('mm')

    got_scaffolds = list(Ldb['scaffold'].unique())

    # Make the windows
    dbs = []
    tally = 0
    breaks = []
    for scaffold in scaffolds:
        if scaffold not in got_scaffolds:
            tally += s2l[scaffold]
            breaks.append(tally)
            continue
        else:
            ldb = Ldb[Ldb['scaffold'] == scaffold]

        for mm, ani in zip(mms, ANI_levels):
            db = ldb[ldb['mm'] <= int(mm)].drop_duplicates(subset=['scaffold', 'position'], keep='last')
            cov = db.set_index('position')['ref_base'].sort_index()

            db = _gen_windowed_cov(cov, window_len, sLen=s2l[scaffold], full_len='count')
            db['scaffold'] = scaffold
            db['mm'] = mm
            db['ANI'] = ani

            db['adjusted_start'] = db['start'] + tally
            db['adjusted_end'] = db['end'] + tally

            dbs.append(db)
        tally += s2l[scaffold]
        breaks.append(tally)

    if len(dbs) > 0:
        Wdb = pd.concat(dbs)
    else:
        Wdb = pd.DataFrame()
    return Wdb

# def multiple_coverage_plot(Wdb, breaks, thing='coverage'):
#     '''
#     Make the multiple coverage plot
#     '''
#     # Get set up for multiple rows
#     i = len(Wdb['name'].unique())
#     if i > 1:
#         fig, ax = plt.subplots(i, 1, sharex=True)
#     else:
#         ax = {}
#         ax[0] = plt.gca()

#     i = 0
#     for name, wdb in Wdb.groupby('name'):
#         med = wdb['avg_cov'].median()

#         # Rotate colors:
#         colors = ['red', 'blue', 'black']
#         c = 0
#         for mm, ddb in wdb.groupby('mm'):
#             ax[i].plot(ddb['adjusted_start'], ddb['avg_cov'], c=colors[c], label=mm)
#             c += 1

#         # Set stuff up
#         if thing == 'coverage':
#             ax[i].set_ylim([0,med*2])
#         ax[i].set_title("{0}".format(name))
#         ax[i].grid(False)

#         if i == 0:
#             ax[i].legend(loc='upper left', title='mm level')

#         # Add breaks
#         for b in breaks:
#             ax[i].axvline(b, ls='-', c='lightgrey', zorder=-1)

#         i += 1

#     plt.gcf().set_size_inches(12, 2*len(Wdb['name'].unique()))
#     plt.xlabel('genome position')
#     #plt.ylabel('coverage')
#     plt.xlim(0, Wdb['adjusted_start'].max())

def _get_mm(IS, ANI, rLen = None):
    '''
    Get the mm corresponding to an ANI level in an IS
    '''
    if ANI > 1:
        ANI = ANI / 100

    if rLen == None:
        rLen = IS.get_read_length()
        #rLen = IS.get('mapping_info')['mean_pair_length'].tolist()[0]
    mm = int(round((rLen - (rLen * ANI))))
    return mm

def _gen_windowed_cov(cov, window_len, sLen=False, full_len=True):
    '''
    From a series of coverage values, return windows
    '''
    if sLen == False:
        sLen = cov.index.max()

    table = defaultdict(list)
    i = 0
    numberChunks = sLen // window_len + 1
    for cov in iterate_chunks_series(cov, chunkSize=window_len, sLen=sLen):
        if i + 1 == numberChunks:
            mLen = sLen - (i * window_len)
            if mLen == 0:
                continue
        else:
            mLen = window_len

        table['start'].append(i*window_len)
        table['end'].append(i*window_len + mLen)

        if full_len == True:
            table['avg_cov'].append(cov.sum() / mLen)
        elif full_len == False:
            table['avg_cov'].append(cov.mean())
        elif full_len == 'count':
            table['avg_cov'].append(len(cov)/mLen)
        i += 1

    return pd.DataFrame(table)

def _get_basewise_clons3(clonT, MM, fill_zeros=False):
    p2c = {}
    mms = sorted([int(mm) for mm in list(clonT.keys()) if int(mm) <= int(MM)])
    for mm in mms:
        p2c.update(clonT[mm].to_dict())

    inds = []
    vals = []
    for ind in sorted(p2c.keys()):
        inds.append(ind)
        vals.append(p2c[ind])

    counts = pd.Series(data = vals, index = np.array(inds).astype('int'))

    if fill_zeros:
        counts = counts.append(pd.Series(np.zeros(fill_zeros - len(counts))))

    return counts

def iterate_chunks_series(d, chunkSize=100, sLen=False):
    '''
    Break up Ndbs into chunks
    '''
    if sLen == False:
        sLen = d.index.max()

    numberChunks = sLen // chunkSize + 1
    for i in range(numberChunks):
        #print(i, int(i*chunkSize), int((i+1)*chunkSize), len(d))
        start = int(i*chunkSize)
        end = int((i+1)*chunkSize)
        yield (d.loc[start:end])

def plot_RC_SNPs(hd5_loc, scaffolds, s2l, mm_level=None, samples=None):
    '''
    Make an all-vs-all comparison here
    '''
    # Load the object
    pairs = []
    for name1 in samples:
        for name2 in samples:
            pair = '-vs-'.join(sorted([name1, name2]))
#             if pair == 'typeStrainsAlpha_v1.fasta-vs-N5_246_000G1.sorted.bam-vs-typeStrainsAlpha_v1.fasta-vs-N5_247_000G1.sorted.bam':
#                 print("got it!")
            pairs.append(pair)
    scaff2pair2mm2SNPs = load_scaff2pair2mm2SNPs(hd5_loc, scaffolds=scaffolds,
                        pairs=pairs)

    # Get set up for multiple rows
    rows = (len(samples) * len(samples)) + (len(samples) - 1)
    assert rows > 1
    fig, ax = plt.subplots(rows, 1, sharex=True)

    # Sort the scaffolds
    scaffolds = [w for w in sorted(s2l, key=s2l.get, reverse=True) if w in scaffolds]

    # Do the plots
    i = 0
    for sample1 in samples:
        for sample2 in samples:
            #print("plotting {0} vs {1}".format(sample1, sample2))
            fdb, breaks = prepare_rugplot(scaff2pair2mm2SNPs, scaffolds, '-vs-'.join(sorted([sample1, sample2])), s2l, mm_level=mm_level)
            _plot_pair(fdb, breaks, ax[i], '-vs-'.join([sample1, sample2]))
            i += 1

        try:
            ax[i].grid(False)
            ax[i].axis('off')
            i += 1
        except:
            pass

    plt.xlim(0, breaks[-1])
    plt.gcf().set_size_inches(12, 1*rows)

def prepare_rugplot(scaff2pair2mm2SNPs, scaffolds, pair, s2l, mm_level=None):
    table = defaultdict(list)
    breaks = []
    adjusted_loc = 0

    for scaffold in scaffolds:
        if scaffold in scaff2pair2mm2SNPs:
            if pair in scaff2pair2mm2SNPs[scaffold]:
                snps = scaff2pair2mm2SNPs[scaffold][pair][_get_mm(mm_level, scaff2pair2mm2SNPs[scaffold][pair])]
            else:
                snps = []
        else:
            snps = []

#         print("{0} - {1} snps".format(scaffold, len(snps)))
        for snp in snps:
            table['loc'].append(snp)
            table['adjusted_loc'].append(snp + adjusted_loc)
            table['scaffold'].append(scaffold)
        adjusted_loc += s2l[scaffold]
        breaks.append(adjusted_loc)

    return pd.DataFrame(table), breaks

# def _get_mm(mm_level, mm2SNPs):
#     if mm_level != None:
#         return mm_level
#     else:
#         return max([int(x) for x in list(mm2SNPs.keys())])

def _plot_pair(fdb, breaks, ax, pair):
    if len(fdb) > 0:
        sample = fdb['adjusted_loc']
        ax.plot(sample, [0.01]*len(sample), '|', color='k')

    # Add breaks
    for b in breaks:
        ax.axvline(b, ls='-', c='lightgrey', zorder=-1)
        #sns.distplot(fdb['adjusted_loc'], ax=ax)
    ax.axhline(0.01, ls='-', c='lightgrey', zorder=-1)

    # Adjust
    ax.set_title("{0}".format(pair))
    ax.set_ylim(0, 0.02)
    ax.grid(False)
    ax.axis('off')

def load_scaff2pair2mm2SNPs(location, scaffolds=[], pairs=[]):
    scaff2pair2mm2SNPs = {}
    f = h5py.File(location, 'r')

    for thing in list(f.keys()):
        scaff, pair, mm = thing.split('::')

        if scaffolds != []:
            if scaff not in scaffolds:
                continue

        if pairs != []:
            if pair not in pairs:
                continue

        dset = list(f[thing])
        mm = int(mm)

        if scaff not in scaff2pair2mm2SNPs:
            scaff2pair2mm2SNPs[scaff] = {}

        if pair not in scaff2pair2mm2SNPs[scaff]:
            scaff2pair2mm2SNPs[scaff][pair] = {}

        scaff2pair2mm2SNPs[scaff][pair][mm] = dset # convert from 2d array to series

    return scaff2pair2mm2SNPs


def genomeWide_microdiveristy_metrics_plot(Wdb, breaks, title=''):
    '''
    Make the multiple metrics plot
    '''
    # Get set up for multiple rows
    i = len(Wdb['metric'].unique())
    if i > 1:
        fig, ax = plt.subplots(i, 1, sharex=True)
    else:
        ax = {}
        ax[0] = plt.gca()

    i = 0
    for metric in ['linkage', 'snp_density', 'coverage', 'nucl_diversity']:
    #for metric, wdb in Wdb.groupby('metric'):
        if metric not in set(Wdb['metric'].tolist()):
            continue

        wdb = Wdb[Wdb['metric'] == metric]
        med = wdb['value'].median()

        # Rotate colors:
        colors = ['red', 'blue', 'black']
        c = 0
        for mm, ddb in wdb.groupby('ANI'):
            ax[i].plot(ddb['midpoint'], ddb['value'], c=colors[c], label=mm, marker='o', ms=1)#, ls='')
            c += 1

        ax[i].set_title("{0}".format(metric))
        ax[i].grid(False)

        if i == 0:
            ax[i].legend(loc='upper left', title='Min read ANI (%)')

        # Add breaks
        for b in breaks:
            ax[i].axvline(b, ls='-', c='lightgrey', zorder=-1)

        i += 1

    plt.xlabel('genome position')
    plt.xlim(0, Wdb['midpoint'].max())
    plt.suptitle(title, y=0.999)
    plt.subplots_adjust(hspace=0.3)

def read_ani_dist_plot(db, title=None):
    # Plot it
    plt.plot(db['ANI_level'], db['read_pairs'])

    # Adjust x-tickets
    if db['ANI_level'].max() != db['ANI_level'].min():
        plt.gca().set_xlim(db['ANI_level'].max(), db['ANI_level'].min())

    # Other stuff
    rl = int(db['read_length'].tolist()[0])
    plt.xlabel('Read ANI level')
    plt.ylabel("Numbner of read pairs (average length {0}bp)".format(rl))
    plt.title(title)


def prepare_read_ani_dist_plot(IS):
    # Make a non-cumulative scaffold table
    covTS = IS.get('covT')
    s2l = IS.get('scaffold2length')
    table = defaultdict(list)
    for scaffold, covT in covTS.items():
        for mm, counts in covT.items():
            lengt = s2l[scaffold]
            counts = counts.append(pd.Series(np.zeros(lengt - len(counts))))

            coverage = np.mean(counts)

            table['scaffold'].append(scaffold)
            table['mm'].append(mm)
            table['coverage'].append(coverage)
            table['length'].append(s2l[scaffold])

    # Make it genome wide
    db = pd.DataFrame(table)
    stb = IS.get('scaffold2bin')
    b2l = IS.get('bin2length')
    gdb = inStrain.genomeUtilities._add_stb(db, stb)

    table = defaultdict(list)
    for mm, Odb in gdb.groupby('mm'):
        for genome, df in Odb.groupby('genome'):
            table['mm'].append(mm)
            table['genome'].append(genome)
            table['length'].append(int(b2l[genome]))

            for col in ['coverage']:
                table[col].append(sum(x * y for x, y in zip(df[col], df['length'])) / b2l[genome])
    db = pd.DataFrame(table)

    # Add the number of read-pairs
    #readLen = int(IS.get('mapping_info')['mean_pair_length'].tolist()[0])
    readLen = int(IS.get_read_length())
    db['read_length'] = readLen
    db['mm'] = db['mm'].astype(int)
    db['read_pairs'] = [int((x*y) / (readLen * 2)) for x, y in zip(db['coverage'], db['length'])]
    db['ANI_level'] = [(readLen - mm)/readLen for mm in db['mm']]

    return db

def major_allele_freq_plot(db, title=''):
    db['major_allele_freq'] = [max(x, y) for x, y in zip(db['var_freq'], db['ref_freq'])]
    sns.distplot(db['major_allele_freq'], bins=np.arange(0.5, 1, 0.01), kde=False)

    plt.xlim(0.5, 1)
    plt.title(title)
    plt.xlabel('Major allele frequency')
    plt.ylabel('Number of SNPs')

def linkage_decay_plot(db, chunkSize=5, min_vals=5, title=''):
    COLS = ['r2', 'r2_normalized', 'd_prime', 'd_prime_normalized']

    # Make chunks
    max_d = db['distance'].max()
    table = defaultdict(list)
    numberChunks = max_d // chunkSize + 1
    db['distance'] = db['distance'].astype(int)
    for i in range(numberChunks):
        d = db[(db['distance'] >= int(i*chunkSize)) & (db['distance'] < int((i+1)*chunkSize))]

        table['d_start'].append(int(i*chunkSize))
        table['d_end'].append(int((i+1)*chunkSize))
        table['values'].append(len(d))
        for col in COLS:
            table[col].append(d[col].mean())
            table[col + '_values'].append(len(d[~d[col].isna()]))
    df = pd.DataFrame(table)
    df['distance'] = [np.mean([x, y]) for x, y in zip(df['d_start'], df['d_end'])]
    df = df.sort_values('distance')

    for col in COLS:
        df[col] = [c if v >= min_vals else np.nan for c, v in zip(df[col], df[col + '_values'])]

    for col in COLS:
        sns.lineplot(df['distance'], df[col], label=col, marker='o')

    plt.title(title)
    plt.xlabel('Distance between SNPs (bp)\nAveraged over {0}bp windows; plotting windows with at least {1} values'.format(
                chunkSize, min_vals))
    plt.ylabel("SNP linkage")
    return df

def read_filtering_plot(db, title=''):
    # Prepare data
    keep_cols = [x for x in db.columns if 'pass' in x] \
                 + ['unfiltered_reads', 'unfiltered_pairs', 'filtered_pairs']
    db = db.melt(id_vars=['genome'], value_vars=keep_cols)

    # Rename
    c2c = {'unfiltered_reads':'Total mapping reads (divided by 2)',
          'unfiltered_pairs':'Total mapped pairs',
          'pass_min_mapq':'Pairs passing mapQ threshold',
          'pass_max_insert':'Pairs passing max insert size threshold',
          'pass_min_insert':'Pairs passing min insert size threshold',
          'pass_filter_cutoff':'Pairs passing ANI threshold',
          'filtered_pairs':'Total filtered pairs'}
    db['variable'] = [c2c[x] if x in c2c else x for x in db['variable']]
    db['value'] = [int(x/2) if y == 'Total mapping reads (divided by 2)' else x for x, y in zip(
                db['value'], db['variable'])]

    # Set up colors
    v2c = {v:'grey' for v in db['variable'].unique()}
    v2c['Total filtered pairs'] = 'green'

    db = db.sort_values(['value', 'variable'], ascending=False)
    ax = sns.barplot(db['value'], db['variable'], palette=v2c)
    plt.xlabel("Number of read pairs")
    plt.ylabel("")

    # Annotate every single Bar with its value, based on it's width
    offset = db['value'].max() / 12
    total = db[db['variable'] == 'Total mapped pairs']['value'].tolist()[0]
    if total > 0:
        for i, p in enumerate(ax.patches):
            if i == 0:
                continue
            width = p.get_width()
            plt.text(offset + p.get_width(), p.get_y()+0.55*p.get_height(),
                     '{:1.0f}%'.format((width/total)*100),
                     ha='center', va='center')

    plt.title(title)

def scaffold_inspection_plot(Wdb, breaks, midpoints, title=''):
    '''
    Make the multiple metrics plot
    '''
    sns.set_style('whitegrid')

    # Get set up for multiple rows
    i = len(Wdb['metric'].unique())
    if i > 1:
        fig, ax = plt.subplots(1, i, sharey=True)
    else:
        ax = {}
        ax[0] = plt.gca()

    i = 0
    for metric in ['linkage', 'snp_density', 'coverage', 'nucl_diversity']:
    #for metric, wdb in Wdb.groupby('metric'):
        if metric not in set(Wdb['metric'].tolist()):
            continue

        wdb = Wdb[Wdb['metric'] == metric]

        # Rotate colors:
        colors = ['red', 'blue', 'black']
        c = 0
        for mm, ddb in wdb.groupby('ANI'):
            ax[i].plot(ddb['value'], ddb['midpoint'], c=colors[c], label=mm, marker='o', ms=5)#, ls='')
            c += 1

        ax[i].set_title("{0}".format(metric))
        ax[i].yaxis.grid(False)
        #ax[i].grid(False)
        #ax[i].set_ylim(Wdb['midpoint'].max(), 0)
        #ax[i].set_xlim([0, med*2])

        if i == 0:
            ax[i].legend(loc='upper left', title='Min read ANI (%)')

        # Add breaks
        for b in breaks:
            ax[i].axhline(b, ls='-', c='lightgrey', zorder=-1)

        i += 1

    # Annotate the scaffolds
    locations = []
    scaffolds = []
    for scaff, l in midpoints.items():
        locations.append(l)
        scaffolds.append(scaff)
    plt.yticks(locations, scaffolds, fontsize=1)
    plt.ylim(Wdb['midpoint'].max(), 0)

    plt.suptitle(title, y=0.999)
    plt.subplots_adjust(hspace=0.3)
    #plt.gca().xaxis.grid(False)

    fig = plt.gcf()
    ylim = _calc_ylim(Wdb['midpoint'].max())
    fig.set_size_inches(8, ylim)

def _calc_ylim(max_size):
    return min(max(5, max_size/5000), 100)

def linkage_decay_type(Odb, chunkSize=5, min_vals=2, title=''):
    COLS = ['r2', 'r2_normalized', 'd_prime', 'd_prime_normalized']

    # Make chunks
    max_d = Odb['distance'].max()
    table = defaultdict(list)
    numberChunks = max_d // chunkSize + 1
    for lt in ['S=S', 'N-N', 'all']:
        if lt == 'all':
            db = Odb
        else:
            db = Odb[Odb['link_type'] == lt]

        db['distance'] = db['distance'].astype(int)
        for i in range(numberChunks):
            d = db[(db['distance'] >= int(i*chunkSize)) & (db['distance'] < int((i+1)*chunkSize))]

            table['d_start'].append(int(i*chunkSize))
            table['d_end'].append(int((i+1)*chunkSize))
            table['values'].append(len(d))
            table['link_type'].append(lt)
            for col in COLS:
                table[col].append(d[col].mean())
                table[col + '_values'].append(len(d[~d[col].isna()]))
    df = pd.DataFrame(table)
    df['distance'] = [np.mean([x, y]) for x, y in zip(df['d_start'], df['d_end'])]
    df = df.sort_values('distance')

    for thing in ['S-S', 'N-N', 'all']:
        if thing == 'all':
            fdb = df
        else:
            fdb = df[df['link_type'] == thing]
        sns.lineplot(data=fdb, x='distance', y='r2', label=thing, marker='o')

    plt.title(title)
    plt.xlabel('Distance between SNPs (bp)\nAveraged over {0}bp windows; plotting windows with at least {1} values'.format(
                chunkSize, min_vals))
    plt.ylabel("SNP linkage")
    return df

def calc_link_type(row, k2t):
    k1 = "{0}:{1}".format(row['scaffold'], row['position_A'])
    k2 = "{0}:{1}".format(row['scaffold'], row['position_B'])

    if (k1 in k2t) & (k2 in k2t):
        k1 = k2t[k1]
        k2 = k2t[k2]

        return "{0}-{1}".format(k1, k2)

    else:
        return np.nan

def gene_histogram_plot(db, title=''):
    COLS = ['coverage', 'nucl_diversity', 'SNPs_per_bp']

     # Get set up for multiple rows
    i = len(set(db.columns).intersection(set(COLS)))
    if i > 1:
        fig, ax = plt.subplots(i, 1, sharex=True)
    else:
        ax = {}
        ax[0] = plt.gca()

    i = 0
    for metric in COLS:
        if metric not in set(db.columns):
            continue

        df = db[[metric]].sort_values(metric, ascending=False).reset_index(drop=True).reset_index()
        ax[i].plot(df['index'], df[metric], marker='o', ms=1)
        ax[i].set_ylabel("{0}".format(metric))
        i += 1

    plt.xlabel('gene index')
    plt.suptitle(title, y=0.999)
    #plt.xlim(-1, df['index'].max())

def plot_readComparerer_dendrograms(gdb, title, cluster_method='single', thresh=0.001, gthresh = 0.01):
    # Set up the dataframe
    gdb = add_av_RC(gdb)

    # Set up the figure
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'wspace': 0.5})

    #for thing in ['av_cov', 'av_ani']:
    for thing in ['av_ani', 'av_cov']:
        t2a = {'av_ani':ax1, 'av_cov':ax2}
        t2c = {'av_ani':'Average Nucleotide Identity (ANI)', 'av_cov':'Shared Genome Coverage (%)'}
        t2t = {'av_ani':thresh, 'av_cov':gthresh}

        gdb['dist'] = 1 - gdb[thing]
        db = gdb.pivot("name1", "name2", 'dist')

        names = db.columns
        arr =  np.asarray(db)
        arr = ssd.squareform(arr, checks=True)

        # Cluster
        linkage = scipy.cluster.hierarchy.linkage(arr, method=cluster_method)
        fclust = scipy.cluster.hierarchy.fcluster(linkage, thresh,
                            criterion='distance')
        Cdb = drep.d_cluster._gen_cdb_from_fclust(fclust, names)

        if thing == 'av_ani':
            name2cluster = Cdb.set_index('genome')['cluster'].to_dict()
            name2color = drep.d_analyze.gen_color_dictionary(names, name2cluster)

        # Plot
        plt.sca(t2a[thing])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _make_RC_dendrogram(linkage, names, xlabel=t2c[thing], subtitle_string=title, threshold=t2t[thing],
                               name2color=name2color)

    # Fix axis labels 1
    for thing in ['av_cov', 'av_ani']:
        axes = t2a[thing]
        labels = axes.xaxis.get_majorticklocs()
        for i, label in enumerate(labels):
            labels[i] = round((1 - float(label)) * 100, 3)
        axes.set_xticklabels(labels)

def _make_RC_dendrogram(linkage, names, **kwargs):
    '''
    Make the dendrogram used in plot 2
    names can be gotten like:
        db = db.pivot("reference","querry","ani")
        names = list(db.columns)
    Args:
        linkage: result of scipy.cluster.hierarchy.linkage
        names: names of the linkage
    Kwargs:
        name2cluster: dict
        self_thresh: x-axis for soft line
        threshold: x-axis for hard line
        title_sting: title of the plot
        subtitle_string: subtitle of the plot
        winners: list of "winning" genomes (to be marked with star)
        genome2taxonomy: dictionary to add taxonomy information
    Returns:
        Matplotlib primed with a figure
    '''
    # Load possible kwargs
    #name2cluster = kwargs.get('name2cluster',False)
    name2color = kwargs.get('name2color',False)
    self_thresh = kwargs.get('self_thresh',False)
    threshold = kwargs.get('threshold',False)
    threshold = False
    title_string = kwargs.get('title_string','')
    subtitle_string = kwargs.get('subtitle_string','')
    winners = kwargs.get('winners',False)
    genome2taxonomy = kwargs.get('genome2taxonomy',False)
    xlabel = kwargs.get('xlabel','Average Nucleotide Identity (ANI)')

#     # Make special things
#     if name2cluster != False:
#         name2color = drep.d_analyze.gen_color_dictionary(names, name2cluster)
#     else:
#         name2color = False

    # Make the dendrogram
    sns.set_style('whitegrid')
    g = drep.d_analyze.fancy_dendrogram(linkage,names,name2color,threshold=threshold,self_thresh =\
                        self_thresh)

    # Add the title and subtitle
    plt.suptitle(title_string, y=1, fontsize=18)
    plt.title(subtitle_string, fontsize=10)

#     # Adjust the x-axis
    plt.xlabel(xlabel)
    if threshold != False:
        plt.xlim([0,3*threshold])
    plt.tick_params(axis='x', which='major', labelsize=8)
    plt.tick_params(axis='y', which='major', labelsize=12)
#     axes = plt.gca()
#     labels = axes.xaxis.get_majorticklocs()
#     for i, label in enumerate(labels):
#         labels[i] = round((1 - float(label)) * 100, 2)
#     axes.set_xticklabels(labels)
    plt.gca().yaxis.grid(False)

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(20,drep.d_analyze._x_fig_size(len(names),factor=.5))
    plt.subplots_adjust(left=0.5)

    # Mark winning ones
    if type(winners) is not bool:
        ax = plt.gca()
        labels = [item.get_text() for item in ax.get_yticklabels()]
        for i, label in enumerate(labels):
            if label in winners: labels[i] = label + ' *'
        ax.set_yticklabels(labels)

    # Add taxonomy
    if genome2taxonomy != False:
        g2t = genome2taxonomy
        axes = plt.gca()
        labels = [item.get_text() for item in axes.get_yticklabels()]
        for i, label in enumerate(labels):
            labels[i] = "{0}\n{1}".format(label, g2t[label.replace(' *','')])
        axes.set_yticklabels(labels)


def add_av_RC(db):
    '''
    add a column titled 'av_ani' to the passed in dataframe
    dataframe must have rows reference, querey, and ani
    Args:
        db: dataframe
    '''

    combo2value = defaultdict(lambda: np.nan)
    combo2value2 = defaultdict(lambda: np.nan)
    for i, row in db.iterrows():
        combo2value["{0}-vs-{1}".format(row['name1'], row['name2'])] \
            = row['popANI']
        combo2value2["{0}-vs-{1}".format(row['name1'], row['name2'])] \
            = row['coverage_overlap']

    table = defaultdict(list)
    samples = set(db['name1'].tolist()).union(set(db['name2'].tolist()))
    for samp1 in samples:
        for samp2 in samples:
            if samp1 == samp2:
                table['name1'].append(samp1)
                table['name2'].append(samp2)
                table['av_ani'].append(1)
                table['av_cov'].append(1)

            else:
                table['name1'].append(samp1)
                table['name2'].append(samp2)
                table['av_ani'].append(np.nanmean([combo2value["{0}-vs-{1}".format(samp1,samp2)], \
                            combo2value["{0}-vs-{1}".format(samp2,samp1)]]))
                table['av_cov'].append(np.nanmean([combo2value2["{0}-vs-{1}".format(samp1,samp2)], \
                            combo2value2["{0}-vs-{1}".format(samp2,samp1)]]))

    return pd.DataFrame(table)

def _shorten_name(name):
    name = name.replace('.sorted.bam', '')
    if '-vs-' in name:
        name = name.split('-vs-')[1]
    if len(name) > 15:
        name = '\n'.join(name[n:n + 15] for n in range(0, len(name), 15))
    return name

def plot_genome(genome, IS, **kwargs):
    '''
    Decide if this genome should be plotted based on the filters in the kwargs. Return True if so

    GWdb is the result of Mdb = inStrain.genomeUtilities.genomeWideFromIS(IS, 'scaffold_info', mm_level=True)
    '''
    GWdb = kwargs.get('GWdb', False)

    # FILTER BY BREADTH
    mb = float(kwargs.get('minimum_breadth', 0))
    if mb > 0:
        if GWdb is False:
            GWdb = inStrain.genomeUtilities.genomeWideFromIS(IS, 'scaffold_info', mm_level=False)
        if 'mm' in GWdb.columns:
            GWdb = GWdb.sort_values('mm').drop_duplicates(subset='genome', keep='last')
        try:
            breadth = GWdb[GWdb['genome'] == genome]['breadth'].tolist()[0]
        except:
            breadth = 0
        if float(breadth) < mb:
            return False

    # FILTER BY GENOME LIST
    genomes = kwargs.get('genomes', None)
    if genomes is not None:
        if genome not in genomes:
            return False

    return True

'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS IS WHERE THE MAIN METHODS ARE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''

def main(args):
    '''
    Controller for the inStrain plot operation
    '''
    # Parse arguments
    args = validate_input(args)
    if args.ret:
        return

    kwargs = vars(args)
    debug = kwargs.get('debug', False)

    # Load the IS object
    IS = kwargs.pop('IS')

    # Figure out if this is an RC or a IS
    bl = IS.get('bam_loc')
    if bl is None:
        IS_TYPE = 'RC'
        options = ['10']
    else:
        IS_TYPE = 'IS'
        options = ['1','2','3','4','5','6','7','8','9']

    # Figure out what plots to make
    to_plot = _parse_plot_options(options, kwargs.get('plots', None))
    logging.info("making plots {0}".format(', '.join(to_plot)))

    # Get the plot directory and basename
    plot_dir = IS.get_location('figures') + os.path.basename(IS.get('location')) + '_'

    # Cache needed data
    if IS_TYPE == 'IS':
        try:
            kwargs['GWdb'] = IS.get('genome_level_info')
        except:
            logging.error("Cannot cache scaffold info - you don't have all required information. You need to run inStrain genome_wide first")
            if debug:
                traceback.print_exc()

    # Keep cacheing
    if (('2' in to_plot) | ('7' in to_plot)):
        kwargs['covT'] = IS.get('covT')
        kwargs['clonT'] = IS.get('clonT')
        kwargs['raw_linkage_table'] = IS.get('raw_linkage_table')
        kwargs['cumulative_snv_table'] = IS.get('cumulative_snv_table')

    # 1) MM plot
    if '1' in to_plot:
        try:
            mm_plot_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #1: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '2' in to_plot:
        try:
            genome_plot_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #2: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '3' in to_plot:
        try:
            ANI_dist_plot_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #3: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '4' in to_plot:
        try:
            allele_freq_plot_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #4: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '5' in to_plot:
        try:
            linkage_decay_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #5: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '6' in to_plot:
        try:
            read_filtering_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #6: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '7' in to_plot:
        try:
            scaffold_inspection_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #7: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '8' in to_plot:
        try:
            linkage_decay_type_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #8: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '9' in to_plot:
        try:
            gene_histogram_from_IS(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #9: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    if '10' in to_plot:
        try:
            dendrograms_from_RC(IS, plot_dir=plot_dir, **kwargs)
        except BaseException as e:
            logging.error('Failed to make plot #10: ' + str(e))
            if debug:
                traceback.print_exc()
                logging.debug(traceback.format_exc())

    logging.debug("Plotting plots finished")

def validate_input(args):
    '''
    Validate and mess with the arguments a bit
    '''
    # Make sure the IS object is OK
    assert os.path.exists(args.IS)
    args.IS = inStrain.SNVprofile.SNVprofile(args.IS)

    # Set up the logger
    log_loc = args.IS.get_location('log') + 'log.log'
    inStrain.controller.setup_logger(log_loc)

    # Make sure this IS object has an stb
    stb = args.IS.get('scaffold2bin')
    if stb == None:
        logging.error("This IS object does not have an .stb file; cant use it to make plots")
        args.ret = True
    else:
        args.ret = False

    return args

def _parse_plot_options(options, args):
    '''
    Read user input and figure out a list of plots to make

    Args:
        options: list of possible plots to make (default [1-6])
        args: the command line passed in

    Returns:
        list: list of ints in the args
    '''
    to_plot = []

    if args[0] in ['all','a']:
        to_plot += options

    elif args == None:
        logging.error("No plots given!")
        sys.exit()
        return None

    else:
        for arg in args:
            if arg in options:
                to_plot.append(arg)

    return to_plot

def mm_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        Mdb = kwargs.get('GWdb', False)
        assert len(Mdb) > 0

        # Add the number of read-pairs
        #readLen = int(IS.get('mapping_info')['mean_pair_length'].tolist()[0])
        readLen = int(IS.get_read_length())
        Mdb['read_length'] = readLen
        Mdb['mm'] = Mdb['mm'].astype(int)
        Mdb['ANI_level'] = [(readLen - mm)/readLen for mm in Mdb['mm']]
    except:
        logging.error("Skipping plot 1 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 1")
    name = 'CoverageAndBreadth_vs_readMismatch.pdf'
    pp = PdfPages(plot_dir + name)
    #print(Mdb.head())

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        mm_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def genome_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        stb = IS.get('scaffold2bin')
        b2s = defaultdict(list)
        for s, b in stb.items():
            b2s[b].append(s)
        assert len(b2s.keys()) > 0

        # Load the cache
        covTs = kwargs.get('covT')#, IS.get('covT'))
        clonTs = kwargs.get('clonT')#, IS.get('clonT'))
        raw_linkage_table = kwargs.get('raw_linkage_table')#, IS.get('raw_linkage_table'))
        cumulative_snv_table = kwargs.get('cumulative_snv_table')#, IS.get('cumulative_snv_table'))
        scaffold2length = IS.get('scaffold2length')
        rl = IS.get_read_length()
        profiled_scaffolds = set(scaffold2length.keys())

    except:
        logging.error("Skipping plot 2 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 2")
    name = 'genomeWide_microdiveristy_metrics.pdf'
    pp = PdfPages(plot_dir + name)


    for genome, scaffolds in b2s.items():
        if not plot_genome(genome, IS, **kwargs):
            continue
        present_scaffolds = list(set(scaffolds).intersection(set(profiled_scaffolds)))
        Wdb, breaks, midpoints = load_windowed_metrics(present_scaffolds,
                                scaffold2length,
                                rl,
                                report_midpoints=True,
                                covTs=covTs, clonTs=clonTs,
                                raw_linkage_table=raw_linkage_table,
                                cumulative_snv_table=cumulative_snv_table)
        genomeWide_microdiveristy_metrics_plot(Wdb, breaks, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(8, 5)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def ANI_dist_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        Mdb = prepare_read_ani_dist_plot(IS)
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 3 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 3")
    name = 'readANI_distribution.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        db = read_ani_dist_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def allele_freq_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        db = IS.get('cumulative_snv_table')
        if len(db) == 0:
            return
        db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                    .sort_index().drop(columns=['mm'])
        db = db[(db['cryptic'] == False)]
        if 'allele_count' in db.columns:
            db = db[db['allele_count'] >= 2]
        if 'morphia' in db.columns:
            db = db[db['morphia'] >= 2]

        stb = IS.get('scaffold2bin')
        Mdb = inStrain.genomeUtilities._add_stb(db, stb)
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 4 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 4")
    name = 'MajorAllele_frequency_plot.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        db = major_allele_freq_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def linkage_decay_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        db = IS.get('raw_linkage_table')
        db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])

        stb = IS.get('scaffold2bin')
        Mdb = inStrain.genomeUtilities._add_stb(db, stb)
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 5 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 5")
    name = 'LinkageDecay_plot.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        linkage_decay_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

# def read_filtering_from_IS(IS, plot_dir=False, **kwargs):
#     # Load the required data
#     try:
#         # Prepare
#         Mdb = inStrain.genomeUtilities.genomeWideFromIS(IS, 'mapping_info', mm_level=True)
#         print(Mdb)
#         assert len(Mdb) > 0
#     except:
#         logging.error("Skipping plot 6 - you don't have all required information. You need to run inStrain genome_wide first")
#         traceback.print_exc()
#         return
#
#     # Make the plot
#     logging.info("Plotting plot 6")
#     name = 'ReadFiltering_plot.pdf'
#     pp = PdfPages(plot_dir + name)
#
#     for genome, mdb in Mdb.groupby('genome'):
#         if not plot_genome(genome, IS, **kwargs):
#             continue
#         read_filtering_plot(mdb, title=genome)
#         fig = plt.gcf()
#         fig.set_size_inches(6, 4)
#         fig.tight_layout()
#         pp.savefig(fig)#, bbox_inches='tight')
#         #plt.show()
#         plt.close(fig)
#
#     # Save the figure
#     pp.close()
#     #plt.show()
#     plt.close('all')

def read_filtering_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        # Prepare
        Mdb = IS.get('mapping_info')
        Mdb = Mdb[Mdb['scaffold'] == 'all_scaffolds']
        Mdb['genome'] = 'all_scaffolds'
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 6 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 6")
    name = 'ReadFiltering_plot.pdf'
    pp = PdfPages(plot_dir + name)

    read_filtering_plot(Mdb, title='all scaffolds')

    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    fig.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

    # Save the figure
    pp.close()
    plt.close('all')

def scaffold_inspection_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        stb = IS.get('scaffold2bin')
        b2s = defaultdict(list)
        for s, b in stb.items():
            b2s[b].append(s)
        assert len(b2s.keys()) > 0

        # Load the cache
        covTs = kwargs.get('covTs', IS.get('covT'))
        clonTs = kwargs.get('clonTs', IS.get('clonT'))
        raw_linkage_table = kwargs.get('raw_linkage_table', IS.get('raw_linkage_table'))
        cumulative_snv_table = kwargs.get('cumulative_snv_table', IS.get('cumulative_snv_table'))
        scaffold2length = IS.get('scaffold2length')
        rl = IS.get_read_length()
        profiled_scaffolds = set(scaffold2length.keys())

    except:
        logging.error("Skipping plot 7 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 7")
    name = 'ScaffoldInspection_plot.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, scaffolds in b2s.items():
        if not plot_genome(genome, IS, **kwargs):
            continue
        present_scaffolds = list(set(scaffolds).intersection(set(profiled_scaffolds)))
        Wdb, breaks, midpoints = load_windowed_metrics(present_scaffolds,
                                scaffold2length,
                                rl,
                                report_midpoints=True,
                                covTs=covTs, clonTs=clonTs,
                                raw_linkage_table=raw_linkage_table,
                                cumulative_snv_table=cumulative_snv_table)
        scaffold_inspection_plot(Wdb, breaks, midpoints, title=genome)
        fig = plt.gcf()
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def linkage_decay_type_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        # Prepare
        db = IS.get('raw_linkage_table')
        db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])
        stb = IS.get('scaffold2bin')
        Mdb = inStrain.genomeUtilities._add_stb(db, stb)

        SNdb = IS.get('SNP_mutation_types')
        assert SNdb is not None
        if len(SNdb) == 0:
            return
        SNdb['key'] = ["{0}:{1}".format(s, p) for s, p in zip(SNdb['scaffold'], SNdb['position'])]
        k2t = SNdb.set_index('key')['mutation_type'].to_dict()

        Mdb['link_type'] = Mdb.apply(calc_link_type, axis=1, k2t=k2t)
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 8 - you don't have all required information. You need to run inStrain profile_genes first")
        if kwargs.get('debug', False):
            traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 8")
    name = 'LinkageDecay_types_plot.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        db = linkage_decay_type(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def gene_histogram_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        # Prepare
        db = inStrain.GeneProfile.get_gene_info(IS)
        stb = IS.get('scaffold2bin')
        Gdb = inStrain.genomeUtilities._add_stb(db, stb)
        if 'clonality' in Gdb.columns:
            Gdb['nucl_diversity'] = 1 - Gdb['clonality']
        assert len(Gdb) > 0
    except:
        logging.error("Skipping plot 9 - you don't have all required information. You need to run inStrain profile_genes first")
        if kwargs.get('debug', False):
            traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 9")
    name = 'GeneHistogram_plot.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Gdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        db = gene_histogram_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(8, 5)
        fig.tight_layout()
        pp.savefig(fig)#, bbox_inches='tight')
        #plt.show()
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

def dendrograms_from_RC(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        db = IS.get('comparisonsTable')
        stb = IS.get('scaffold2bin')
        b2l = IS.get('bin2length')
        gdb = inStrain.genomeUtilities._add_stb(db, stb)
        Mdb = inStrain.genomeUtilities._genome_wide_readComparer(gdb, stb, b2l, mm_level=False)

        Mdb['name1'] = [_shorten_name(x) for x in Mdb['name1']]
        Mdb['name2'] = [_shorten_name(x) for x in Mdb['name2']]

        Mdb = Mdb.sort_values(['genome', 'name1', 'name2'])
        assert len(Mdb) > 0
    except:
        logging.error("Skipping plot 10 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Plot
    logging.info("Plotting plot 10")
    name = 'inStrainCompare_dendrograms.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        # if not plot_genome(genome, IS, **kwargs):
        #     continue
        plot_readComparerer_dendrograms(mdb, genome, cluster_method='average')
        fig = plt.gcf()
        #fig.tight_layout()
        pp.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # Save the figure
    pp.close()
    #plt.show()
    plt.close('all')

    print('Done!')
