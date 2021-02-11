import logging
import traceback
import pandas as pd
import numpy as np
import seaborn as sns
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42

import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import inStrain.plotting.utilities
from inStrain.plotting.utilities import plot_genome
from inStrain.plotting.utilities import estimate_breadth

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
        if len(Wdb) == 0:
            logging.debug(f"{genome} could not have windowed metrics loaded")
            continue
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
        if len(Wdb) == 0:
            logging.debug(f"{genome} could not have windowed metrics loaded")
            continue
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
    ylim = inStrain.plotting.utilities._calc_ylim(Wdb['midpoint'].max())
    fig.set_size_inches(8, ylim)