import logging
import traceback
import pandas as pd
import seaborn as sns

from collections import defaultdict
import numpy as np

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


def mm_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        Mdb = kwargs.get('GWdb', False)
        assert len(Mdb) > 0

        if 'mm' not in Mdb:
            raise Exception('Plot 1 cannot be created when run with --database_mode or --skip_mm_profiling')

        # Add the number of read-pairs
        readLen = int(IS.get_read_length())
        Mdb['read_length'] = readLen
        Mdb['mm'] = Mdb['mm'].astype(int)
        Mdb.loc[:,'ANI_level'] = [(readLen - mm) / readLen for mm in Mdb['mm']]
    except:
        logging.error(
            "Skipping plot 1 - you don't have all required information. You need to run inStrain genome_wide first")
        traceback.print_exc()
        return

    # Make the plot
    logging.info("Plotting plot 1")
    name = 'CoverageAndBreadth_vs_readMismatch.pdf'
    pp = PdfPages(plot_dir + name)

    for genome, mdb in Mdb.groupby('genome'):
        if not plot_genome(genome, IS, **kwargs):
            continue
        mm_plot(mdb, title=genome)
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        pp.savefig(fig)
        plt.close(fig)

    # Save the figure
    pp.close()
    plt.close('all')

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

def ANI_dist_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        Mdb = prepare_read_ani_dist_plot(IS)
        if len(Mdb['ANI_level'].unique()) == 1:
            raise Exception('Plot 3 cannot be created when run with --database_mode or --skip_mm_profiling')
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

def mm_plot(db, left_val='breadth', right_val='coverage', title='',
            minANI=0.9):
    """
    The input db for this is "mm_genome_info" from "makeGenomeWide" in genomeUtilities.py
    """
    db = db.sort_values('ANI_level')
    sns.set_style('white')

    # breadth
    fig, ax1 = plt.subplots()
    ax1.plot(db['ANI_level'], db[left_val], ls='-', color='blue')
    if left_val == 'breadth':
        ax1.plot(db['ANI_level'], estimate_breadth(db['coverage']), ls='--', color='lightblue')
    ax1.set_ylabel(left_val, color='blue')
    ax1.set_xlabel('Minimum read ANI level')
    ax1.set_ylim(0, 1)

    # coverage
    ax2 = ax1.twinx()
    ax2.plot(db['ANI_level'], db[right_val], ls='-', color='red')
    ax2.set_ylabel(right_val, color='red')
    ax2.set_ylim(0, )

    # asthetics
    plt.xlim(1, max(minANI, db['ANI_level'].min()))

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
    gdb = inStrain.genomeUtilities._add_stb(db, stb, verbose=False)

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
    db.loc[:,'read_pairs'] = [int((x*y) / (readLen * 2)) for x, y in zip(db['coverage'], db['length'])]
    db.loc[:,'ANI_level'] = [(readLen - mm)/readLen for mm in db['mm']]

    return db

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
    db.loc[:,'variable'] = [c2c[x] if x in c2c else x for x in db['variable']]
    db.loc[:,'value'] = [int(x/2) if y == 'Total mapping reads (divided by 2)' else x for x, y in zip(
                db['value'], db['variable'])]

    # Set up colors
    v2c = {v:'grey' for v in db['variable'].unique()}
    v2c['Total filtered pairs'] = 'green'

    db = db.sort_values(['value', 'variable'], ascending=False)
    ax = sns.barplot(x=db['value'], y=db['variable'], palette=v2c)
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

    # Remove top and right axes
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.title(title)
