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

def linkage_decay_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        db = IS.get('raw_linkage_table')
        db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])

        stb = IS.get('scaffold2bin')
        Mdb = inStrain.genomeUtilities._add_stb(db, stb, verbose=False)
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

def linkage_decay_type_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        # Prepare
        db = IS.get('raw_linkage_table')
        db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_index().drop(columns=['mm'])
        stb = IS.get('scaffold2bin')
        Mdb = inStrain.genomeUtilities._add_stb(db, stb, verbose=False)

        SNdb = IS.get('SNP_mutation_types')
        assert SNdb is not None
        if len(SNdb) == 0:
            return
        SNdb.loc[:,'key'] = ["{0}:{1}".format(s, p) for s, p in zip(SNdb['scaffold'], SNdb['position'])]
        k2t = SNdb.set_index('key')['mutation_type'].to_dict()

        Mdb.loc[:,'link_type'] = Mdb.apply(calc_link_type, axis=1, k2t=k2t)
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

def linkage_decay_plot(db, chunkSize=5, min_vals=5, title=''):
    COLS = ['r2', 'r2_normalized', 'd_prime', 'd_prime_normalized']
    sns.set_style('white')

    # Make chunks
    max_d = db['distance'].max()
    table = defaultdict(list)
    numberChunks = max_d // chunkSize + 1
    #db['distance'] = db['distance'].astype(int)
    db = db.astype({'distance': int})
    for i in range(numberChunks):
        d = db[(db['distance'] >= int(i*chunkSize)) & (db['distance'] < int((i+1)*chunkSize))]

        table['d_start'].append(int(i*chunkSize))
        table['d_end'].append(int((i+1)*chunkSize))
        table['values'].append(len(d))
        for col in COLS:
            table[col].append(d[col].mean())
            table[col + '_values'].append(len(d[~d[col].isna()]))
    df = pd.DataFrame(table)
    df.loc[:,'distance'] = [np.mean([x, y]) for x, y in zip(df['d_start'], df['d_end'])]
    df = df.sort_values('distance')

    for col in COLS:
        df[col] = [c if v >= min_vals else np.nan for c, v in zip(df[col], df[col + '_values'])]

    for col in COLS:
        sns.lineplot(x=df['distance'], y=df[col], label=col, marker='o')

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

def linkage_decay_type(Odb, chunkSize=5, min_vals=2, title=''):
    sns.set_style('white')

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

        #db['distance'] = db['distance'].astype(int)
        db = db.astype({'distance':int})
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
    df.loc[:,'distance'] = [np.mean([x, y]) for x, y in zip(df['d_start'], df['d_end'])]
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