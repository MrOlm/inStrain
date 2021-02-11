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

def gene_histogram_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        # Prepare
        db = inStrain.GeneProfile.get_gene_info(IS)
        stb = IS.get('scaffold2bin')
        Gdb = inStrain.genomeUtilities._add_stb(db, stb, verbose=False)
        if 'clonality' in Gdb.columns:
            Gdb.loc[:,'nucl_diversity'] = 1 - Gdb['clonality']
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

def gene_histogram_plot(db, title=''):
    COLS = ['coverage', 'nucl_diversity', 'SNPs_per_bp']
    sns.set_style('white')

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
        ax[i].axvline(0, c='black')
        ax[i].axhline(0, c='black')
        ax[i].plot(df['index'], df[metric], marker='o', ms=1)
        ax[i].set_ylabel("{0}".format(metric))
        i += 1

    plt.xlabel('gene index')
    plt.suptitle(title, y=0.999)
    plt.axvline(0)
    #plt.xlim(-1, df['index'].max())