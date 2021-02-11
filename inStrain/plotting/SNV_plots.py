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

def allele_freq_plot_from_IS(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        if not hasattr(sns, 'histplot') and callable(getattr(sns, 'histplot')):
            raise Exception("Cannot make plot 4 because your seaborn is out of date- need v0.11+")

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
        Mdb = inStrain.genomeUtilities._add_stb(db, stb, verbose=False)
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

def major_allele_freq_plot(db, title=''):
    db.loc[:, 'major_allele_freq'] = [max(x, y) for x, y in zip(db['var_freq'], db['ref_freq'])]
    sns.histplot(db['major_allele_freq'], bins=np.arange(0.5, 1, 0.01), kde=False, binwidth=0.005)

    plt.xlim(0.5, 1)
    plt.title(title)
    plt.xlabel('Major allele frequency')
    plt.ylabel('Number of SNPs')