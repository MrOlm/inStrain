import logging
import traceback
import pandas as pd
import seaborn as sns

from collections import defaultdict
import numpy as np

import matplotlib
import warnings
import scipy
import scipy.spatial.distance as ssd

matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42

import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import inStrain.compare_utils
import inStrain.plotting.utilities
from inStrain.plotting.utilities import plot_genome
from inStrain.plotting.utilities import estimate_breadth

def dendrograms_from_RC(IS, plot_dir=False, **kwargs):
    # Load the required data
    try:
        db = IS.get('comparisonsTable')
        stb = IS.get('scaffold2bin')
        b2l = IS.get('bin2length')
        gdb = inStrain.genomeUtilities._add_stb(db, stb)
        Mdb = inStrain.genomeUtilities._genome_wide_readComparer(gdb, stb, b2l, mm_level=False)

        Mdb['name1'] = [inStrain.plotting.utilities._shorten_name(x) for x in Mdb['name1']]
        Mdb['name2'] = [inStrain.plotting.utilities._shorten_name(x) for x in Mdb['name2']]

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

        # Evaluate
        if not inStrain.compare_utils.evalute_genome_dist_matrix(mdb, genome):
            continue

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
        Cdb = inStrain.compare_utils._gen_cdb_from_fclust(fclust, names)

        if thing == 'av_ani':
            name2cluster = Cdb.set_index('genome')['cluster'].to_dict()
            name2color = inStrain.plotting.utilities.gen_color_dictionary(names, name2cluster)

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
    g = fancy_dendrogram(linkage,names,name2color,threshold=threshold,self_thresh =\
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
    fig.set_size_inches(20,inStrain.plotting.utilities._x_fig_size(len(names),factor=.5))
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

def fancy_dendrogram(linkage,names,name2color=False,threshold=False,self_thresh=False):
    '''
    Make a fancy dendrogram
    '''
    # Make the dendrogram
    if threshold == False:
        scipy.cluster.hierarchy.dendrogram(linkage,labels=names,orientation='right')
    else:
        scipy.cluster.hierarchy.dendrogram(linkage,labels=names, color_threshold=threshold,\
                                            orientation='right')

    # Color the names
    if name2color != False:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            lbl.set_color(name2color[lbl.get_text()])

    # Add the threshold
    if threshold:
        plt.axvline(x=threshold, c='k', linestyle='dotted')
    if self_thresh != False:
        plt.axvline(x=self_thresh, c='red', linestyle='dotted', lw=1)

    g = plt.gcf()
    return g