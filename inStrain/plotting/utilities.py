import sys
import logging
import numpy as np
from matplotlib import pyplot as plt

import inStrain
import inStrain.genomeUtilities

def plot_genome(genome, IS, **kwargs):
    '''
    Decide if this genome should be plotted based on the filters in the kwargs. Return True if so

    GWdb is the result of Mdb = inStrain.genomeUtilities.genomeWideFromIS(IS, 'scaffold_info', mm_level=True)
    '''
    GWdb = kwargs.get('GWdb', False)

    # FILTER BY GENOME LIST
    genomes = kwargs.get('genomes', None)
    if genomes is not None:
        if genome not in genomes:
            return False
        else:
            return True

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

    return True

def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage

    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1

def _calc_ylim(max_size):
    return min(max(5, max_size/5000), 100)

def _shorten_name(name):
    name = name.replace('.sorted.bam', '')
    if '-vs-' in name:
        name = name.split('-vs-')[1]
    if len(name) > 15:
        name = '\n'.join(name[n:n + 15] for n in range(0, len(name), 15))
    return name

def gen_color_dictionary(names, name2cluster):
    '''
    Make the dictionary name2color
    Args:
        names: key in the returned dictionary
        name2cluster: a dictionary of name to it's cluster
    Returns:
        dict: name -> color
    '''
    #cm = _rand_cmap(len(set(name2cluster.values()))+1,type='bright')
    vals = np.linspace(0,1,len(set(name2cluster.values()))+1)
    np.random.shuffle(vals)
    cm = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

    # 1. generate cluster to color
    cluster2color = {}
    clusters = set(name2cluster.values())
    NUM_COLORS = len(clusters)
    for cluster in clusters:
        try:
            cluster2color[cluster] = cm(1.*int(cluster)/NUM_COLORS)
        except:
            cluster2color[cluster] = cm(1.*int(str(cluster).split('_')[1])/NUM_COLORS)

    #2. name to color
    name2color = {}
    for name in names:
        name2color[name] = cluster2color[name2cluster[name]]

    return name2color

def _x_fig_size(points, factor= .07, min= 8):
    '''
    Calculate how big the x of the figure should be
    '''
    size = points * factor
    return max([size,min])