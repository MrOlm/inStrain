#!/usr/bin/env python3
"""
Self-contained script to re-cluster strains from inStrain compare
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

import scipy.cluster
import scipy.spatial.distance

def recluster_instrain(db, cluster_method='single', ANI=99.999):
    '''
    Cluster inStrain compare based on popANI

    If the ANI is set, it will overwrite the "thresh" parameter

    v1.0 - 7/27/21
    '''

    def _gen_cdb_from_fclust(fclust, names):
        '''
        Make Cdb from the result of scipy.cluster.hierarchy.fcluster
        Args:
            fclust: result of scipy.cluster.hierarchy.fcluster
            names: list(db.columns) of the input dataframe
        Returns:
            DataFrame: Cdb
        '''
        Table = {'cluster': [], 'genome': []}
        for i, c in enumerate(fclust):
            Table['cluster'].append(c)
            Table['genome'].append(names[i])

        return pd.DataFrame(Table)

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
                table['name1'].append(samp1)
                table['name2'].append(samp2)

                if samp1 == samp2:
                    table['av_ani'].append(1)
                    table['av_cov'].append(1)

                else:
                    table['av_ani'].append(np.nanmean([combo2value["{0}-vs-{1}".format(samp1, samp2)], \
                                                       combo2value["{0}-vs-{1}".format(samp2, samp1)]]))
                    table['av_cov'].append(np.nanmean([combo2value2["{0}-vs-{1}".format(samp1, samp2)], \
                                                       combo2value2["{0}-vs-{1}".format(samp2, samp1)]]))
        return pd.DataFrame(table)

    thresh = 1 - (ANI / 100)

    # Set up the dataframe
    gdb = add_av_RC(db)
    gdb['dist'] = 1 - gdb['av_ani']

    missing_comps = len(gdb[gdb['dist'].isna()])
    if missing_comps > 0:
        print(
            f"WARNING! {missing_comps} / {len(gdb)} ({(missing_comps / len(gdb)) * 100:.2f}%) of comparisons are missing")

    gdb['dist'] = gdb['dist'].fillna(1)
    db = gdb.pivot(index="name1", columns="name2", values='dist')

    # Set up some more
    names = db.columns
    try:
        arr = np.asarray(db)
        arr = scipy.spatial.distance.squareform(arr, checks=True)
    except:
        print(db)
        raise Exception

    # Cluster
    # print("thresh is {0}".format(thresh))
    linkage = scipy.cluster.hierarchy.linkage(arr, method=cluster_method)
    fclust = scipy.cluster.hierarchy.fcluster(linkage, thresh,
                                              criterion='distance')
    Cdb = _gen_cdb_from_fclust(fclust, names)
    Cdb = Cdb.rename(columns={'genome': 'sample'})

    return Cdb

def main(args):
    """
    Quick wrapper for this
    """
    odb = pd.read_csv(args.db, sep='\t')
    counter = 1
    dbs = []
    for genome, odb in odb.groupby('genome'):
        cdb = recluster_instrain(odb, cluster_method=args.clusterAlg, ANI=args.ani)
        cdb['genome'] = genome
        cdb['cluster'] = [f"{counter}_{c}" for c in cdb['cluster']]
        counter += 1
        dbs.append(cdb)
    pd.concat(dbs).reset_index(drop=True).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("-d", "--db", help="Location of genomeWide_compare.tsv file (from inStrain compare output)", required=True)
    parser.add_argument("-o", "--output", help="output filename", default='strain_clusters_reClustered.tsv')
    parser.add_argument("-a", "--ani", help="popANI threshold (as a percentage, not fraction; default is 99.999)", default=99.999, type=float)
    parser.add_argument("--clusterAlg", help="Algorithm used to cluster genomes (passed\
                            to scipy.cluster.hierarchy.linkage", default='single',
                           choices={'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'})

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))


    args = parser.parse_args()
    main(args)