#!/usr/bin/env python

import os
import copy
import time
import logging
import pandas as pd
from tqdm import tqdm
import multiprocessing
import traceback
from collections import defaultdict
import numpy as np
import scipy.spatial.distance

import inStrain
import inStrain.logUtils
import inStrain.readComparer

def run_compare_multiprocessing(SC_queue, result_queue, null_model, num_to_run, **kwargs):
    """
    Run the multiprocessing associated with the "compare" operation
    """
    # Get kwargs
    p = int(kwargs.get('processes', 6))

    # Set up progress bar
    results = []
    pbar = tqdm(desc='Comparing scaffolds: ', total=num_to_run)

    # Launch the multiprocessing
    processes = []
    inStrain.logUtils.log_checkpoint("Compare", "multiprocessing", "start")
    if p > 1:
        ctx = multiprocessing.get_context('spawn')
        for i in range(0, p):
            processes.append(ctx.Process(target=run_SC_objects, args=(SC_queue, result_queue, null_model, kwargs)))
        for proc in processes:
            proc.start()
    else:
        run_SC_objects(SC_queue, result_queue, null_model, kwargs, single_thread=True)

    # Get the results
    recieved_groups = 0
    while recieved_groups < num_to_run:
        result, log = result_queue.get()
        logging.debug(log)
        if result is not None:
            results.append(result)
        recieved_groups += 1
        pbar.update(1)

    # Close multi-processing
    for proc in processes:
        proc.terminate()

    # Close progress bar
    pbar.close()
    inStrain.logUtils.log_checkpoint("Compare", "multiprocessing", "end")

    return results

def run_SC_objects(cmd_queue, result_queue, null_model, kwargs, single_thread=False):
    """
    The actually paralellized unit for running compare
    """
    # Apply multiprocessing patch
    inStrain.controller.patch_mp_connection_bpo_17560()

    # Continually push as you run
    while True:
        # Get an SC object
        if single_thread:
            try:
                cmd = cmd_queue.get(timeout=5)
            except:
                return
        else:
            cmd = cmd_queue.get(True)

        # Process that SC object
        SC_result, log = SC_object_wrapper(cmd, null_model, kwargs)
        result_queue.put((SC_result, log))

        # Clean up memory
        for r in SC_result:
            del r

def SC_object_wrapper(SC, null_model, kwargs):
    '''
    Take a Scaffold Compare object and profile the scaffold
    '''
    debug = kwargs.get('debug', False)
    try:
        results, log = inStrain.readComparer.compare_scaffold(SC.scaffold, SC.names, SC.SNPtables, SC.covTs,
                                                                      SC.length, null_model, **kwargs)

    except Exception as e:
        if debug:
            print(e)
            traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(SC.scaffold)))
        t = time.strftime('%m-%d %H:%M')
        log_message = "\n{1} DEBUG FAILURE CompareScaffold {0} {2}\n" \
            .format(SC.scaffold, t, str(SC.names))
        results = None
        log = log_message

    return (results, log)


def subset_SNP_table(db, scaffold):
    if len(db) > 0:
        db = db[db['scaffold'] == scaffold]
        if len(db) == 0:
            db = pd.DataFrame()
        else:
            db = db.sort_values('mm')
            #db = db[['position', 'mm', 'con_base', 'ref_base', 'var_base', 'position_coverage', 'A', 'C', 'T', 'G', 'allele_count']]
    else:
        db = pd.DataFrame()

    return db

def find_relevant_scaffolds(input, bts, kwargs):
    """
    Return a list of scaffolds in the input based on the parameters of kwargs
    """
    GIdb = inStrain.SNVprofile.SNVprofile(input).get('genome_level_info')
    if GIdb is None:
        logging.error(f"Profile {input} does not have genome-level information; needed to run compare in database mode")
        raise Exception
    if 'mm' in GIdb:
        GIdb = GIdb.sort_values('mm', ascending=True).drop_duplicates(subset=['genome'], keep='last')

    min_breadth = kwargs.get('breadth', 0.5)
    genomes = GIdb[GIdb['breadth_minCov'] >= min_breadth]['genome'].tolist()

    scaffolds = []
    for genome in genomes:
        if genome in bts:
            scaffolds += bts[genome]
        else:
            logging.error(f'{genome} is in input {input} but not the provided stb file!')
            raise Exception(f'{genome} is in input {input} but not the provided stb file!')

    message = f'{input} has {len(genomes)} genomes detected and {len(scaffolds)} scaffolds'
    print(message)
    logging.info(message)

    return set(scaffolds)

def add_av_RC(db, v1='popANI', v2='coverage_overlap', n1='av_ani', n2='av_cov'):
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
            = row[v1]
        combo2value2["{0}-vs-{1}".format(row['name1'], row['name2'])] \
            = row[v2]

    table = defaultdict(list)
    samples = set(db['name1'].tolist()).union(set(db['name2'].tolist()))
    for samp1 in samples:
        for samp2 in samples:
            if samp1 == samp2:
                table['name1'].append(samp1)
                table['name2'].append(samp2)
                table[n1].append(1)
                table[n2].append(1)

            else:
                table['name1'].append(samp1)
                table['name2'].append(samp2)
                table[n1].append(np.nanmean([combo2value["{0}-vs-{1}".format(samp1,samp2)],
                            combo2value["{0}-vs-{1}".format(samp2,samp1)]]))
                table[n2].append(np.nanmean([combo2value2["{0}-vs-{1}".format(samp1,samp2)],
                            combo2value2["{0}-vs-{1}".format(samp2,samp1)]]))

    return pd.DataFrame(table)

def cluster_genome_strains(Mdb, kwargs):
    """
    Perform dRep-style clustering on the genomes
    """
    cluster_method = kwargs.get('clusterAlg')
    thresh = 1 - kwargs.get('ani_threshold')
    cov_thresh = kwargs.get('coverage_treshold')

    cdbs = []

    cluster_num = 1
    for genome, gdb in Mdb.groupby('genome'):

        if not evalute_genome_dist_matrix(gdb, genome):
            continue

        # Average popANI values
        gdb = add_av_RC(gdb, v2='percent_compared', n2='av_cov')
        gdb['dist'] = 1 - gdb['av_ani']

        # Remove values with low coverage
        gdb['dist'] = [1 if c < cov_thresh else d for d, c in zip(gdb['dist'], gdb['av_cov'])]

        # Make squareform
        db = gdb.pivot("name1", "name2", 'dist')
        names = db.columns
        arr = np.asarray(db)
        arr = scipy.spatial.distance.squareform(arr, checks=True)

        # Cluster
        linkage = scipy.cluster.hierarchy.linkage(arr, method=cluster_method)
        fclust = scipy.cluster.hierarchy.fcluster(linkage, thresh,
                                                  criterion='distance')

        # Get Cdb
        cdb = _gen_cdb_from_fclust(fclust, names)
        cdb = cdb.rename(columns={'genome':'sample'})
        cdb['genome'] = genome

        # Parse cdb
        cdb['cluster'] = [f'{cluster_num}_{x}' for x in cdb['cluster']]
        cluster_num += 1

        # store
        cdbs.append(cdb)

    return pd.concat(cdbs).reset_index(drop=True)

def evalute_genome_dist_matrix(mdb, genome):
    """
    Check whether the distance matrix is OK. Return True if OK
    """
    edb = mdb[mdb['compared_bases_count'] == 0]
    if len(edb) > 0:
        logging.error(f"Cannot cluster genome {genome}; {len(edb)} of {len(mdb)} comaprisons involve no genomic overlap at all: see log for more")
        for i, row in edb.iterrows():
            logging.debug(f"clustering failure {genome} \t {row['name1']} vs {row['name2']}")
        return False
    else:
        return True

def _gen_cdb_from_fclust(fclust,names):
    '''
    Make Cdb from the result of scipy.cluster.hierarchy.fcluster
    Args:
        fclust: result of scipy.cluster.hierarchy.fcluster
        names: list(db.columns) of the input dataframe
    Returns:
        DataFrame: Cdb
    '''
    Table={'cluster':[],'genome':[]}
    for i, c in enumerate(fclust):
        Table['cluster'].append(c)
        Table['genome'].append(names[i])

    return pd.DataFrame(Table)

