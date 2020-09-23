# Import packages
import os
import gc
import sys
import time
import pysam
import scipy
import psutil
import random
import logging
import resource
import traceback
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing
import concurrent.futures
from subprocess import call
from concurrent import futures
from collections import defaultdict
from numba import jit

import inStrain.profile.linkage
import inStrain.SNVprofile
import inStrain.readComparer
import inStrain.logUtils
import inStrain.GeneProfile
import inStrain.controller

from inStrain import __version__

'''
The following are globals used during processing
'''
P2C = {'A':0, 'C':1, 'T':2, 'G':3} # base -> position
C2P = {0:'A', 1:'C', 2:'T', 3:'G'} # position -> base

def split_profile_worker(split_cmd_queue, Sprofile_dict, log_list,
                        null_model, bam, single_thread=False):
    '''
    A worker to profile splits that can be multi-processed

    Kind of a wrapper for the function "split_profile_wrapper_groups"

    Args:
        split_cmd_queue: A queue of split commands to grab
        Sprofile_dict: A dictionary to store processed splits
        log_list: A queue to put logs

        null_model: Used for SNP profiling (needs to be pre-generated)
        bam: Location of .bam file

    Keyword Args:
        singe_thread: Run in single_thread mode
    '''
    # Initilize the .bam file
    bam_init = samfile = pysam.AlignmentFile(bam)

    # Apply patch
    inStrain.controller.patch_mp_connection_bpo_17560()

    j = 0
    pid = os.getpid()
    while True:
        j += 1
        group_name = "{0}_{1}".format(pid, j)
        group_log = inStrain.logUtils.get_group_log('SplitProfile', group_name, 'start')

        # Get command
        if not single_thread:
            cmds = split_cmd_queue.get(True)
        else:
            try:
                cmds = split_cmd_queue.get(timeout=5)
            except:
                return

        # Process split
        Splits = split_profile_wrapper_groups(cmds, null_model, bam_init)
        LOG = ''
        for Split in Splits:
            # This is an error
            if type(Split) == type('string'):
                LOG += Split
                continue
            Sprofile_dict[Split.scaffold + '.' + str(Split.split_number)] = Split
            LOG += Split.log + '\n'

        group_log += inStrain.logUtils.get_group_log('SplitProfile', group_name, 'end')
        LOG += group_log + '\n'
        log_list.put(LOG)

def split_profile_wrapper_groups(cmds, null_model, bam_init):
    '''
    A wrapper to profile a list of splits and handle exceptions

    This is really based around the method "profile_split"
    '''
    results = []
    for cmd in cmds:
        try:
            results.append(profile_split(bam_init, cmd.scaffold, cmd.start,
                            cmd.end, cmd.split_number, cmd.sequence, cmd.R2M,
                            null_model, bam_name=cmd.samfile, **cmd.arguments))
        except Exception as e:
            print(e)
            traceback.print_exc()

            t = time.strftime('%m-%d %H:%M')
            log_message = "\n{1} DEBUG FAILURE SplitException {0} {2}\n"\
                            .format(cmd.scaffold, t, cmd.split_number)
            results.append(log_message)
    return results


def profile_split(samfile, scaffold, start, end, split_number,
                    seq, R2M, null_model, **kwargs):
    '''
    The heart of inStrain split profiling

    Args:
        samfile: A pre-initialized pysam.AlignmentFile
        scaffold: The name of the scaffold to profile
        start: Where to start profiling the scaffold (inclusive 0-based)
        end: Where to end profiling the scaffold (inclusive 0-based)
        split_number: Int representation of the split number
        seq: The sequence of the scaffold
        R2M: A dictionary of read -> number of mismatches
        null_model: An initialized null_model

    Returns:
        Sprofile: A SplitObject()
    '''
    # Log
    unit = "{0}.{1}".format(scaffold, split_number)
    log_message = inStrain.logUtils.get_worker_log('SplitProfile', unit, 'start')

    # For testing purposes
    if ((scaffold == 'FailureScaffoldHeaderTesting') & (split_number == 1) & (kwargs.get('debug', False))):
        assert False

    # Get kwargs
    min_cov = int(kwargs.get('min_cov', 5))
    min_covR = int(kwargs.get('rarefied_coverage', 5))
    min_freq = float(kwargs.get('min_freq', .05))
    min_snp = int(kwargs.get('min_snp', 10))
    store_everything = kwargs.get('store_everything', False)

    # Set up the .bam iterator
    try:
        iter = samfile.pileup(scaffold, truncate=True, max_depth=100000,
                                stepper='nofilter', compute_baq=True,
                                ignore_orphans=True, ignore_overlaps=True,
                                min_base_quality=30, start=start, stop=end+1)
    except ValueError:
        logging.error("scaffold {0} is not in the .bam file {1}!".format(scaffold, samfile))
        return None, log_message

    # Initialize
    mLen = len(seq) # Length of sequence
    covT = {} # Dictionary of mm -> positional coverage along the genome
    clonT = {} # Diciontary of mm -> clonality
    clonTR = {} # Diciontary of mm -> rarefied clonality
    p2c = {} # dictionary of position -> cryptic SNPs
    read_to_snvs = defaultdict(_dlist) # Dictionary of mm -> read variant -> count
    snv2mm2counts = {} # dictionary of position to mm to counts for SNP positions
    Stable = defaultdict(list) # Holds SNP information
    if store_everything:
        pileup_counts = np.zeros(shape=(mLen, 4), dtype=int) # Holds all pileup counts - alexcc 5/9/2019
    else:
        pileup_counts = None

    # Do per-site processing on all these objects
    process_bam_sites(scaffold, seq, iter, covT, clonT, clonTR, p2c,
                      read_to_snvs, snv2mm2counts, Stable, pileup_counts,
                      mLen, null_model, R2M, start=start, **kwargs)

    # Shrink these dictionaries into index pandas Series
    covT = shrink_basewise(covT, 'coverage', start=start, len=mLen)
    clonT = shrink_basewise(clonT, 'clonality', start=start, len=mLen)
    clonTR = shrink_basewise(clonTR, 'clonality', start=start, len=mLen)

    # Make the SNP table
    SNPTable = inStrain.profile.snv_utilities.generate_snp_table(Stable, scaffold, p2c)
    if len(SNPTable) > 0:
        SNPTable['position'] = SNPTable['position'] + start

    # Make linkage network and calculate linkage
    mm_to_position_graph = inStrain.profile.linkage.calc_mm_SNV_linkage_network(read_to_snvs, scaff=scaffold)
    LDdb = inStrain.profile.linkage.calculate_ld(mm_to_position_graph, min_snp, snv2mm2counts=snv2mm2counts, scaffold=scaffold)
    if len(LDdb) > 0:
        for p in ['position_A', 'position_B']:
            LDdb[p] = LDdb[p] + start

    # Make a Scaffold profile to return
    Sprofile = SplitObject()
    Sprofile.scaffold = scaffold
    Sprofile.split_number = split_number
    Sprofile.bam = kwargs.get('bam_name')
    Sprofile.length = mLen
    Sprofile.raw_snp_table = SNPTable
    Sprofile.raw_linkage_table = LDdb
    Sprofile.covT = covT
    Sprofile.clonT = clonT
    Sprofile.clonTR = clonTR
    Sprofile.min_freq = min_freq

    # store extra things if required
    if store_everything:
        for att in ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts']:
            setattr(Sprofile, att, eval(att))

    # Return
    log_message += inStrain.logUtils.get_worker_log('SplitProfile', unit, 'end')
    Sprofile.log = log_message

    return Sprofile

def process_bam_sites(scaffold, seq, iter, covT, clonT, clonTR, p2c,
                        read_to_snvs, snv2mm2counts, Stable, pileup_counts,
                        mLen, null_model, R2M, start=0, **kwargs):
    '''
    Iterate through the iterator and fill in the information

    Args:
        scaffold: Name of the scaffold
        seq: Sequence of the scaffold
        iter: The "samfile.pileup" object
        covT, clonT, clonTR, p2c, read_to_snvs, snv2mm2counts, Stable, pileup_counts: Pre-initialized arrays
        mLen: Lnegth of sequence
        R2M: A dictionary of read -> number of mismatches
        null_model: An initialized null_model
    '''
    # Get kwargs
    min_cov = int(kwargs.get('min_cov', 5))
    min_covR = int(kwargs.get('rarefied_coverage', 5))
    min_freq = float(kwargs.get('min_freq', .05))
    min_snp = int(kwargs.get('min_snp', 10))
    store_everything = kwargs.get('store_everything', False)

    for pileupcolumn in iter:
        # NOTE! If you have paired reads in a pileup column, only one will come back starting in pysam v0.15 (which is usually the intended functionality) - M.O. 3.19.19

        TruePosition =  pileupcolumn.pos
        RelPosition = pileupcolumn.pos - start

        # Get basic base counts
        MMcounts = get_base_counts_mm(pileupcolumn, R2M)
        update_covT(covT, MMcounts, RelPosition, mLen)

        # Call SNPs
        snp, bases, total_counts = inStrain.profile.snv_utilities.update_snp_table(
                                    Stable, clonT, clonTR, MMcounts, p2c,
                                    RelPosition, scaffold, mLen, seq[RelPosition],
                                    null_model, min_cov=min_cov, min_covR=min_covR,
                                    min_freq=min_freq)

        # add the counts for this position to the numpy array alexcc
        if store_everything:
            pileup_counts[RelPosition] = total_counts

        # get linked reads
        if snp:
            inStrain.profile.linkage.update_linked_reads(read_to_snvs,
                                    pileupcolumn, MMcounts, RelPosition,
                                    bases, R2M, scaffold=scaffold)
            snv2mm2counts[RelPosition] = MMcounts

def get_base_counts_mm(pileupcolumn, R2M):
    '''
    From a pileupcolumn object, return a dictionary of readMismatches ->
        list with the counts of [A, C, T, G]
    '''
    table = defaultdict(_4zeros)
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            try:
                if type(R2M) == type({}):
                    table[R2M[pileupread.alignment.query_name]]\
                    [P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
                elif pileupread.alignment.query_name in R2M:
                    table[0]\
                    [P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
            except KeyError: # This would be like an N or something not A/C/T/G
                pass
    return table

def update_covT(covT, MMcounts, position, mLen):
    '''
    Update covT at this position
    '''
    for mm, count in MMcounts.items():
        if mm not in covT:
            covT[mm] = np.zeros(mLen, dtype=int)
        covT[mm][position] = sum(count)

def mm_counts_to_counts(MMcounts, maxMM=100):
    '''
    Take mm counts and return just counts
    '''
    counts = None
    for mm, count in [(mm, count) for mm, count in MMcounts.items() if mm <= maxMM]:
        if counts is None:
            counts = count
        else:
            counts = np.add(counts, count)

    if counts is None:
        return np.zeros(4, dtype=int)

    else:
        return counts

# def mm_counts_to_counts(MMcounts, maxMM=100):
#     '''
#     Take mm counts and return just counts
#     '''
#     counts = np.zeros(4, dtype=int)
#
#     mms = np.array(list(MMcounts.keys()), dtype='int32')
#     covs = np.array(list(MMcounts.values()))
#
#     return _mm_counts_to_counts_fast(mms, covs, counts, maxMM)
#
# @jit(nopython=True)
# def _mm_counts_to_counts_fast(mms, covs, counts, maxMM):
#     '''
#     Fast implementation
#     '''
#     i = 0
#     for mm in mms:
#         if mm <= maxMM:
#             counts = np.add(counts, covs[i])
#         i += 1
#     return counts

def shrink_basewise(mm2array, name, start=0, len=0):
    NAME2TYPE = {'coverage':'int32', 'clonality':'float32', 'snpCounted':'bool'}

    mm2bases = {}
    for mm in list(mm2array.keys()):
        if start != 0:
            db = pd.Series(mm2array[mm], index=pd.RangeIndex(start=start, stop=start+len), dtype=NAME2TYPE[name]).dropna()
        else:
            db = pd.Series(mm2array[mm], dtype=NAME2TYPE[name]).dropna()
        db = db[db > 0]
        mm2bases[mm] = db
        del mm2array[mm]

    return mm2bases

def merge_profile_worker(sprofile_cmd_queue, Sprofile_dict, Sprofiles,
                        null_model, single_thread=False):
    '''
    Worker to merge_splits
    '''
    # Apply patch
    inStrain.controller.patch_mp_connection_bpo_17560()

    j = 0
    pid = os.getpid()
    while True:
        j += 1
        group_name = "{0}_{1}".format(pid, j)
        group_log = inStrain.logUtils.get_group_log('MergeProfile', group_name, 'start')

        if not single_thread:
            cmds = sprofile_cmd_queue.get(True)
        else:
            try:
                cmds = sprofile_cmd_queue.get(timeout=5)
            except:
                return

        last_c = len(cmds) - 1
        result_list = []
        for c, cmd in enumerate(cmds):
            Sprofile = cmd
            Sprofile.null_model = null_model

            for i in range(Sprofile.number_splits):
                Sprofile = Sprofile.update_splits(i, Sprofile_dict.pop(
                                                Sprofile.scaffold + '.' + str(i)))

            Sprofile = Sprofile.merge()

            if Sprofile is not None:
                if Sprofile.profile_genes:
                    try:
                        Sprofile.run_profile_genes()
                    except:
                        t = time.strftime('%m-%d %H:%M')
                        Sprofile.merge_log += "\n{1} DEBUG FAILURE GeneException {0}".format(str(Sprofile.scaffold), t)

            if c == last_c:
                try:
                    group_log += inStrain.logUtils.get_group_log('MergeProfile', group_name, 'end')
                    Sprofile.merge_log += "\n{0}".format(group_log)
                except:
                    pass
            result_list.append(Sprofile)

        Sprofiles.put(result_list)

def merge_basewise(mm2array_list):
    NAME2TYPE = {'coverage':'int32', 'clonality':'float32', 'snpCounted':'bool'}

    mm2bases = {}

    mms = set()
    for mm2bases in mm2array_list:
        mms = mms.union(set(mm2bases.keys()))

    for mm in mms:
        mm2bases[mm] = pd.concat([mm2array[mm] for mm2array in mm2array_list if mm in mm2array], verify_integrity=True)
    return mm2bases

def merge_special(vals, thing):
    #TODO ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts']
    return True

def make_coverage_table(covT, clonT, clonTR, lengt, scaff, SNPTable, null_model,
                        min_freq=0.05, debug=False):
    '''
    Add information to the table
    Args:
        covT: list of coverage values
        lengt: length of scaffold
        scaff: name of scaffold
        Wdb: windows to process over (NOT YET SUPPORTED)
    '''
    table = defaultdict(list)
    for mm in sorted(list(covT.keys())):

        covs = mm_counts_to_counts_shrunk(covT, mm, fill_zeros=int(lengt))

        if len(covs) == 0:
            covs = pd.Series([0]*lengt)

        nonzeros = np.count_nonzero(covs)

        # Get clonalities
        clons = get_basewise_clons(clonT, mm)
        Rclons = get_basewise_clons(clonTR, mm)

        counted_bases = len(clons)
        rarefied_bases = len(Rclons)
        SNS_count, SNV_count, div_site_count, con_snps, pop_snps = \
                    inStrain.profile.snv_utilities.calc_snps(
                                        SNPTable, mm, null_model,
                                        min_freq=min_freq)

        assert len(covs) == lengt, [covs, lengt, mm]

        # fill in all coverage information
        table['scaffold'].append(scaff)
        table['length'].append(lengt)
        table['breadth'].append(nonzeros/lengt)
        table['coverage'].append(np.mean(covs))
        table['coverage_median'].append(int(np.median(covs)))
        table['coverage_std'].append(np.std(covs))
        table['coverage_SEM'].append(scipy.stats.sem(covs))

        if len(clons) > 0:
            mean_c = np.mean(clons)
            median_c = np.median(clons)
            table['nucl_diversity'].append(1-mean_c)
            table['nucl_diversity_median'].append(1-median_c)
        else:
            table['nucl_diversity'].append(np.nan)
            table['nucl_diversity_median'].append(np.nan)

        if len(Rclons) > 0:
            mean_c = np.mean(Rclons)
            median_c = np.median(Rclons)
            table['nucl_diversity_rarefied'].append(1-mean_c)
            table['nucl_diversity_rarefied_median'].append(1-median_c)
        else:
            table['nucl_diversity_rarefied'].append(np.nan)
            table['nucl_diversity_rarefied_median'].append(np.nan)

        table['breadth_minCov'].append(counted_bases / lengt)
        table['breadth_rarefied'].append(rarefied_bases / lengt)
        table['breadth_expected'].append(estimate_breadth(table['coverage'][-1]))

        table['divergent_site_count'].append(div_site_count) # divergent_sites

        table['SNS_count'].append(SNS_count) # SNS_count
        table['SNV_count'].append(SNV_count)

        table['consensus_divergent_sites'].append(con_snps)
        table['population_divergent_sites'].append(pop_snps)

        if counted_bases == 0:
            table['conANI_reference'].append(0)
            table['popANI_reference'].append(0)
        else:
            table['conANI_reference'].append((counted_bases - con_snps)/ counted_bases)
            table['popANI_reference'].append((counted_bases - pop_snps)/ counted_bases)

        table['mm'].append(mm)

    return pd.DataFrame(table)

def mm_counts_to_counts_shrunk(MMcounts, maxMM=100, fill_zeros=False):
    '''
    Converts shunken mm-level items (e.g. covT) to a 1-d array

    Arguments:
        MMcounts: covT or similar object for a single scaffold
        maxMM: The maximum mm level to be allowed into the final array
        fill_zeros: If an int, fill 0s to make the array this length

    Returns:
        A 1-d array of counts

    '''
    if fill_zeros:
        counts = pd.Series(index=np.arange(fill_zeros))
    else:
        counts = pd.Series()

    for mm, count in [(mm, count) for mm, count in MMcounts.items() if mm <= maxMM]:
        counts = counts.add(count, fill_value=0)

    if fill_zeros:
        counts = counts.fillna(0)

    return counts

def get_basewise_clons(clonT, MM, fill_zeros=False):
    p2c = {}
    mms = sorted([int(mm) for mm in list(clonT.keys()) if int(mm) <= int(MM)])
    for mm in mms:
        p2c.update(clonT[mm].to_dict())

    counts = list(p2c.values())

    if fill_zeros:
        counts = counts.append(pd.Series(np.zeros(fill_zeros - len(counts))))

    return counts

def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage

    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    import numpy as np
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1

def _clonT_to_table(clonT):
    '''
    Based on the new clonT

    You want this to have "clonality, mm, and position"

    Called by SNVProfile
    '''
    dbs = []
    for mm, d in clonT.items():
        db = d.to_frame(name='clonality')#pd.DataFrame(d, columns=['clonality'])
        db['mm'] = mm
        db = db.reset_index()
        db = db.rename(columns={'index':'position'})
        dbs.append(db)

    if len(dbs) > 0:
        return pd.concat(dbs)

    else:
        return pd.DataFrame()

def _make_snp_table(Stable):
    '''
    Called by class scaffold_profile().make_cumulative_tables().
    Run on raw_snp_table.
    Very unclear what it actually does
    '''
    if Stable is not False:
        try:
            Sdb = pd.DataFrame(Stable)
            Sdb['scaffold'] = Sdb['scaffold'].astype('category')
            Sdb['con_base'] = Sdb['con_base'].astype('category')
        except KeyError:
            #logging.info("No SNPs detected!")
            Sdb = pd.DataFrame()
    else:
        Sdb = pd.DataFrame()

    return Sdb

def _parse_Sdb(sdb):
    '''
    Add some information to sdb
    '''
    if len(sdb) == 0:
        return sdb

    sdb['var_freq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['var_base'], sdb['position_coverage'])]
    sdb['con_freq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['con_base'], sdb['position_coverage'])]
    sdb['ref_freq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s if v in ['A','C','T','G'] else np.nan for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['ref_base'], sdb['position_coverage'])]

    return sdb

def gen_snv_profile(Sprofiles, ISP_loc=None, **kwargs):
    '''
    Take a bunch of scaffold profiles and make an SNVprofile
    '''
    if ISP_loc == None:
        location = kwargs.get('output', False)
    else:
        location = ISP_loc
        
    store_everything = kwargs.get('store_everything', False)

    # Merge things
    scaffold_list = []
    bam_list = []

    raw_snp_dbs = []
    raw_link_dbs = []

    cumu_scaff_dbs = []
    cumu_snv_dbs = []

    scaffold_2_mm_2_read_2_snvs = {}

    gene_level_results = defaultdict(list)

    for Sprof in Sprofiles:

        scaffold_list.append(Sprof.scaffold)
        bam_list.append(Sprof.bam)

        raw_snp_dbs.append(Sprof.raw_snp_table)
        raw_link_dbs.append(Sprof.raw_linkage_table)

        cumu_scaff_dbs.append(Sprof.cumulative_scaffold_table)
        cumu_snv_dbs.append(Sprof.cumulative_snv_table)

        if hasattr(Sprof, 'mm_reads_to_snvs'):
            scaffold_2_mm_2_read_2_snvs[Sprof.scaffold] = Sprof.mm_reads_to_snvs

        if hasattr(Sprof, 'gene_results'):
            for i, name in enumerate(['coverage', 'clonality', 'SNP_density', 'SNP_mutation_types']):
                gene_level_results[name].append(Sprof.gene_results[i])

    # Make some dataframes
    raw_snp_table = pd.concat(raw_snp_dbs).reset_index(drop=True)
    if len(raw_snp_table) > 0:
        for col in ['A', 'C', 'G', 'T', 'mm', 'position']:
            raw_snp_table[col] = raw_snp_table[col].astype(int)
            raw_snp_table['scaffold'] = raw_snp_table['scaffold'].astype('category')
            raw_snp_table['con_base'] = raw_snp_table['con_base'].astype('category')


    # Make the profile
    Sprofile = inStrain.SNVprofile.SNVprofile(location)

    # Add some things
    Sprofile.store('object_type', 'profile', 'value', 'Type of SNVprofile (profile or compare)')
    Sprofile.store('bam_loc', bam_list[0], 'value', 'Location of .bam file')
    Sprofile.store('scaffold_list', scaffold_list, 'list', '1d list of scaffolds that were profiled')
    Sprofile.store('raw_linkage_table', pd.concat(raw_link_dbs).reset_index(drop=True),
                    'pandas', 'Raw table of linkage information')
    Sprofile.store('raw_snp_table', raw_snp_table,
                    'pandas', 'Contains raw SNP information on a mm level')
    Sprofile.store('cumulative_scaffold_table', pd.concat(cumu_scaff_dbs).reset_index(drop=True),
                    'pandas', 'Cumulative coverage on mm level. Formerly scaffoldTable.csv')
    Sprofile.store('cumulative_snv_table', pd.concat(cumu_snv_dbs).reset_index(drop=True),
                    'pandas', 'Cumulative SNP on mm level. Formerly snpLocations.pickle')

    if scaffold_2_mm_2_read_2_snvs is not {}:
        Sprofile.store('scaffold_2_mm_2_read_2_snvs', scaffold_2_mm_2_read_2_snvs,
                        'pickle', 'crazy nonsense needed for linkage')

    if len(gene_level_results['coverage']) > 0:
        Sprofile.store('genes_coverage', pd.concat(gene_level_results['coverage']).reset_index(drop=True),
                 'pandas', 'Coverage of individual genes')
        Sprofile.store('genes_clonality', pd.concat(gene_level_results['clonality']).reset_index(drop=True),
                 'pandas', 'Clonality of individual genes')
        Sprofile.store('genes_SNP_count', pd.concat(gene_level_results['SNP_density']).reset_index(drop=True),
                 'pandas', 'SNP density and counts of individual genes')
        Sprofile.store('SNP_mutation_types', pd.concat(gene_level_results['SNP_mutation_types']).reset_index(drop=True),
                 'pandas', 'The mutation types of SNPs')

    # Store extra things
    att2descr = {'covT':'Scaffold -> mm -> position based coverage',
                 #'snpsCounted':"Scaffold -> mm -> position based True/False on if a SNPs is there",
                 'clonT':"Scaffold -> mm -> position based clonality",
                 'clonTR':"Scaffold -> mm -> rarefied position based clonality",
                 'read_to_snvs':'?',
                 'mm_to_position_graph':'?'}
    #for att in ['covT', 'snpsCounted', 'clonT', 'read_to_snvs', 'mm_to_position_graph']:
    for att in ['covT', 'clonT', 'clonTR', 'read_to_snvs', 'mm_to_position_graph']:
        if hasattr(Sprofiles[0], att):
            thing = {S.scaffold:getattr(S, att) for S in Sprofiles}
            Sprofile.store(att, thing, 'special', att2descr[att])

    if store_everything:
        counts_table = list() #alexcc 5/9/2019
        for Sprof in Sprofiles:
            counts_table.append(Sprof.pileup_counts) #add pileup counts to list
        counts_table = np.array(counts_table)
        Sprofile.store('counts_table', counts_table, 'numpy', '1d numpy array of 2D counts tables for each scaffold')

    return Sprofile

class ScaffoldSplitObject():
    '''
    This is a temporary object that holds all of the splits in a scaffold, and can eventually merge them into a scaffold_profile
    '''
    def __init__(self, number_splits):
        '''
        Initialize based on the number of expected splits to calculate
        '''
        self.number_splits = int(number_splits)
        self.split_dict = {i:None for i in range(number_splits)}
        assert self.number_splits == len(self.split_dict.keys())

    def ready(self):
        '''
        Return True if all splits are calculated
        '''
        calced = sum([True if v is not None else False for n, v in self.split_dict.items()])
        return calced == self.number_splits

    def print_status(self):
        '''
        Print a status report
        '''
        calced = sum([True if v is not None else False for n, v in self.split_dict.items()])
        print("{0} splits; {1} done".format(self.number_splits, calced))

    def update_splits(self, number, obj):
        '''
        Return a copy of the object with the split_dict updated
        '''
        self.split_dict[number] = obj
        return self

    def merge(self):
        '''
        Merge the split dict in order to make a scaffold profile
        '''
        log_message = inStrain.logUtils.get_worker_log('MergeProfile', self.scaffold, 'start')

        try:
            if self.number_splits == 1:
                Sprofile = self.split_dict[0].merge_single_profile(self)
                log_message += inStrain.logUtils.get_worker_log('MergeProfile', self.scaffold, 'end')
                Sprofile.merge_log = log_message
                return Sprofile

            else:
                Sprofile = scaffold_profile()
                Sprofile.null_model = self.null_model
                # Handle gene stuff
                for att in ['profile_genes', 'gene_database', 'gene2sequence']:
                    if hasattr(self, att):
                        setattr(Sprofile, att, getattr(self, att))

                # Handle value objects
                for att in ['scaffold', 'bam', 'min_freq']:
                    vals = set([getattr(Split, att) for num, Split in self.split_dict.items()])
                    assert len(vals) == 1, vals
                    setattr(Sprofile, att, list(vals)[0])

                # Handle sum objects
                for att in ['length']:
                    vals = [getattr(Split, att) for num, Split in self.split_dict.items()]
                    setattr(Sprofile, att, sum(vals))

                # Handle dataframes
                for att in ['raw_snp_table', 'raw_linkage_table']:
                    val = pd.concat([getattr(Split, att) for num, Split in self.split_dict.items()]).reset_index(drop=True)
                    setattr(Sprofile, att, val)

                # Handle series
                for att in ['covT', 'clonT', 'clonTR']:
                    vals = [getattr(Split, att) for num, Split in self.split_dict.items()]
                    setattr(Sprofile, att, merge_basewise(vals))

                # Handle extra things
                for att in ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts']:
                    if hasattr(self.split_dict[0], att):
                        vals = [getattr(Split, att) for num, Split in self.split_dict.items()]
                        setattr(Sprofile, att, merge_special(vals, att))

                Sprofile.make_cumulative_tables()

                log_message += inStrain.logUtils.get_worker_log('MergeProfile', self.scaffold, 'end')
                Sprofile.merge_log = log_message

                return Sprofile
        except:
            traceback.print_exc()

            t = time.strftime('%m-%d %H:%M')
            log_message = "\n{1} DEBUG FAILURE MergeError {0}\n"\
                            .format(self.scaffold, t)
            print(log_message)

            return None

    def delete_self(self):
        '''
        Delete everything to free up RAM
        '''
        pass


class SplitObject():
    '''
    Holds the profile of an individual split
    '''
    def __init__(self):
        pass

    def merge_single_profile(self, ScaffoldSplitObject):
        '''
        Convert self into scaffold_profile object
        '''
        Sprofile = scaffold_profile()
        Sprofile.scaffold = self.scaffold
        Sprofile.null_model = ScaffoldSplitObject.null_model
        Sprofile.bam = self.bam
        Sprofile.length = self.length
        Sprofile.raw_snp_table = self.raw_snp_table
        Sprofile.raw_linkage_table = self.raw_linkage_table
        Sprofile.covT = self.covT
        Sprofile.clonT = self.clonT
        Sprofile.clonTR = self.clonTR
        Sprofile.min_freq = self.min_freq

        for att in ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts'] + ['profile_genes', 'gene_database', 'gene2sequence']:
            if hasattr(self, att):
                setattr(Sprofile, att, getattr(self, att))

        for att in ['profile_genes', 'gene_database', 'gene2sequence']:
            if hasattr(ScaffoldSplitObject, att):
                setattr(Sprofile, att, getattr(ScaffoldSplitObject, att))

        # Make cummulative tables
        Sprofile.make_cumulative_tables()

        # Return
        return Sprofile

class scaffold_profile():
    '''
    This is a temporary object that represents the profile of a single scaffold
    '''
    def __init__(self, **kwargs):
        '''
        initialize
        '''
        self.version = __version__

    def make_cumulative_tables(self):
        '''
        Make cumulative tables (still raw-looking)
        This is all on the scaffold level
        '''
        if (self.raw_snp_table is not None):
            self.cumulative_snv_table = _make_snp_table(self.raw_snp_table)
            self.cumulative_snv_table = _parse_Sdb(self.cumulative_snv_table)

        self.cumulative_scaffold_table = make_coverage_table(self.covT, self.clonT,
                            self.clonTR, self.length, self.scaffold, self.raw_snp_table,
                            self.null_model, min_freq=self.min_freq)

    def run_profile_genes(self):
        """
        Run gene-level profile
        """
        if hasattr(self, 'gene_database'):

            scaffold = self.scaffold
            Gdb = self.gene_database
            covT = self.covT
            clonT = self.clonT
            gene2sequence = self.gene2sequence
            cumulative_snv_table = self.cumulative_snv_table

            self.gene_results = inStrain.GeneProfile.profile_genes_from_profile(scaffold, Gdb, covT, clonT, cumulative_snv_table, gene2sequence)

def _dlist():
    return defaultdict(list)

def _4zeros():
    return np.zeros(4, dtype=int)
