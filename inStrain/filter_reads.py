#!/usr/bin/env python

import copy
import gzip
import logging
import multiprocessing
import os
import random
import time
import traceback
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from numba import jit
from tqdm import tqdm

import inStrain.logUtils
import inStrain.profile.fasta

global i2o
global v2o

i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
v2o = {'min_read_ani':0, 'max_insert':1, 'min_insert':2, 'min_mapq':3}

class Controller():

    def main_from_profile(self, ISP, bam, **kwargs):
        '''
        The main method when called from the "profile" module

        Args:
            ISP = pre-initialized inStrain profile
            bam = location of .bam file
            args = the rest of the command line arguments

        Returns:
            Ridc = dictionary of read -> mismatches
            RR = pandas dataframe describing filtering
        '''
        detailed_report = kwargs.get('detailed_mapping_info', False)

        # Set up and parse .fasta file
        inStrain.logUtils.log_checkpoint("FilterReads", "load_fasta", "start")
        fasta_db, scaff2sequence, s2l = inStrain.profile.fasta.load_fasta(**kwargs)
        scaffolds = list(fasta_db['scaffold'].unique())
        inStrain.logUtils.log_checkpoint("FilterReads", "load_fasta", "end")

        # Filter the reads and store read reports
        if detailed_report:
            Rdic, RR, dRR = load_paired_reads(bam, scaffolds, **kwargs)

            # Store and delete the detailed report
            ISP.store('detailed_mapping_info', dRR, 'pandas', "Details report on reads")
            del dRR

        else:
            Rdic, RR = load_paired_reads(bam, scaffolds, **kwargs)

        # Return the Rdic and ReadReport
        return Rdic, RR, fasta_db, scaff2sequence, s2l

    def main(self, args):
        '''
        The main method when called explicitly (as its own module)
        '''
        bam = args.bam
        vargs = vars(args)
        del vargs['bam']

        detailed_report = vargs.get('detailed_mapping_info', False)
        generate_sam = vargs.get('generate_sam', False)
        out_folder = vargs.get('output', False)

        # Set up the output folder
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)

        # Set up .fasta file
        FAdb, s2s = load_fasta(args.fasta)

        # Get the paired reads
        scaffolds = list(FAdb['scaffold'].unique())
        if detailed_report:
            Rdic, RR, dRR = load_paired_reads(bam, scaffolds, **vargs)
        else:
            Rdic, RR = load_paired_reads(bam, scaffolds, **vargs)
            dRR = None

        # Make a .sam
        if generate_sam:
            print("The ability to make .sam files is not finished yet; sorry!")

        # Save results
        self.write_results(out_folder, RR, dRR, **vargs)

    def write_results(self, out_folder, RR, dRR, **kwargs):
        '''
        Save the results in a folder for the filter_reads module
        '''
        assert os.path.isdir(out_folder)

        RR_loc = os.path.join(out_folder, 'mapping_info.csv')
        write_mapping_info(RR, RR_loc, **kwargs)

        if dRR is not None:
            RR_loc = os.path.join(out_folder, 'detailed_mapping_info.csv')
            dRR.to_csv(RR_loc, index=False, sep='\t')

def read_profile_worker(read_cmd_queue, read_result_queue, bam, single_thread=False):
    '''
    Worker to filter reads
    '''
    # Initilize the .bam file
    bam_init = samfile = pysam.AlignmentFile(bam)

    while True:
        if not single_thread:
            cmds = read_cmd_queue.get(True)
        else:
            try:
                cmds = read_cmd_queue.get(timeout=5)
            except:
                return

        dicts, log = scaffold_profile_wrapper(cmds, bam_init)
        read_result_queue.put((dicts, log))

        # Clean up memory
        for d in dicts:
            del d
        del log
        del dicts

def load_paired_reads(bam, scaffolds, **kwargs):
    '''
    Load paired reads to be profiled

    You have this method do a lot of things because all of these things take lots of RAM, and you want them all to be cleared as soon as possible

    Return a dictionary of results. Some things that could be in it are:
        pair2infoF: A filtered dictionary of read name to number of mismatches
        RR: A summary read reaport
        RR_detailed: A detailed read report
    '''
    # Parse the kwargs
    detailed_report = kwargs.get('detailed_mapping_info', False)
    priority_reads_loc = kwargs.get('priority_reads', None)

    # Establish tallys to keep track of numbers
    tallys = {}

    # Get the pairs
    inStrain.logUtils.log_checkpoint("FilterReads", "get_paired_reads_multi", "start")
    scaff2pair2info = get_paired_reads_multi(bam, scaffolds, **kwargs)
    inStrain.logUtils.log_checkpoint("FilterReads", "get_paired_reads_multi", "end")

    if detailed_report:
        dRR = make_detailed_mapping_info(scaff2pair2info)

    # Handle paired-read filtering
    inStrain.logUtils.log_checkpoint("FilterReads", "paired_reads", "start")
    priority_reads = load_priority_reads(priority_reads_loc)
    scaff2pair2info = paired_read_filter(scaff2pair2info, priority_reads_set=priority_reads, tallys=tallys, **kwargs)
    inStrain.logUtils.log_checkpoint("FilterReads", "paired_reads", "end")

    # Filter and make the report
    inStrain.logUtils.log_checkpoint("FilterReads", "filter_reads", "start")
    scaff2pair2infoF, RR = filter_scaff2pair2info(scaff2pair2info, tallys,
                                                priority_reads_set=priority_reads,
                                                **kwargs)
    inStrain.logUtils.log_checkpoint("FilterReads", "filter_reads", "end")

    if detailed_report:
        return scaff2pair2infoF, RR, dRR
    else:
        return scaff2pair2infoF, RR

def filter_scaff2pair2info(scaff2pair2info, tallys={}, priority_reads_set=set(), **kwargs):
    '''
    Filter scaff2pair2info and generate a read report
    '''
    # Set up priority reads
    assert type(kwargs.get('priority_reads', 'none')) != type(set())
    priority_reads = priority_reads_set
    assert type(priority_reads) == type(set()), type(priority_reads)

    #item2order, to make it easier to read
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Calculate max insert
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    median_insert = np.median([value[i2o['insert_distance']] for scaff, pair2info \
                    in scaff2pair2info.items() for pair, value in pair2info.items()\
                    if value[i2o['reads']] == 2])
    max_insert = median_insert * max_insert_relative

    # Get filter values
    values = {}
    values['min_mapq'] = kwargs.get('min_mapq', 2)
    values['max_insert'] = max_insert
    values['min_insert'] = kwargs.get('min_insert', 50)
    values['min_read_ani'] = kwargs.get('min_read_ani', 0.97)
    values['pairing_filter'] = kwargs.get('pairing_filter', 'paired_only')

    # Set up the filtered dictionary
    scaff2pair2mm = {}

    # Make tallys for individual scaffolds
    table = defaultdict(list)

    for scaff, pair2info in scaff2pair2info.items():
        # Do the tallys
        if scaff not in tallys:
            tallys[scaff] = defaultdict(int)

        # Initialize some columns
        for c in ["pass_pairing_filter", "pass_min_read_ani", "pass_max_insert",
                    "pass_min_insert", "pass_min_mapq", "filtered_pairs",
                    "filtered_singletons", "filtered_priority_reads"]:
            tallys[scaff][c] = 0

        scaff2pair2mm[scaff] = {}
        for pair, info in pair2info.items():
            update_tallys(tallys, pair, info, values, scaff, scaff2pair2mm, priority_reads)

        # Make into a table
        table['scaffold'].append(scaff)
        for key, value in tallys[scaff].items():
            table[key].append(value)


        if len(pair2info.keys()) > 0:
            # Do the means
            for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
                table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
            table['mean_PID'].append(np.mean([(1 - (float(info[i2o['nm']]) /
                                                    float(info[i2o['length']]))) for pair, info in pair2info.items()]))

            # Do a the medians
            table['median_insert'].append(np.median([info[i2o['insert_distance']] for pair, info in pair2info.items()]))
        else:
            for att in ['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']:
                table['mean_' + att].append(np.nan)
            table['mean_PID'].append(np.nan)
            table['median_insert'].append(np.nan)

    try:
        Adb = pd.DataFrame(table)
    except:
        for k, v in table.items():
            print(k, len(v))
        assert False

    # Make tallys for all scaffolds
    table = defaultdict(list)
    table['scaffold'].append('all_scaffolds')

    CAdb = Adb[Adb['pass_pairing_filter'] > 0]
    total_reads = CAdb['pass_pairing_filter'].sum()

    for c in list(Adb.columns):
        if c == 'scaffold':
            pass
        elif c.startswith('mean_'):
            table[c].append(sum([v * m for v, m in zip(CAdb[c],
                                                       CAdb['pass_pairing_filter'])])/total_reads)
        elif c.startswith('median_'):
            table[c].append(sum([v * m for v, m in zip(CAdb[c],
                                                       CAdb['pass_pairing_filter'])])/total_reads)
        else:
            table[c].append(int(CAdb[c].sum()))
    adb = pd.DataFrame(table)

    # Concat
    Rdb = pd.concat([adb, Adb]).reset_index(drop=True)

    return  scaff2pair2mm, Rdb

# def update_tallys(tallys, pair, info, values, scaffold, scaff2pair2mm, priority_reads):
#     '''
#     The meat of filter_scaff2pair2info
#     '''
#     # Evaluate this pair
#     tallys[scaffold]['pass_pairing_filter'] += 1
#     f_results = evaluate_pair(info, values)
#
#     # Tally the results for what filteres passed
#     for name, index in v2o.items():
#          tallys[scaffold]['pass_' + name] += f_results[index]
#
#     # Tally the results for if the whole pair passed
#     if f_results.sum() == 4:
#         tallys[scaffold]['filtered_pairs'] += 1
#         scaff2pair2mm[scaffold][pair] = info[0]
#
#         if info[i2o['reads']] == 1:
#             tallys[scaffold]['filtered_singletons'] += 1
#
#         if pair in priority_reads:
#             tallys[scaffold]['filtered_priority_reads'] += 1
#
#
# i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
# v2o = {'min_read_ani':0, 'max_insert':1, 'min_insert':2, 'min_mapq':3}
# def evaluate_pair(info, values):
#     '''
#     Return a list of the filters that this pair passes and fails
#     Argumnets:
#         info: np array listing info about this pair in the i2o order
#         values: dictionary listing the filters to use when evaluting this pair
#     Returns:
#         f_resilts: np array listing which filters pass (1) and fail (0) in v2o order
#     '''
#     # Initialize results for this pair
#     f_results = np.zeros(4)
#
#     # Handle PID
#     PID = 1 - (float(info[i2o['nm']]) / float(info[i2o['length']]))
#     if PID > values['min_read_ani']:
#         f_results[v2o['min_read_ani']] = 1
#
#     # Handle mapQ
#     if info[i2o['mapq']] > values['min_mapq']:
#         f_results[v2o['min_mapq']] = 1
#
#     # If this is a pair check insert distance:
#     if ((info[i2o['reads']] == 2) & (info[i2o['insert_distance']] != -1)):
#         if info[i2o['insert_distance']] > values['min_insert']:
#             f_results[v2o['min_insert']] = 1
#         if info[i2o['insert_distance']] < values['max_insert']:
#             f_results[v2o['max_insert']] = 1
#
#     # Otherwise give those a pass
#     else:
#         f_results[v2o['min_insert']] = 1
#         f_results[v2o['max_insert']] = 1
#
#     return f_results

def update_tallys(tallys, pair, info, values, scaffold, scaff2pair2mm, priority_reads):
    '''
    The meat of filter_scaff2pair2info
    '''
    # Evaluate this pair
    tallys[scaffold]['pass_pairing_filter'] += 1
    f_results = evaluate_pair(info, np.zeros(4), values['min_read_ani'], values['min_mapq'], values['min_insert'],
                              values['max_insert'])

    # Tally the results for what filteres passed
    for name, index in v2o.items():
         tallys[scaffold]['pass_' + name] += f_results[index]

    # Tally the results for if the whole pair passed
    if f_results.sum() == 4:
        tallys[scaffold]['filtered_pairs'] += 1
        scaff2pair2mm[scaffold][pair] = info[0]

        if info[i2o['reads']] == 1:
            tallys[scaffold]['filtered_singletons'] += 1

        if pair in priority_reads:
            tallys[scaffold]['filtered_priority_reads'] += 1

@jit(nopython=True)
def evaluate_pair(info, f_results, min_read_ani, min_mapq, min_insert, max_insert):
    '''
    Return a list of the filters that this pair passes and fails

    Argumnets:
        info: np array listing info about this pair in the i2o order
        values: dictionary listing the filters to use when evaluting this pair

    Returns:
        f_resilts: np array listing which filters pass (1) and fail (0) in v2o order

    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
    v2o = {'min_read_ani':0, 'max_insert':1, 'min_insert':2, 'min_mapq':3}
    '''
    # Initialize results for this pair
    #f_results = np.zeros(4)

    # Handle PID
    PID = 1 - (float(info[0]) / float(info[3]))
    if PID > min_read_ani:
        f_results[0] = 1

    # Handle mapQ
    if info[2] > min_mapq:
        f_results[3] = 1

    # If this is a pair check insert distance:
    if ((info[4] == 2) & (info[1] != -1)):
        if info[1] > min_insert:
            f_results[2] = 1
        if info[1] < max_insert:
            f_results[1] = 1

    # Otherwise give those a pass
    else:
        f_results[1] = 1
        f_results[2] = 1

    return f_results

def load_priority_reads(file_loc):
    '''
    Loads a file of reads and returns a set of their names
    '''
    # is it None?
    if file_loc is None:
        return set()

    # Is it zipped?
    if file_loc[-3:] == '.gz':
        reader = gzip.open(file_loc, 'rt')
    else:
        reader = open(file_loc, 'r')

    # Figure out the type
    for line in reader.readlines():
        if line[0] == '@':
            TYPE = 'fastq'
        else:
            TYPE = 'list'
        break

    reader.close()
    if file_loc[-3:] == '.gz':
        reader = gzip.open(file_loc, 'rt')
    else:
        reader = open(file_loc, 'r')

    reads = set()
    if TYPE == 'fastq':
        for line in reader.readlines():
            if line[0] != '@':
                continue
            reads.add(line[1:].strip())

    elif TYPE == 'list':
        for line in reader.readlines():
            reads.add(line.strip())

    reader.close()

    return reads

def paired_read_filter(scaff2pair2info, priority_reads_set=set(), tallys=None, **kwargs):
    '''
    Filter scaff2pair2info to keep / remove paired / unpaired reads
    '''
    assert type(kwargs.get('priority_reads', 'none')) != type(set())
    priority_reads = priority_reads_set
    pairing_filter = kwargs.get('pairing_filter', 'paired_only')
    scaff2pair2infoF = {}
    pair2scaffold = {}
    assert type(priority_reads) == type(set()), type(priority_reads)

    for scaff, p2i in scaff2pair2info.items():
        # Initilize this scaffold
        scaff2pair2infoF[scaff] = {}
        if tallys is not None:
            tallys[scaff] = defaultdict(int)
            for v in ['unfiltered_reads', 'unfiltered_pairs', 'unfiltered_singletons',
                      'unfiltered_priority_reads']:
                tallys[scaff][v] = 0

        for p, i in p2i.items():
            # Update tallys; info[4] = number of reads
            if tallys is not None:
                tallys[scaff]['unfiltered_reads'] += i[4]
                if i[4] == 2:
                    tallys[scaff]['unfiltered_pairs'] += 1
                if i[4] == 1:
                    tallys[scaff]['unfiltered_singletons'] += 1
                if p in priority_reads:
                    tallys[scaff]['unfiltered_priority_reads'] += 1

            # Determine if it's going to make it into the final set
            if pairing_filter == 'paired_only':
                if ((i[4] == 2) | (p in priority_reads)):
                    scaff2pair2infoF[scaff][p] = i

            elif pairing_filter == 'non_discordant':
                # Add it if it's not already in there
                if ((p not in pair2scaffold) | (p in priority_reads)):
                    scaff2pair2infoF[scaff][p] = i
                    pair2scaffold[p] = scaff

                # If it is already in there, that means its concordant, so delete it
                else:
                    del scaff2pair2infoF[pair2scaffold[p]][p]

            elif pairing_filter == 'all_reads':
                if p in pair2scaffold:
                    # Average the pairs
                    mi = _merge_info(i, scaff2pair2infoF[pair2scaffold[p]][p])
                    scaff2pair2infoF[scaff][p] = mi
                    scaff2pair2infoF[pair2scaffold[p]][p] = mi

                else:
                    pair2scaffold[p] = scaff
                    scaff2pair2infoF[scaff][p] = i

            else:
                logging.error("Do not know paired read filter \"{0}\"; crashing now".format(pairing_filter))
                raise Exception

    return scaff2pair2infoF

def _merge_info(i1, i2):
    #{'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
    return np.array([i1[0] + i2[0],
            -2,
            max([i1[2] + i2[2]]),
            i1[3] + i2[3],
            i1[4] + i2[4],
            -1,
            -1], dtype="int64")

def make_detailed_mapping_info(scaff2pair2info, pairTOinfo=None, version=2):
    '''
    Make a detailed pandas dataframe from pair2info
    '''
    if pairTOinfo is None:
        pairTOinfo = dict()
    if version == 2:
        i2o = {'mm':0, 'insert_dist':1, 'mapq':2, 'length':3, 'reads':4,
                'start':5, 'stop':6}
    elif version == 1:
        i2o = {'mm':0, 'insert_dist':1, 'mapq':2, 'length':3,}

    keepers = pairTOinfo.keys()
    report_keepers =  (len(keepers) > 0)

    table = defaultdict(list)
    for scaff, pair2info in scaff2pair2info.items():
        for pair, array in pair2info.items():
            table['read_pair'].append(pair)
            table['scaffold'].append(scaff)
            if report_keepers:
                table['pass_filters'].append(pair in keepers)

            for item, location in i2o.items():
                table[item].append(array[location])

    return pd.DataFrame(table)

def load_fasta(fasta_file):
    '''
    Load the sequences to be profiled

    Return a table listing scaffold name, start, end
    '''
    # PROFILE ALL SCAFFOLDS IN THE .FASTA FILE
    scaff2sequence = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")) # set up .fasta file
    s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length
    Fdb = pd.DataFrame(list(s2l.items()), columns=['scaffold', 'end'])
    Fdb['start'] = 0

    return Fdb, scaff2sequence # also return s2l - alexcc 5/9/2019: Nah, make it scaff2sequence (s2s) (M.O. 6/10/19)

def filter_paired_reads_dict2(pair2info, **kwargs):
    '''
    Filter the dictionary of paired reads, end with read -> mm
    '''
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Get kwargs
    min_read_ani = kwargs.get('min_read_ani', 0.97)
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    min_insert = kwargs.get('min_insert', 50)
    min_mapq = kwargs.get('min_mapq', 2)

    # Get max insert
    max_insert = np.median([value[1] for key, value in pair2info.items() if value[i2o['reads']] == 2]) * max_insert_relative

    # Return dictionary of pairs
    return {copy.deepcopy(key):copy.deepcopy(value[0])
            for key, value in pair2info.items() if _evaluate_pair2(value,
            min_read_ani=min_read_ani,
            max_insert=max_insert,
            min_insert=min_insert,
            min_mapq=min_mapq)}

def makeFilterReport2(scaff2pair2info, pairTOinfo=False, priority_reads_set=None, **kwargs):
    '''
    Make a scaffold-level report on when reads are filtered using get_paired_reads_multi2

    If you've already filtered out pairs as you'd like, pass in pairTOinfo
    '''
    if priority_reads_set is None:
        priority_reads_set = set()
    assert type(kwargs.get('priority_reads', 'none')) != type(set())
    priority_reads = priority_reads_set
    profile_scaffolds = kwargs.get('scaffold_level_mapping_info', None)

    #item2order
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Calculate max insert
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    median_insert = np.median([value[i2o['insert_distance']] for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if value[i2o['reads']] == 2])
    max_insert = median_insert * max_insert_relative

    # Get values
    values = {}
    values['min_read_ani'] = kwargs.get('min_read_ani', 0.97)
    values['max_insert'] = max_insert
    values['min_insert'] = kwargs.get('min_insert', 50)
    values['min_mapq'] = kwargs.get('min_mapq', 2)

    # Make report on all scaffolds
    logging.debug('running on all reads')
    table = defaultdict(list)
    table['scaffold'].append('all_scaffolds')
    table['unfiltered_reads'].append(sum([value[i2o['reads']] for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items()]))
    table['unfiltered_pairs'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if value[i2o['reads']] == 2]))
    table['unfiltered_singletons'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if (value[i2o['reads']] == 1)]))
    table['unfiltered_priority_reads'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if (pair in priority_reads)]))

    if pairTOinfo != False:
        keepers = set(pairTOinfo.keys())
        infos = [info for scaff, pair2info in scaff2pair2info.items() for pair, info in pair2info.items() if pair in keepers]
        table['pass_pairing_filter'].append(len(infos))
    else:
        infos = [info for scaff, pair2info in scaff2pair2info.items() for pair, info in pair2info.items()]

    for att, v in values.items():
        kwargs={att:v}
        table['pass_' + att].append(len([True for info in infos if (_evaluate_pair2(info, **kwargs))]))
    table['filtered_pairs'].append(len([True for info in infos if (_evaluate_pair2(info, **values))]))
    table['filtered_singletons'].append(len([True for info in infos if ((info[i2o['reads']] == 1) & (_evaluate_pair2(info, **values)))]))
    table['filtered_priority_reads'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, info in pair2info.items() if ((pair in priority_reads) & (_evaluate_pair2(info, **values)))]))

    for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
        table['mean_' + att].append(np.mean([info[i] for info in infos]))
    table['median_insert'].append(np.median([info[i2o['insert_distance']] for info in infos]))
    table['mean_PID'].append(np.mean([(1 - (float(info[i2o['nm']]) / float(info[i2o['length']]))) for info in infos]))
    Adb = pd.DataFrame(table)

    table = defaultdict(list)
    logging.debug('running on individual scaffolds')
    for scaff, pair2info in scaff2pair2info.items():
        table['scaffold'].append(scaff)

        if pairTOinfo != False:
            #keepers = set(pairTOinfo.keys()) This is calculated above; dont need twice
            infos = [info for pair, info in pair2info.items() if pair in keepers]
            table['pass_pairing_filter'].append(len(infos))
        else:
            infos = [info for pair, info in pair2info.items()]

        table['filtered_pairs'].append(len([True for info in infos if (_evaluate_pair2(info, **values))]))

        if profile_scaffolds == True:
            table['unfiltered_reads'].append(sum([value[i2o['reads']] for pair, value in pair2info.items()]))
            table['unfiltered_pairs'].append(len([True for pair, value in pair2info.items() if value[i2o['reads']] == 2]))
            table['unfiltered_singletons'].append(len([True for pair, info in pair2info.items() if (info[i2o['reads']] == 1)]))
            table['unfiltered_priority_reads'].append(len([True for pair, info in pair2info.items() if (pair in priority_reads)]))

            for att, v in values.items():
                kwargs={att:v}
                table['pass_' + att].append(len([True for info in infos if (_evaluate_pair2(info, **kwargs))]))
            table['filtered_singletons'].append(len([True for info in infos if ((info[i2o['reads']] == 1) & (_evaluate_pair2(info, **values)))]))
            table['filtered_priority_reads'].append(len([True for pair, info in pair2info.items() if ((pair in priority_reads) & (_evaluate_pair2(info, **values)))]))

            for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
                table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
            table['median_insert'].append(np.median([value[1] for key, value in pair2info.items()]))
            table['mean_PID'].append(np.mean([(1 - (float(info[i2o['nm']]) / float(info[i2o['length']]))) for pair, info in pair2info.items()]))
    Sdb = pd.DataFrame(table)

    return pd.concat([Adb, Sdb]).reset_index(drop=True)

def write_mapping_info(RR, location, **kwargs):
    # Get header materials
    values = {}
    values['min_read_ani'] = kwargs.get('min_read_ani', 0.97)
    values['max_insert_relative'] = kwargs.get('max_insert_relative', 3)
    values['min_insert'] = kwargs.get('min_insert', 50)
    values['min_mapq'] = kwargs.get('min_mapq', 2)
    values['pairing_filter'] = kwargs.get('pairing_filter', 'paired_only')

    if location is None:
        return values

    # Write header
    os.remove(location) if os.path.exists(location) else None
    f = open(location, 'a')
    f.write("# {0}\n".format(' '.join(["{0}:{1}".format(k, v) for k, v in values.items()])))

    # Write csv
    RR.to_csv(f, index=False, sep='\t')

    # Close
    f.close()

def _evaluate_pair2(value, max_insert=1000000, min_read_ani=-1, min_insert=-1,
                        min_mapq=-1):
    # calculate PID for this pair
    PID = 1 - (float(value[i2o['nm']]) / float(value[i2o['length']]))

    # If this is a pair:
    if ((value[i2o['reads']] == 2) & (value[i2o['insert_distance']] != -1)):
        # See if it passes filtering
        if ((PID > min_read_ani) & (value[i2o['insert_distance']] > min_insert) & (value[i2o['insert_distance']] < max_insert)
                & (value[i2o['mapq']] > min_mapq)):
            return True
        else:
            return False

    # If this is not a pair
    else:
        if ((PID > min_read_ani) & (value[i2o['mapq']] > min_mapq)):
            return True
        else:
            return False

def get_paired_reads_multi(bam, scaffolds, **kwargs):
    '''
    Returns scaffold 2 read 2 info

    Modified version of (get_paired_reads_multi2)2
    '''
    # Initialize dictionary
    #dicts = []
    allPair2info = {}
    p = int(kwargs.get('processes', 6))
    ret_total = kwargs.get('ret_total', False)

    # Make commands
    cmd_groups = [x for x in iterate_command_groups(scaffolds, kwargs)]

    # Lets go with spawn; see if that reduces memory usage
    ctx = multiprocessing.get_context('spawn')

    # Make queues
    read_cmd_queue = ctx.Queue()
    read_result_queue = ctx.Queue()

    for cmd_group in cmd_groups:
        read_cmd_queue.put(cmd_group)

    inStrain.logUtils.log_checkpoint("FilterReads", "multiprocessing", "start")

    # Do the multiprocessing
    if p > 1:
        logging.debug("Establishing processes")
        processes = []
        for i in range(0, p):
            processes.append(ctx.Process(target=read_profile_worker, args=(read_cmd_queue, read_result_queue, bam)))
        for proc in processes:
            proc.start()

        # Set up progress bar
        pbar = tqdm(desc='Filtering reads: ', total=len(cmd_groups))

        # Get the results
        recieved_groups = 0
        while recieved_groups < len(cmd_groups):
            result_group = read_result_queue.get()
            recieved_groups += 1
            pbar.update(1)

            assert len(result_group) == 2
            results, log = result_group
            for s, pair2info in results:
                if pair2info != False:
                    allPair2info[s] = {}
                    for k, v in pair2info.items():
                        allPair2info[s][k] = v
                #dicts.append(result)
            logging.debug(log)

        # Close multi-processing
        for proc in processes:
            proc.terminate()

        # Close progress bar
        pbar.close()

    else:
        read_profile_worker(read_cmd_queue, read_result_queue, bam, single_thread=True)
        logging.info("Done profiling genes")

        # Get the results
        recieved_groups = 0
        while recieved_groups < len(cmd_groups):
            result_group = read_result_queue.get(timeout=5)
            recieved_groups += 1

            assert len(result_group) == 2
            results, log = result_group
            for s, pair2info in results:
                if pair2info != False:
                    allPair2info[s] = {}
                    for k, v in pair2info.items():
                        allPair2info[s][k] = v
                #dicts.append(result)
            logging.debug(log)

    inStrain.logUtils.log_checkpoint("FilterReads", "multiprocessing", "end")

    return allPair2info

class read_command():
    '''
    This class just holds all of the mumbo-jumbo needed to profile a scaffold
    '''
    def __init__(self):
        pass

def iterate_read_commands(scaffolds, bam, profArgs):
    '''
    Make and iterate profiling commands
    Doing it in this way makes it use way less RAM
    '''
    for scaff in scaffolds:
        yield [scaff, bam]

def iterate_command_groups(scaffolds, profArgs):
    '''
    Break these scaffolds into a series of scaffolds

    A command group is a list of scaffolds
    '''
    cmds = []
    number_groups = int(profArgs.get('ReadGroupSize', 5000))

    if number_groups > len(scaffolds):
        for scaff in scaffolds:
            yield [scaff]

    else:
        scaffs = random.sample(scaffolds, len(scaffolds))
        for n in range(number_groups):
            yield scaffs[n::number_groups]

def scaffold_profile_wrapper(cmd_group, bam_init):
    '''
    Take a command group and get the reads
    '''
    results = []
    log = ''

    scaffolds = cmd_group
    for scaff in scaffolds:
        try:
            pair2info, cur_log = get_paired_reads(bam_init, scaff)
            results.append([scaff, pair2info])
            log += cur_log

        except Exception as e:
            print(e)
            traceback.print_exc()
            logging.error("read filtering whole scaffold exception- {0}".format(scaff))
            results.append([scaff, None])

    return results, log

def get_paired_reads(samfile, scaff):
    '''
    Filter reads from a .bam file

    Rhis gets all reads and lets you filter out pairs (or not) later.

    Returns:
        pair2info: dictionary of read pair -> (mismatches, insert distance (-1 if number of reads is not 2), mapq score (highest), combined length, number of reads,
                                                reference_position_0, reference_position_1 (only used for small things))
    '''
    log_message = inStrain.logUtils.get_worker_log('GetPairedReads', scaff, 'start')

    # item2order, to make things more readable
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Initialize
    pair2info = {} # Information on pairs

    try:
        iter = samfile.fetch(scaff)
    except ValueError:
        t = time.strftime('%m-%d %H:%M')
        log_message = "\n{1} DEBUG FAILURE FilterReads {0} is not in .bam file\n".format(scaff, t)
        return {}, log_message

    for read in iter:
        # Dont use unmapped reads
        if read.get_reference_positions() == []:
            continue

        # Add a read the first time you see it
        elif read.query_name not in pair2info:
            pair2info[read.query_name] = np.array([read.get_tag('NM'),
                                            -1,
                                            read.mapping_quality,
                                            read.infer_query_length(),
                                            1,
                                            read.get_reference_positions()[0],
                                            read.get_reference_positions()[-1]], dtype="int64")

        # If we've seen this read's pair before
        elif read.query_name in pair2info:
            # Adjust the number of mismatches
            pair2info[read.query_name][i2o['nm']] = int(pair2info[read.query_name][i2o['nm']]) + int(read.get_tag('NM'))

            # Adjust the number of reads
            pair2info[read.query_name][i2o['reads']] += 1

            # Adjust the combined length
            pair2info[read.query_name][i2o['length']] += read.infer_query_length()

            # Adjust the mapQ
            pair2info[read.query_name][i2o['mapq']] = max(pair2info[read.query_name][i2o['mapq']], read.mapping_quality)

            # If the number of reads is 2, calculate insert size
            if pair2info[read.query_name][i2o['reads']] == 2:
                if read.get_reference_positions()[-1] > pair2info[read.query_name][i2o['start']]:
                    pair2info[read.query_name][i2o['insert_distance']] = read.get_reference_positions()[-1] - pair2info[read.query_name][i2o['start']]
                else:
                    pair2info[read.query_name][i2o['insert_distance']] = pair2info[read.query_name][i2o['stop']] - read.get_reference_positions()[0]

            # Otherwise get rid of insert distance
            else:
                pair2info[read.query_name][i2o['insert_distance']] = -1

            # Get rid of start and stop; not needed anymore
            pair2info[read.query_name][i2o['start']] = 0
            pair2info[read.query_name][i2o['stop']] = 0

    log_message += inStrain.logUtils.get_worker_log('GetPairedReads', scaff, 'end')

    return pair2info, log_message