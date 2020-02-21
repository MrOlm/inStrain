#!/usr/bin/env python

import os
import sys
import gzip
import pysam
import logging
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict
import concurrent.futures
from concurrent import futures
import traceback
import itertools
from collections import ChainMap

class Controller():

    def main(self, args):
        '''
        The main method
        '''
        # Set up .fasta file
        FAdb, s2s = load_fasta(args.fasta)

        # Make the read report
        if args.read_report != False:
            read_report_wrapper(args, FAdb)
            sys.exit()

        # inStrain.filter_reads.filter_reads(args.bam, positions[0], positions[1],
        #             args.filter_cutoff, args.max_insert_relative, args.min_insert,
        #             args.min_mapq, write_data = args.write, write_bam=args.generate_sam)

def load_paired_reads2(bam, scaffolds, **kwargs):
    '''
    Load paired reads to be profiled

    Return a dictionary of results. Some things that could be in it are:
        pair2infoF: A filtered dictionary of read name to number of mismatches
        RR: A summary read reaport
        RR_detailed: A detailed read report
    '''
    # Get the kwargs
    scaff2pair2info = get_paired_reads_multi2(bam, scaffolds, **kwargs)

    # Handle paired-read filtering
    priority_reads_loc = kwargs.get('priority_reads', None)
    priority_reads = load_priority_reads(priority_reads_loc)
    pair2info = paired_read_filter(scaff2pair2info, priority_reads_set=priority_reads, **kwargs)

    # Make read report
    logging.info('Making read report')
    RR = makeFilterReport2(scaff2pair2info, pairTOinfo=pair2info, priority_reads_set=priority_reads, **kwargs)

    # Filter the dictionary
    logging.info('Filtering reads')
    pair2infoF = filter_paired_reads_dict2(pair2info,
        **kwargs)

    return pair2infoF, RR

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

def paired_read_filter(scaff2pair2info, priority_reads_set=set(), **kwargs):
    '''
    Filter scaff2pair2info to keep / remove paired / unpaired reads
    '''
    assert type(kwargs.get('priority_reads', 'none')) != type(set())
    priority_reads = priority_reads_set
    pairing_filter = kwargs.get('pairing_filter', 'paired_only')
    pair2info = {}
    assert type(priority_reads) == type(set()), type(priority_reads)

    if pairing_filter == 'paired_only':
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                if ((i[4] == 2) | (p in priority_reads)):
                    pair2info[p] = i

    elif pairing_filter == 'non_discordant':
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                # Add it if it's not already in there
                if ((p not in pair2info) | (p in priority_reads)):
                    pair2info[p] = i

                # If it is already in there, that means its concordant, so delete it
                else:
                    del pair2info[p]

    elif pairing_filter == 'all':
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                if p in pair2info:
                    pair2info[p] = _merge_info(i, pair2info[p])
                else:
                    pair2info[p] = i

    return pair2info

def _merge_info(i1, i2):
    #{'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
    return np.array([i1[0] + i2[0],
            -2,
            max([i1[2] + i2[2]]),
            i1[3] + i2[3],
            i1[4] + i2[4],
            -1,
            -1], dtype="int64")

def make_detailed_read_report(scaff2pair2info, version=2):
    '''
    Make a detailed pandas dataframe from pair2info
    '''
    if version == 2:
        i2o = {'mm':0, 'insert_dist':1, 'mapq':2, 'length':3, 'reads':4,
                'start':5, 'stop':6}
    elif version == 1:
        i2o = {'mm':0, 'insert_dist':1, 'mapq':2, 'length':3,}

    table = defaultdict(list)
    for scaff, pair2info in scaff2pair2info.items():
        for pair, array in pair2info.items():
            table['read_pair'].append(pair)
            table['scaffold'].append(scaff)

            for item, location in i2o.items():
                table[item].append(array[location])

    return pd.DataFrame(table)

def get_fasta(fasta_file = None):
    positions = []
    total_length = 0
    for rec in SeqIO.parse(fasta_file, "fasta"):
        start = 0
        total_length += len(rec.seq)
        positions.append([str(rec.id),start, len(rec.seq)])

    return [positions, total_length]

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

def filter_paired_reads(bam, scaffolds, **kwargs):
    '''
    Filter reads according to specifications
    '''
    # Get kwargs
    filter_cutoff = kwargs.get('filter_cutoff', 0.97)
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    min_insert = kwargs.get('min_insert', 50)
    min_mapq = kwargs.get('min_mapq', 2)

    # Get paired reads
    pair2info = get_paired_reads(bam, scaffolds)

    # Get max insert
    max_insert = np.median([value[1] for key, value in pair2info.items()]) * max_insert_relative

    # Return names of pairs
    return [key for key, value in pair2info.items() if _evaluate_pair(value,
            filter_cutoff=filter_cutoff,
            max_insert=max_insert,
            min_insert=min_insert,
            min_mapq=min_mapq)]

def filter_paired_reads_dict(pair2info, **kwargs):
    '''
    Filter the dictionary of paired reads, end with read -> mm
    '''
    # Get kwargs
    filter_cutoff = kwargs.get('filter_cutoff', 0.97)
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    min_insert = kwargs.get('min_insert', 50)
    min_mapq = kwargs.get('min_mapq', 2)

    # Get max insert
    max_insert = np.median([value[1] for key, value in pair2info.items()]) * max_insert_relative

    # Return dictionary of pairs
    return {key:value[0] for key, value in pair2info.items() if _evaluate_pair(value,
            filter_cutoff=filter_cutoff,
            max_insert=max_insert,
            min_insert=min_insert,
            min_mapq=min_mapq)}

def filter_paired_reads_dict2(pair2info, **kwargs):
    '''
    Filter the dictionary of paired reads, end with read -> mm
    '''
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Get kwargs
    filter_cutoff = kwargs.get('filter_cutoff', 0.97)
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    min_insert = kwargs.get('min_insert', 50)
    min_mapq = kwargs.get('min_mapq', 2)

    # Get max insert
    max_insert = np.median([value[1] for key, value in pair2info.items() if value[i2o['reads']] == 2]) * max_insert_relative

    # Return dictionary of pairs
    return {key:value[0] for key, value in pair2info.items() if _evaluate_pair2(value,
            filter_cutoff=filter_cutoff,
            max_insert=max_insert,
            min_insert=min_insert,
            min_mapq=min_mapq)}

def makeFilterReport(scaff2pair2info, scaff2total, pairTOinfo=False, **kwargs):
    '''
    Make a scaffold-level report on when reads are filtered

    info = (mismatches, insert distance, mapq score, combined length)
    '''
    # Get kwargs
    max_insert_relative = kwargs.get('max_insert_relative', 3)

    # Get max insert
    median_insert = np.median(list(itertools.chain.from_iterable([[value[1] for key, value in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()])))
    max_insert = median_insert * max_insert_relative

    # Get values
    values = {}
    values['filter_cutoff'] = kwargs.get('filter_cutoff', 0.97)
    values['max_insert'] = max_insert
    values['min_insert'] = kwargs.get('min_insert', 50)
    values['min_mapq'] = kwargs.get('min_mapq', 2)

    # Make report on scaffolds
    logging.debug('running on all reads')
    table = defaultdict(list)
    table['scaffold'].append('all_scaffolds')
    table['unfiltered_reads'].append(sum([total for scaff, total in scaff2total.items()]))

    if pairTOinfo == False:
        pair2info = pairTOinfo
        table['unfiltered_pairs'].append(len(list(itertools.chain.from_iterable([[x for x in list(pair2info.keys())] for scaff, pair2info in scaff2pair2info.items()]))))
        for att, v in values.items():
            kwargs={att:v}
            table['pass_' + att].append(len(list(itertools.chain.from_iterable([[True for pair, info in pair2info.items() if (_evaluate_pair(info, **kwargs))] for scaff, pair2info in scaff2pair2info.items()]))))
        table['filtered_pairs'].append(len(list(itertools.chain.from_iterable([[True for pair, info in pair2info.items() if (_evaluate_pair(info, **values))] for scaff, pair2info in scaff2pair2info.items()]))))

        for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
            table['mean_' + att].append(np.mean(list(itertools.chain.from_iterable([[info[i] for pair, info in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()]))))
        table['median_insert'].append(median_insert)
        table['mean_PID'].append(np.mean(list(itertools.chain.from_iterable([[(1 - (float(info[0]) / float(info[3]))) for pair, info in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()]))))

    else:
        table['unfiltered_pairs'].append(len(pair2info.keys()))
        for att, v in values.items():
            kwargs={att:v}
            table['pass_' + att].append(len([True for pair, info in pair2info.items() if (_evaluate_pair(info, **kwargs))]))
        table['filtered_pairs'].append(len([True for pair, info in pair2info.items() if (_evaluate_pair(info, **values))]))

        for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
            table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
        table['median_insert'].append(np.median([value[1] for key, value in pair2info.items()]))
        table['mean_PID'].append(np.mean([(1 - (float(info[0]) / float(info[3]))) for pair, info in pair2info.items()]))

    logging.debug('running on individual scaffolds')
    for scaff, pair2info in scaff2pair2info.items():
        table['scaffold'].append(scaff)
        table['unfiltered_reads'].append(scaff2total[scaff])
        table['unfiltered_pairs'].append(len(pair2info.keys()))
        for att, v in values.items():
            kwargs={att:v}
            table['pass_' + att].append(len([True for pair, info in pair2info.items() if (_evaluate_pair(info, **kwargs))]))
        table['filtered_pairs'].append(len([True for pair, info in pair2info.items() if (_evaluate_pair(info, **values))]))

        for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
            table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
        table['median_insert'].append(np.median([value[1] for key, value in pair2info.items()]))
        table['mean_PID'].append(np.mean([(1 - (float(info[0]) / float(info[3]))) for pair, info in pair2info.items()]))

    return pd.DataFrame(table)

def makeFilterReport2(scaff2pair2info, pairTOinfo=False, priority_reads_set=set(), **kwargs):
    '''
    Make a scaffold-level report on when reads are filtered using get_paired_reads_multi2

    If you've already filtered out pairs as you'd like, pass in pairTOinfo
    '''
    assert type(kwargs.get('priority_reads', 'none')) != type(set())
    priority_reads = priority_reads_set

    #item2order
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Calculate max insert
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    median_insert = np.median([value[i2o['insert_distance']] for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if value[i2o['reads']] == 2])
    max_insert = median_insert * max_insert_relative

    # Get values
    values = {}
    values['filter_cutoff'] = kwargs.get('filter_cutoff', 0.97)
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

    # # You didn't do any other filtering of paired reads
    # if pairTOinfo == False:
    #     table['unfiltered_pairs'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if value[i2o['reads']] == 2]))
    #     for att, v in values.items():
    #         kwargs={att:v}
    #         table['pass_' + att].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if (_evaluate_pair2(value, **kwargs))]))
    #     table['filtered_pairs'].append(len([True for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if (_evaluate_pair2(value, **values))]))
    #
    #     # for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
    #     #     table['mean_' + att].append(np.mean([value[i2o['insert_distance']] for scaff, pair2info in scaff2pair2info.items() for pair, value in pair2info.items() if value[i2o['reads']] != 2]))
    #     #     #table['mean_' + att].append(np.mean(list(itertools.chain.from_iterable([[info[i] for pair, info in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()]))))
    #     # table['median_insert'].append(median_insert)
    #     # table['mean_PID'].append(np.mean(list(itertools.chain.from_iterable([[(1 - (float(info[0]) / float(info[3]))) for pair, info in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()]))))
    #
    # # You did do filtering of paired reads
    # else:
    #     pair2info = pairTOinfo
    #     table['unfiltered_pairs'].append(len(pair2info.keys()))
    #     table['unfiltered_singletons'].append(len([True for pair, info in pair2info.items() if (info[i2o['reads']] == 1)]))
    #     table['unfiltered_priority_reads'].append(len([True for pair, info in pair2info.items() if (pair in priority_reads)]))
    #     table['pass_pairing_filter'].append(len(pair2info.keys()))
    #     for att, v in values.items():
    #         kwargs={att:v}
    #         table['pass_' + att].append(len([True for pair, info in pair2info.items() if (_evaluate_pair2(info, **kwargs))]))
    #     table['filtered_pairs'].append(len([True for pair, info in pair2info.items() if (_evaluate_pair2(info, **values))]))
    #     table['filtered_singletons'].append(len([True for pair, info in pair2info.items() if ((info[i2o['reads']] == 1) & (_evaluate_pair2(info, **values)))]))
    #     table['filtered_priority_reads'].append(len([True for pair, info in pair2info.items() if ((pair in priority_reads) & (_evaluate_pair2(info, **values)))]))
    #
    #     for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
    #         table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
    #     table['median_insert'].append(np.median([value[1] for key, value in pair2info.items()]))
    #     table['mean_PID'].append(np.mean([(1 - (float(info[i2o['nm']]) / float(info[i2o['length']]))) for pair, info in pair2info.items()]))

    logging.debug('running on individual scaffolds')
    for scaff, pair2info in scaff2pair2info.items():
        table['scaffold'].append(scaff)
        table['unfiltered_reads'].append(sum([value[i2o['reads']] for pair, value in pair2info.items()]))
        table['unfiltered_pairs'].append(len([True for pair, value in pair2info.items() if value[i2o['reads']] == 2]))
        table['unfiltered_singletons'].append(len([True for pair, info in pair2info.items() if (info[i2o['reads']] == 1)]))
        table['unfiltered_priority_reads'].append(len([True for pair, info in pair2info.items() if (pair in priority_reads)]))

        if pairTOinfo != False:
            keepers = set(pairTOinfo.keys())
            infos = [info for pair, info in pair2info.items() if pair in keepers]
            table['pass_pairing_filter'].append(len(infos))
        else:
            infos = [info for pair, info in pair2info.items()]

        for att, v in values.items():
            kwargs={att:v}
            table['pass_' + att].append(len([True for info in infos if (_evaluate_pair2(info, **kwargs))]))
        table['filtered_pairs'].append(len([True for info in infos if (_evaluate_pair2(info, **values))]))
        table['filtered_singletons'].append(len([True for info in infos if ((info[i2o['reads']] == 1) & (_evaluate_pair2(info, **values)))]))
        table['filtered_priority_reads'].append(len([True for pair, info in pair2info.items() if ((pair in priority_reads) & (_evaluate_pair2(info, **values)))]))

        for i, att in enumerate(['mistmaches', 'insert_distance', 'mapq_score', 'pair_length']):
            table['mean_' + att].append(np.mean([info[i] for pair, info in pair2info.items()]))
        table['median_insert'].append(np.median([value[1] for key, value in pair2info.items()]))
        table['mean_PID'].append(np.mean([(1 - (float(info[i2o['nm']]) / float(info[i2o['length']]))) for pair, info in pair2info.items()]))

    return pd.DataFrame(table)

def write_read_report(RR, location, **kwargs):
    # Get header materials
    values = {}
    values['filter_cutoff'] = kwargs.get('filter_cutoff', 0.97)
    values['max_insert_relative'] = kwargs.get('max_insert_relative', 3)
    values['min_insert'] = kwargs.get('min_insert', 50)
    values['min_mapq'] = kwargs.get('min_mapq', 2)

    # Write header
    os.remove(location) if os.path.exists(location) else None
    f = open(location, 'a')
    f.write("# {0}\n".format(' '.join(["{0}:{1}".format(k, v) for k, v in values.items()])))

    # Write csv
    RR.to_csv(f, index=False, sep='\t')

    # Close
    f.close()

def filter_paired_reads_dict_scaff(scaff2pair2info, **kwargs):
    '''
    Filter the dictionary of paired reads, end with read -> mm
    '''
    # Get kwargs
    filter_cutoff = kwargs.get('filter_cutoff', 0.97)
    max_insert_relative = kwargs.get('max_insert_relative', 3)
    min_insert = kwargs.get('min_insert', 50)
    min_mapq = kwargs.get('min_mapq', 2)

    # Get max insert
    max_insert = np.median(list(itertools.chain.from_iterable([[value[1] for key, value in pair2info.items()] for scaff, pair2info in scaff2pair2info.items()]))) * max_insert_relative

    # Return dictionary of pairs
    ret = {}
    for scaff, pair2info in scaff2pair2info.items():
        for key, value in pair2info.items():
            if _evaluate_pair(value,
                filter_cutoff=filter_cutoff,
                max_insert=max_insert,
                min_insert=min_insert,
                min_mapq=min_mapq):
                ret[key] = value[0]
    return ret

def _evaluate_pair(value, max_insert=1000000, filter_cutoff=-1, min_insert=-1,
                                                        min_mapq=-1):
    # calculate PID for this pair
    PID = 1 - (float(value[0]) / float(value[3]))

    # See if it passes filtering
    if ((PID > filter_cutoff) & (value[1] > min_insert) & (value[1] < max_insert)\
        & (value[2] > min_mapq)):
        return True
    else:
        return False

i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}
def _evaluate_pair2(value, max_insert=1000000, filter_cutoff=-1, min_insert=-1,
                        min_mapq=-1):
    # calculate PID for this pair
    PID = 1 - (float(value[i2o['nm']]) / float(value[i2o['length']]))

    # If this is a pair:
    if ((value[i2o['reads']] == 2) & (value[i2o['insert_distance']] != -1)):
        # See if it passes filtering
        if ((PID > filter_cutoff) & (value[i2o['insert_distance']] > min_insert) & (value[i2o['insert_distance']] < max_insert)\
            & (value[i2o['mapq']] > min_mapq)):
            return True
        else:
            return False

    # If this is not a pair
    else:
        if ((PID > filter_cutoff) & (value[i2o['mapq']] > min_mapq)):
            return True
        else:
            return False

def get_paired_reads_multi(bam, scaffolds, **kwargs):
    '''
    Returns scaffold 2 read 2 info
    '''
    # Initialize dictionary
    dicts = []
    p = int(kwargs.get('processes', 6))
    ret_total = kwargs.get('ret_total', False)

    # Do the multiprocessing
    if p > 1:
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=p)

        total_cmds = len([x for x in iterate_read_commands(scaffolds, bam, kwargs)])

        wait_for = [executor.submit(scaffold_profile_wrapper, cmd) for cmd in iterate_read_commands(scaffolds, bam, kwargs)]

        for f in tqdm(futures.as_completed(wait_for), total=total_cmds, desc='Getting read pairs: '):
            try:
                results = f.result()
                dicts.append(results)
            except:
                logging.error("We had a failure! Not sure where!")

    else:
        for cmd in tqdm(iterate_read_commands(scaffolds, bam, kwargs), desc='Getting read pairs: ', total=len(scaffolds)):
            dicts.append(scaffold_profile_wrapper(cmd))

    # Try and save any ones that failed
    failed_scaffs = set(scaffolds) - set([dict[0] for dict in dicts])
    if len(failed_scaffs) > 0:
        logging.error("The following scaffolds failed- I'll try again {0}".format('\n'.join(failed_scaffs)))
        for cmd in tqdm(iterate_read_commands(failed_scaffs, bam, kwargs), desc='Getting read pairs for previously failed scaffolds: ', total=len(failed_scaffs)):
            dicts.append(scaffold_profile_wrapper(cmd))

    allPair2info = {}
    scaff2total = {}
    for s, d, t in dicts:
        if d != False:
            allPair2info[s] = {}
            for k, v in d.items():
                allPair2info[s][k] = v
            scaff2total[s] = t

    if ret_total:
        return allPair2info, scaff2total
    else:
        return allPair2info

def get_paired_reads_multi2(bam, scaffolds, **kwargs):
    '''
    Returns scaffold 2 read 2 info
    '''
    # Initialize dictionary
    dicts = []
    p = int(kwargs.get('processes', 6))
    ret_total = kwargs.get('ret_total', False)

    # Do the multiprocessing
    if p > 1:
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=p)

        total_cmds = len([x for x in iterate_read_commands(scaffolds, bam, kwargs)])

        wait_for = [executor.submit(scaffold_profile_wrapper2, cmd) for cmd in iterate_read_commands(scaffolds, bam, kwargs)]

        for f in tqdm(futures.as_completed(wait_for), total=total_cmds, desc='Getting read pairs: '):
            try:
                results = f.result()
                dicts.append(results)
            except:
                logging.error("We had a failure! Not sure where!")

    else:
        for cmd in tqdm(iterate_read_commands(scaffolds, bam, kwargs), desc='Getting read pairs: ', total=len(scaffolds)):
            dicts.append(scaffold_profile_wrapper2(cmd))

    # Try and save any ones that failed
    failed_scaffs = set(scaffolds) - set([dict[0] for dict in dicts])
    if len(failed_scaffs) > 0:
        logging.error("The following scaffolds failed- I'll try again {0}".format('\n'.join(failed_scaffs)))
        for cmd in tqdm(iterate_read_commands(failed_scaffs, bam, kwargs), desc='Getting read pairs for previously failed scaffolds: ', total=len(failed_scaffs)):
            dicts.append(scaffold_profile_wrapper2(cmd))

    allPair2info = {}
    for s, pair2info in dicts:
        if pair2info != False:
            allPair2info[s] = {}
            for k, v in pair2info.items():
                allPair2info[s][k] = v

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

def scaffold_profile_wrapper(cmd):
    '''
    Take a scaffold and get the reads

    cmd[0] = scaffold, cmd[1] = bam
    '''
    logging.debug('running {0} reads'.format(cmd[0]))
    try:
        x, y =  get_paired_reads(cmd[1], [cmd[0]], ret_total=True)
        logging.debug('returning {0} reads'.format(cmd[0]))
        return cmd[0], x, y
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(cmd[0])))
        return cmd[0], False

def scaffold_profile_wrapper2(cmd):
    '''
    Take a scaffold and get the reads

    cmd[0] = scaffold, cmd[1] = bam, cmd[2] = kwargs
    '''
    logging.debug('running {0} reads'.format(cmd[0]))
    try:
        pair2info =  get_paired_reads2(cmd[1], cmd[0])
        logging.debug('returning {0} reads'.format(cmd[0]))
        return cmd[0], pair2info
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(cmd[0])))
        return cmd[0], False

def get_paired_reads(bam, scaffolds, ret_total=False):
    '''
    Filter reads from a .bam file

    Edits 8.14.19 for huge RAM efficiency boost:
    - Making read_data not a nested dictinoary
    - Removing scaffold check (not necessary)

    Returns:
        pair2info - dictionary of read pair -> (mismatches, insert distance, mapq score, combined length)

    Documentation 2.11.20
    - read_data = dictionary of query_name -> np.array([read.get_tag('NM'),
                                    read.infer_query_length(),
                                    read.mapping_quality,
                                    read.get_reference_positions()[0],
                                    read.get_reference_positions()[-1]], dtype="int64")
    - the first time a name is encountered, it's added to read data
    '''
    #item2order
    i2o = {'nm':0, 'len':1, 'mapq':2, 'start':3, 'stop':4}

    # Initialize dictionaries
    pair2info = {} # Information on pairs
    total_unpaired = 0

    samfile = pysam.AlignmentFile(bam)
    for scaff in scaffolds:
    #for scaff in tqdm(scaffolds, desc='Getting read pairs: '):
        read_data = {} # Information on the first pair of each read

        try:
            iter = samfile.fetch(scaff)
        except ValueError:
            logging.error("{0} is not in .bam file".format(scaff))
            continue
            # if ret_total:
            #     return pair2info, total_unpaired
            # else:
            #     return pair2info

        for read in iter:
            total_unpaired += 1
            # If we've seen this read's pair before
            if read.query_name in read_data:
                # Make sure that the pair is on the same scaffold and that it's mapped at all
                if (read.get_reference_positions() != []):
                    # Add this read_pair to pair2info
                    pairMM = int(read_data[read.query_name][i2o['nm']]) + int(read.get_tag('NM'))

                    # Below is the old way
                    mapped_read_lengths = read_data[read.query_name][i2o['len']] \
                                            + read.infer_query_length()

                    # This is the new way
                    #mapped_read_lengths = len(set(read.get_reference_positions() + read_data[read.query_name]['ref_pos']))


                    if read.get_reference_positions()[-1] > read_data[read.query_name][i2o['start']]:
                        pair_inserts = read.get_reference_positions()[-1] - read_data[read.query_name][i2o['start']]
                    else:
                        pair_inserts = read_data[read.query_name][i2o['stop']] - read.get_reference_positions()[0]
                    pair_mapq = max(read.mapping_quality, read_data[read.query_name][i2o['mapq']])

                    pair2info[read.query_name] = (pairMM, pair_inserts, pair_mapq, mapped_read_lengths)

            # Add this in search of its mate
            elif read.get_reference_positions() != []: # don't use unmapped reads:
                read_data[read.query_name] = np.array([read.get_tag('NM'),
                                                read.infer_query_length(),
                                                read.mapping_quality,
                                                read.get_reference_positions()[0],
                                                read.get_reference_positions()[-1]], dtype="int64")
    if ret_total:
        return pair2info, total_unpaired

    else:
        return pair2info

def get_paired_reads2(bam, scaff):
    '''
    Filter reads from a .bam file

    As opposed to get_paired_reads, this gets all reads and lets you filter out pairs (or not) later. This just means adding a bit to pair2info

    Returns:
        pair2info: dictionary of read pair -> (mismatches, insert distance (-1 if number of reads is not 2), mapq score (highest), combined length, number of reads,
                                                reference_position_0, reference_position_1 (only used for small things))
    '''
    # item2order, to make things more readable
    i2o = {'nm':0, 'insert_distance':1, 'mapq':2, 'length':3, 'reads':4, 'start':5, 'stop':6}

    # Initialize
    pair2info = {} # Information on pairs
    samfile = pysam.AlignmentFile(bam)

    try:
        iter = samfile.fetch(scaff)
    except ValueError:
        logging.error("{0} is not in .bam file".format(scaff))
        return {}

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

    return pair2info

def read_report_wrapper(args, FAdb):
    # Get paired reads
    scaffolds = list(FAdb['scaffold'].tolist())
    scaff2pair2info, scaff2total = get_paired_reads_multi(args.bam, scaffolds, processes=args.processes, ret_total=True)

    # Make report
    RR = makeFilterReport(scaff2pair2info, scaff2total, **vars(args))

    # Save report
    write_read_report(RR, args.read_report, **vars(args))

def filter_reads(bam, positions, fasta_length, filter_cutoff = 0.97, max_insert_relative = 3, min_insert = 50, min_mapq = 2, write_data = None, write_bam = False):

    # read sets
    observed_read1s = set()
    observed_read2s = set()
    mapped_pairs = set()
    final_reads = set()

    # counters
    total_read_count = 0
    total_read_pairs = 0
    total_mapped_pairs = 0
    mapped_read_lengths = 0

    # storing data
    read_data = {}
    pair_mapqs = {}
    pair_mismatch = {}
    pair_inserts = {}

    samfile = pysam.AlignmentFile(bam)

    #for printing out a new bam file
    if write_bam:
        logging.info("Copying header for new bam...")
        samfile_out = pysam.AlignmentFile(bam.split("/")[-1].split(".")[0] + "_filtered.bam", "wb", template=samfile)
        reads_all = defaultdict(list)

    logging.info("READING BAM: " + bam.split("/")[-1])
    logging.info("Using reads with >" + str(filter_cutoff) + "% PID to consensus reference.")

    ## STEP 1: collect paired reads and their information
    for gene in tqdm(positions, desc='Getting read pairs: '):
        for read in samfile.fetch(gene[0], gene[1], gene[2]):
            total_read_count += 1

            #store all reads if we're going to write them back to a new bam file
            if write_bam:
                reads_all[read.query_name].append(read)

            ## If we've seen this read's pair before
            if (read.is_read2 and read.query_name in observed_read1s) or (read.is_read1 and read.query_name in observed_read2s):

                #But if we haven't already seen this complete pair, then we can complete the pair and store the information
                #Also check that the pair is on the same scaffold
                if read.query_name not in mapped_pairs and gene[0] == read_data[read.query_name]['scaf']:
                    total_read_pairs += 1

                    if read.get_reference_positions() != []:
                        total_mapped_pairs += 1
                        mapped_pairs.add(read.query_name) #add to found

                        #for calculating mean read length
                        mapped_read_lengths += float(read_data[read.query_name]['len'])

                        #set mismatch percentage
                        pair_mismatch[read.query_name] = 1- ( ( float(read_data[read.query_name]['nm']) + float(read.get_tag('NM')) ) / ( float(read_data[read.query_name]['len']) + read.infer_query_length()) )
                        #set insert size
                        if read.get_reference_positions()[-1] > read_data[read.query_name]['start']:
                            pair_inserts[read.query_name] = read.get_reference_positions()[-1] - read_data[read.query_name]['start']
                        else:
                            pair_inserts[read.query_name] = read_data[read.query_name]['stop'] - read.get_reference_positions()[0]
                        #set mapq
                        pair_mapqs[read.query_name] = read.mapping_quality
                        if read_data[read.query_name]['mapq'] > read.mapping_quality:
                            pair_mapqs[read.query_name] = read_data[read.query_name]['mapq']

            #this is the first time we see a read from this pair and don't double count
            elif (read.is_read1 and read.query_name not in observed_read1s) or (read.is_read2 and read.query_name not in observed_read2s):
                if read.get_reference_positions() != []: # don't use unmapped reads
                        if read.is_read1:
                            observed_read1s.add(read.query_name)
                        else:
                            observed_read2s.add(read.query_name)
                        #record the data for this read
                        read_data[read.query_name] = {"nm": read.get_tag('NM'), "len": read.infer_query_length(), "mapq": read.mapping_quality, "start": read.get_reference_positions()[0], 'stop': read.get_reference_positions()[-1], 'scaf': gene[0]}


    ## STEP 2: INSERT SIZE CUTOFF, MAPQ CUTOFF, AND MISMATCH CUTOFF
    mapped_read_lengths = mapped_read_lengths / total_mapped_pairs

    max_insert = np.median(list(pair_inserts.values())) * max_insert_relative #insert size should be less than max_insert_relative * median value 
    too_short = 0.0
    too_long = 0.0
    good_length = 0.0
    mapq_good = 0.0
    filter_cutoff_good = 0.0

    logging.info("Filtering reads...")

    for read_pair in mapped_pairs:
        if pair_inserts[read_pair] > min_insert:
            if pair_inserts[read_pair] < max_insert:
                good_length += 2
                if pair_mapqs[read_pair] > min_mapq:
                    mapq_good += 2

                    # Which set does this read go into?
                    if pair_mismatch[read_pair] > filter_cutoff:
                        filter_cutoff_good += 2
                        final_reads.add(read_pair)

                        #write out to new bam file if option selected
                        if write_bam:
                            for read in reads_all[read_pair]:
                                samfile_out.write(read)
            else:
                too_long += 2
        else:
            too_short += 2

    table = defaultdict(list)
    table["total reads found"].append(str(total_read_count))
    table["average mapped read length"].append(str(mapped_read_lengths))
    table["total fasta length"].append(str(fasta_length))
    table["expected possible coverage"].append(str(float(total_read_count)*mapped_read_lengths / fasta_length))
    table["total paired reads"].append(str(total_read_pairs*2))
    table["total paired reads (%)"].append(str(int(100*total_read_pairs*2.0 / total_read_count)))
    table["total same scaffold mapped paired reads"].append(str(total_mapped_pairs*2))
    table["total same scaffold mapped paired reads (%)"].append(str(int(100*total_read_pairs*2.0 / total_read_count)))
    table["median insert size"].append(str(max_insert / max_insert_relative))
    table["paired reads < 50 bp apart"].append(str(too_short))
    table["max insert"].append(str(max_insert))
    table["paired reads > max insert apart"].append(str(too_long))
    table["reads which also pass both pair insert size filters"].append(str(good_length))
    table["reads which also pass both pair insert size filters (%)"].append(str(int(100*float(good_length) / total_read_count)))
    table["minimum mapq threshold"].append(str(min_mapq))
    table["reads which pass minimum mapq threshold"].append(str(mapq_good))
    table["reads which pass minimum mapq threshold (%)"].append(str(int(100*float(mapq_good) / total_read_count)))
    table['minimum PID'].append(str(filter_cutoff))
    table["(final) reads which also pass read pair PID"].append(filter_cutoff_good)
    table["(final) reads which also pass read pair PID (%)"].append(str(int(100*float(filter_cutoff_good) / total_read_count)))
    table["(final) expected coverage"].append(str(float(filter_cutoff_good) * mapped_read_lengths / fasta_length))
    Rdb = pd.DataFrame(table)

    logging.debug("**READ STATSTICS**")
    logging.debug("total reads found: " + str(total_read_count))
    logging.debug("average mapped read length: " + str(mapped_read_lengths))
    logging.debug("total fasta length: " + str(fasta_length))
    logging.debug("expected possible coverage: " + str(float(total_read_count)*mapped_read_lengths / fasta_length))
    logging.debug("total paired reads: " + str(total_read_pairs*2) + " (" + str(int(100*total_read_pairs*2.0 / total_read_count)) + "%)")
    logging.debug("total same scaffold mapped paired reads: " + str(total_mapped_pairs*2) + " (" + str(int(100*total_mapped_pairs*2.0 / total_read_count)) + "%)")
    logging.debug("")
    logging.debug("median insert size: " + str(max_insert / max_insert_relative))
    logging.debug("paired reads < 50 bp apart: " + str(too_short))
    logging.debug("paired reads > " + str(max_insert) + " apart: " + str(too_long))
    logging.debug("reads which also pass both pair insert size filters: " + str(good_length) + " (" + str(int(100*float(good_length) / total_read_count)) + "%)")
    logging.debug("reads which pass minimum mapq threshold of " + str(min_mapq) + ": " + str(mapq_good) + " (" + str(int(100*float(mapq_good) / total_read_count)) +  "%)")
    logging.debug("(final) reads which also pass read pair PID >" + str(filter_cutoff) + "%: " + str(filter_cutoff_good) + " (" + str(int(100*float(filter_cutoff_good) / total_read_count)) + "%)")
    logging.debug("(final) expected coverage: " + str(float(filter_cutoff_good) * mapped_read_lengths / fasta_length))

    ## STEP 3: WRITE DATA IF NEEDED
    if write_data:
        f = open(write_data, 'w+')
        for read_pair in mapped_pairs:
            f.write(read_pair + "\t" + "\t" + str(pair_inserts[read_pair]) + "\t" + str(pair_mapqs[read_pair]) + "\t" + str(pair_mismatch[read_pair]) + "\n")
        f.close()
    ## STEP 4: WRITE NEW BAM IF NEEDED (TODO)

    samfile.close()
    if write_bam:
        samfile_out.close()
        logging.info("sorting new bam")
        pysam.sort("-o", bam.split("/")[-1].split(".")[0] + "_filtered_sort.bam", bam.split("/")[-1].split(".")[0] + "_filtered.bam")
        os.system('rm ' + bam.split("/")[-1].split(".")[0] + "_filtered.bam")

    return final_reads, Rdb
