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
        bam = args.bam
        vargs = vars(args)
        del vargs['bam']

        detailed_report = vargs.get('deatiled_read_report', False)
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
            Rdic, RR, dRR = load_paired_reads2(bam, scaffolds, **vargs)
        else:
            Rdic, RR = load_paired_reads2(bam, scaffolds, **vargs)
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

        RR_loc = os.path.join(out_folder, 'read_report.csv')
        write_read_report(RR, RR_loc, **kwargs)

        if dRR is not None:
            RR_loc = os.path.join(out_folder, 'deatiled_read_report.csv')
            dRR.to_csv(RR_loc, index=False, sep='\t')

def load_paired_reads2(bam, scaffolds, **kwargs):
    '''
    Load paired reads to be profiled

    You have this method do a lot of things because all of these things take lots of RAM, and you want them all to be cleared as soon as possible

    Return a dictionary of results. Some things that could be in it are:
        pair2infoF: A filtered dictionary of read name to number of mismatches
        RR: A summary read reaport
        RR_detailed: A detailed read report
    '''
    # Parse the kwargs
    detailed_report = kwargs.get('deatiled_read_report', False)

    # Get the pairs
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

    if detailed_report:
        dRR = make_detailed_read_report(scaff2pair2info, pairTOinfo=pair2info, version=2)
        return pair2infoF, RR, dRR

    else:
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

    elif pairing_filter == 'all_reads':
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                if p in pair2info:
                    pair2info[p] = _merge_info(i, pair2info[p])
                else:
                    pair2info[p] = i

    else:
        logging.error("Do not know paired read filter \"{0}\"; crashing now".format(pairing_filter))
        raise Exception

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

def make_detailed_read_report(scaff2pair2info, pairTOinfo=dict(), version=2):
    '''
    Make a detailed pandas dataframe from pair2info
    '''
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
    values['pairing_filter'] = kwargs.get('pairing_filter', 'paired_only')

    # Write header
    os.remove(location) if os.path.exists(location) else None
    f = open(location, 'a')
    f.write("# {0}\n".format(' '.join(["{0}:{1}".format(k, v) for k, v in values.items()])))

    # Write csv
    RR.to_csv(f, index=False, sep='\t')

    # Close
    f.close()

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
