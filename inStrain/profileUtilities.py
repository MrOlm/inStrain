# Import packages
import os
import sys
import time
import pysam
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

import inStrain.linkage
import inStrain.SNVprofile
import inStrain.readComparer

from ._version import __version__

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
        if self.number_splits == 1:
            return self.split_dict[0].merge_single_profile()
        else:
            Sprofile = scaffold_profile()
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

            Sprofile.make_cumulative_tables()

            return Sprofile

class SplitObject():
    '''
    Holds the profile of a split
    '''
    def __init__(self):
        pass

    def merge_single_profile(self):
        '''
        Convert self into scaffold_profile object
        '''
        Sprofile = scaffold_profile()
        Sprofile.scaffold = self.scaffold
        Sprofile.bam = self.bam
        Sprofile.length = self.length
        Sprofile.raw_snp_table = self.raw_snp_table
        Sprofile.raw_linkage_table = self.raw_linkage_table
        Sprofile.covT = self.covT
        Sprofile.clonT = self.clonT
        Sprofile.clonTR = self.clonTR
        Sprofile.min_freq = self.min_freq

        for att in ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts']:
            if hasattr(self, att):
                setattr(Sprofile, att, eval(att))

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
                            min_freq=self.min_freq)

class profile_command():
    '''
    This is a stupid object that just holds the arguments to profile a split
    '''
    def __init__(self):
        pass

def profile_contig_worker(available_index_queue, sprofile_cmd_queue, Sprofile_dict, log_list, Sprofiles):
    '''
    Based on https://github.com/merenlab/anvio/blob/18b3c7024e74e0ac7cb9021c99ad75d96e1a10fc/anvio/profiler.py
    '''

    # Process splits
    while not available_index_queue.empty():
        cmd = available_index_queue.get(True)
        try:
            Split = split_profile_wrapper(cmd)
            Sprofile_dict[cmd.scaffold + '.' + str(cmd.split_number)] = Split
            log_list.put(Split.log)
            #Sprofile_dict[cmd.scaffold].print_status()
        except:
            Sprofile_dict[cmd.scaffold + '.' + str(cmd.split_number)] = False
            log_list.put('FAILURE FOR SCAFFOLD {0} SPLIT {1}'.format(cmd.scaffold, cmd.split_number))

        # if (Sprofile_dict[cmd.scaffold].ready()):# & (cmd.scaffold not in available_index_queue):
        #     sprofile_cmd_queue.put(cmd.scaffold)
        # else:
        #     pass
            #Sprofile_dict[cmd.scaffold].print_status()

    # Merge splits into Sprofiles
    #while ((not sprofile_cmd_queue.empty()) | (len(Sprofile_dict) != 0)):
    while True:
        # If there's something in the queue, run it
        ScaffoldSplitObject = sprofile_cmd_queue.get(True)
        Sprofiles.append(ScaffoldSplitObject.merge())

        # See if there are any scaffolds to move into the queue
        # for scaff, obj in Sprofile_dict.items():
        #     if obj.ready():
        #         sprofile_cmd_queue.put(scaff)
        #         break


    print("Done!")
    return

        # if contig is not None:
        #     # We mark these for deletion the next time garbage is collected
        #     for split in contig.splits:
        #         del split.coverage
        #         del split.auxiliary
        #         del split
        #     del contig.splits[:]
        #     del contig.coverage
        #     del contig

def profile_bam(bam, Fdb, r2m, **kwargs):
    '''
    Profile the .bam file really  well.
    This is by far the meat of the program

    arguments:
        bam = location of .bam file
        Fdb = dictionary listing fasta locations to profile
        R2M = dictionary of read pair -> number of mm
    '''
    # get arguments for profiling the scaffolds
    profArgs = kwargs

    # get arguments for the wrapper
    p = int(kwargs.get('processes', 6))
    fdr = kwargs.get('fdr', 1e-6)

    global scaff2sequence
    scaff2sequence = profArgs.pop('s2s')

    # make the dictionary global
    global R2M
    R2M = r2m

    # make a global null model
    global null_model
    null_loc = os.path.dirname(__file__) + '/helper_files/NullModel.txt'
    null_model = generate_snp_model(null_loc, fdr=fdr)

    # make a list of the scaffolds to multiprocess
    scaffolds = list(Fdb['scaffold'].unique())

    # do the multiprocessing
    Sprofiles = []

    # Set up commands
    cmds, Sprofile_dict, s2splits = prepare_commands(Fdb, bam, profArgs)

    print(Fdb[Fdb['scaffold'] == 'N5_271_010G1_scaffold_0'])

    if p > 1:
    # p = 1
    # if p == 1:
        # Set up queues to be synced between the processes
        manager = multiprocessing.Manager()

        split_cmd_queue = manager.Queue() # Holds the commands to profile splits
        sprofile_cmd_queue = manager.Queue() # Holds the commands to merge splits
        log_list = manager.Queue() # Holds the resulting logs
        Sprofile_dict = manager.dict(Sprofile_dict) # Holds a synced directory of splits
        Sprofiles = manager.list(Sprofiles) # Holds the resulting Sprofiles

        # Fill in the queue with commands to profile splits
        for cmd in cmds:
            split_cmd_queue.put(cmd)

        # Start a number of processes to do this work
        processes = []
        for i in range(0, p):
            processes.append(multiprocessing.Process(target=profile_contig_worker, args=(split_cmd_queue, sprofile_cmd_queue, Sprofile_dict, log_list, Sprofiles)))

        for proc in processes:
            proc.start()

        # Set up progress bar
        pbar = tqdm(desc='Profiling splits: ', total=len(scaffolds))

        # Handle the queues
        while len(Sprofiles) < len(scaffolds):
            # See if there are any splits that are ready to be merged
            got = []
            for scaff, splits in s2splits.items():
                num_ready = sum([Sprofile_dict[scaff + '.' + str(i)] is not None for i in range(splits)])
                if num_ready == splits:
                    Sprofile = ScaffoldSplitObject(splits)
                    Sprofile.scaffold = scaff

                    for i in range(splits):
                        Sprofile = Sprofile.update_splits(i, Sprofile_dict.pop(scaff + '.' + str(i)))

                    sprofile_cmd_queue.put(Sprofile)
                    got.append(scaff)

            for g in got:
                pbar.update(1)
                s2splits.pop(g)


            # for scaff, obj in Sprofile_dict.items():
            #     if obj.ready():
            #         sprofile_cmd_queue.put(Sprofile_dict.pop(scaff))
            #     else:
            #         print("{0} not ready".format(obj.scaffold))


            time.sleep(1)
            print("split_cmd_queue {0}; sprofile_cmd_queue {1}; Sprofiles {2}; Sprofile dict {3}".format(split_cmd_queue.qsize(), sprofile_cmd_queue.qsize(), len(Sprofiles), len(Sprofile_dict)))
            # if Sprofile is False:
            #     failed_results += 1
            # else:
            #     Sprofiles.append(Sprofile)

            log_message = log_list.get()
            logging.debug(log_message)


        # Finish up multi-processing
        for proc in processes:
            proc.terminate()
        pbar.close()

    else:
        for cmd in tqdm(iterate_commands(Fdb, bam, profArgs), desc='Profiling scaffolds: ', total=total_cmds):
            results = scaffold_profile_wrapper(cmd)
            for s, log in results:
                logging.debug(log)
                Sprofiles.append(s)

    # Re-run failed scaffolds
    failed_scaffs = set(scaffolds) - set([s.scaffold for s in Sprofiles if s is not None])
    if len(failed_scaffs) > 0:
        logging.error("The following scaffolds failed scaffold profiling- I'll try again in single mode  {0}\n".format('\n'.join(failed_scaffs)))
        fdb = Fdb[Fdb['scaffold'].isin(failed_scaffs)]
        total_cmds = len([x for x in iterate_commands(fdb, bam, profArgs)])
        for cmd in tqdm(iterate_commands(fdb, bam, profArgs), desc='Profiling scaffolds: ', total=total_cmds):
            results = scaffold_profile_wrapper(cmd)
            try:
                for s, log in results:
                    logging.debug(log)
                    Sprofiles.append(s)
            except:
                logging.error("Double failure! Here are commands: {0}".format(cmd))

    # collate results
    logging.debug("Done processing scaffolds- making SNV profile")
    SNVprof = gen_snv_profile([s for s in Sprofiles if s is not None], **kwargs)

    # return results
    return SNVprof

def iterate_Fdbs(df, chunkSize=100):
    '''
    Break up Ndbs into chunks
    '''
    numberChunks = len(df) // chunkSize + 1
    for i in range(numberChunks):
        yield (df[i*chunkSize:(i+1)*chunkSize])

def generate_snp_model(model_file, fdr=1e-6):
    '''
    The model_file contains the columns "coverage", "probability of X coverage base by random chance"
    '''
    f = open(model_file)
    model = defaultdict(int)
    for line in f.readlines():
        if 'coverage' in line:
            continue
        counts = line.split()[1:]
        coverage = line.split()[0]
        i = 0
        for count in counts:
            if float(count) < fdr:
                model[int(coverage)] = i
                break
            i += 1

    model.default_factory = lambda:max(model.values())

    return model

def scaffold_profile_wrapper(cmds):
    '''
    Take a command and profile the scaffold
    '''
    try:
        results = []
        for cmd in cmds:
            results.append(_profile_scaffold(cmd.samfile, cmd.scaffold, cmd.locations, **cmd.arguments))
        return results
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(" ".join([cmd.scaffold for cmd in cmds]))))
        return pd.DataFrame({'Failed':[True]})

def scaffold_profile_wrapper2(cmd):
    '''
    Take a command and profile the scaffold
    '''
    try:
        return _profile_scaffold(cmd.samfile, cmd.scaffold, cmd.start, cmd.end, **cmd.arguments)
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(cmd.scaffold)))
        return pd.DataFrame({'Failed':[True]})

def split_profile_wrapper(cmd):
    try:
        return _profile_split(cmd.samfile, cmd.scaffold, cmd.start, cmd.end, cmd.split_number, **cmd.arguments)
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("split exception- {0} pt {1}".format(str(cmd.scaffold), str(cmd.split_number)))
        return False

def prepare_commands(Fdb, bam, args):
    '''
    Make and iterate profiling commands
    Doing it in this way makes it use way less RAM
    '''
    cmds = []
    Sprofiles = {}
    Sdict = {}
    s2splits = {}

    for scaff, db in Fdb.groupby('scaffold'):
        Sprofile = ScaffoldSplitObject(len(db))
        Sprofile.scaffold = scaff
        s2splits[scaff] = len(db)

        for i, row in db.iterrows():

            # make this command
            cmd = profile_command()
            cmd.scaffold = scaff
            cmd.samfile = bam
            cmd.arguments = args
            cmd.start = row['start']
            cmd.end = row['end']
            cmd.split_number = int(row['split_number'])
            cmds.append(cmd)

            Sdict[scaff + '.' + str(row['split_number'])] = None

        Sprofiles[scaff] = Sprofile


    return cmds, Sdict, s2splits

def calc_estimated_runtime(pairs):
    SLOPE_CONSTANT = 0.0061401594694834305
    return pairs * SLOPE_CONSTANT

def _profile_split(bam, scaffold, start, end, split_number, **kwargs):
    '''
    Run the meat of the program to profile a split and return a split object

    Start and end are both inclusive and 0-based
    '''
    # Log
    log_message = _get_log_message('profile_start', scaffold, split_number=split_number)

    # Get kwargs
    min_cov = int(kwargs.get('min_cov', 5))
    min_covR = int(kwargs.get('rarefied_coverage', 5))
    min_freq = float(kwargs.get('min_freq', .05))
    min_snp = int(kwargs.get('min_snp', 10))
    store_everything = kwargs.get('store_everything', False)

    # Get sequence from global
    seq = scaff2sequence[scaffold][start:end+1] # The plus 1 makes is so the end is inclusive
    print("{0} to {1} is {2} long".format(start, end, len(seq)))

    # Set up the .bam file
    samfile = pysam.AlignmentFile(bam)
    try:
        iter = samfile.pileup(scaffold, truncate=True, max_depth=100000,
                                stepper='nofilter', compute_baq=True,
                                ignore_orphans=True, ignore_overlaps=True,
                                min_base_quality=30, start=start, stop=end+1)
    except ValueError:
        logging.error("scaffold {0} is not in the .bam file {1}!".format(scaffold, bam))
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

    # Do per-site processing
    _process_bam_sites(scaffold, seq, iter, covT, clonT, clonTR, p2c, read_to_snvs, snv2mm2counts, Stable, pileup_counts, mLen, start=start, **kwargs)

    # Shrink these dictionaries into index pandas Series
    covT = shrink_basewise(covT, 'coverage', start=start, len=mLen)
    clonT = shrink_basewise(clonT, 'clonality', start=start, len=mLen)
    clonTR = shrink_basewise(clonTR, 'clonality', start=start, len=mLen)

    # Make the SNP table
    SNPTable = _make_snp_table2(Stable, scaffold, p2c)
    if len(SNPTable) > 0:
        SNPTable['position'] = SNPTable['position'] + start

    # Do per-scaffold aggregating
    #CoverageTable = make_coverage_table(covT, clonT, clonTR, mLen, scaffold, SNPTable, min_freq=min_freq)

    # Make linkage network and calculate linkage
    mm_to_position_graph = inStrain.linkage.calc_mm_SNV_linkage_network(read_to_snvs, scaff=scaffold)
    LDdb = inStrain.linkage.calculate_ld(mm_to_position_graph, min_snp, snv2mm2counts=snv2mm2counts, scaffold=scaffold)
    if len(LDdb) > 0:
        for p in ['position_A', 'position_B']:
            LDdb[p] = LDdb[p] + start

    # Make a Scaffold profile to return
    Sprofile = SplitObject()
    Sprofile.scaffold = scaffold
    Sprofile.split_number = split_number
    Sprofile.bam = bam
    Sprofile.length = mLen
    #Sprofile.cumulative_scaffold_table = CoverageTable
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

    # Make cummulative tables
    #Sprofile.make_cumulative_tables()

    # Return
    log_message +=  _get_log_message('profile_end', scaffold, split_number=split_number)
    Sprofile.log = log_message

    return Sprofile

def _profile_scaffold(bam, scaffold, start, end, **kwargs):
    '''
    Profile a scaffold

    Arguments:
        bam = name of .bam file
        scaffold = name of scaffold
        seq = sequence of .fasta file thats mapped to
    '''
    # Log
    log_message = _get_log_message('profile_start', scaffold)

    # Get kwargs
    min_cov = int(kwargs.get('min_cov', 5))
    min_covR = int(kwargs.get('rarefied_coverage', 5))
    min_freq = float(kwargs.get('min_freq', .05))
    min_snp = int(kwargs.get('min_snp', 10))
    store_everything = kwargs.get('store_everything', False)

    # Get sequence from global
    seq = scaff2sequence[scaffold]

    # Set up the .bam file
    samfile = pysam.AlignmentFile(bam)
    try:
        iter = samfile.pileup(scaffold, truncate = True, max_depth=100000,
                                stepper = 'nofilter', compute_baq= True,
                                ignore_orphans = True, ignore_overlaps = True,
                                min_base_quality = 30, start=start, stop=end)
    except ValueError:
        logging.error("scaffold {0} is not in the .bam file {1}!".format(scaffold, bam))
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

    # Do per-site processing
    _process_bam_sites(scaffold, seq, iter, covT, clonT, clonTR, p2c, read_to_snvs, snv2mm2counts, Stable, pileup_counts, mLen, **kwargs)

    # Shrink these dictionaries into index pandas Series
    covT = shrink_basewise(covT, 'coverage')
    clonT = shrink_basewise(clonT, 'clonality')
    clonTR = shrink_basewise(clonTR, 'clonality')

    # Make the SNP table
    SNPTable = _make_snp_table2(Stable, scaffold, p2c)

    # Do per-scaffold aggregating
    CoverageTable = make_coverage_table(covT, clonT, clonTR, mLen, scaffold, SNPTable, min_freq=min_freq)

    # Make linkage network and calculate linkage
    mm_to_position_graph = inStrain.linkage.calc_mm_SNV_linkage_network(read_to_snvs, scaff=scaffold)
    LDdb = inStrain.linkage.calculate_ld(mm_to_position_graph, min_snp, snv2mm2counts=snv2mm2counts, scaffold=scaffold)

    # Make a Scaffold profile to return
    Sprofile = scaffold_profile()
    Sprofile.scaffold = scaffold
    Sprofile.bam = bam
    Sprofile.length = mLen
    Sprofile.cumulative_scaffold_table = CoverageTable
    Sprofile.raw_snp_table = SNPTable
    Sprofile.raw_linkage_table = LDdb
    Sprofile.covT = covT
    Sprofile.clonT = clonT
    Sprofile.clonTR = clonTR

    # store extra things if required
    if store_everything:
        for att in ['read_to_snvs', 'mm_to_position_graph', 'pileup_counts']:
            setattr(Sprofile, att, eval(att))

    # Make cummulative tables
    Sprofile.make_cumulative_tables()

    # Return
    log_message +=  _get_log_message('profile_end', scaffold)
    return Sprofile, log_message

def _process_bam_sites(scaffold, seq, iter, covT, clonT, clonTR, p2c, read_to_snvs, snv2mm2counts, Stable, pileup_counts, mLen, start=0, **kwargs):
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
        MMcounts = _get_base_counts_mm(pileupcolumn)
        _update_covT(covT, MMcounts, RelPosition, mLen)

        # Call SNPs
        snp, bases, total_counts = _update_snp_table_T(Stable, clonT, clonTR,\
                    MMcounts, p2c,\
                    RelPosition, scaffold, mLen, seq[RelPosition], min_cov=min_cov, min_covR=min_covR, min_freq=min_freq)

        # add the counts for this position to the numpy array alexcc
        if store_everything:
            pileup_counts[RelPosition] = total_counts

        # get linked reads
        if snp:
            _update_linked_reads(read_to_snvs, pileupcolumn, MMcounts, RelPosition, bases, scaffold=scaffold)
            snv2mm2counts[RelPosition] = MMcounts

def _get_log_message(kind, scaffold, split_number=False):
    if kind == 'profile_start':
        pid = os.getpid()
        process  = psutil.Process(os.getpid())
        bytes_used = process.memory_info().rss
        total_available_bytes = psutil.virtual_memory()
        log_message = "\n{4} PID {0} start at {5} with {1} RAM. System has {2} of {3} available".format(
                pid, bytes_used, total_available_bytes[1], total_available_bytes[0],
                scaffold, time.time())
        return log_message

    elif kind == 'profile_end':
        pid = os.getpid()
        process  = psutil.Process(os.getpid())
        bytes_used = process.memory_info().rss
        total_available_bytes = psutil.virtual_memory()
        return "\n{4} PID {0} end at {5} with {1} RAM. System has {2} of {3} available".format(
                pid, bytes_used, total_available_bytes[1], total_available_bytes[0],
                scaffold, time.time())

def _dlist():
    return defaultdict(list)

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
def _get_base_counts_mm(pileupcolumn):
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

def _4zeros():
    return np.zeros(4, dtype=int)

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

def merge_basewise(mm2array_list):
    NAME2TYPE = {'coverage':'int32', 'clonality':'float32', 'snpCounted':'bool'}

    mm2bases = {}

    mms = set()
    for mm2bases in mm2array_list:
        mms = mms.union(set(mm2bases.keys()))

    for mm in mms:
        mm2bases[mm] = pd.concat([mm2array[mm] for mm2array in mm2array_list if mm in mm2array], verify_integrity=True)
    return mm2bases


    # mm2bases = {}
    # for mm in list(mm2array.keys()):
    #     if start != 0:
    #         db = pd.Series(mm2array[mm], index=pd.RangeIndex(start=start, stop=start+len), dtype=NAME2TYPE[name]).dropna()
    #     else:
    #         db = pd.Series(mm2array[mm], dtype=NAME2TYPE[name]).dropna()
    #     db = db[db > 0]
    #     mm2bases[mm] = db
    #     del mm2array[mm]
    #
    # return mm2bases

def _update_covT(covT, MMcounts, position, mLen):
    '''
    Update covT at this position
    '''
    for mm, count in MMcounts.items():
        if mm not in covT:
            covT[mm] = np.zeros(mLen, dtype=int)
        covT[mm][position] = sum(count)

def _update_snp_table_T(Stable, clonT, clonTR, MMcounts, p2c,\
        pos, scaff, mLen, refBase, min_cov=5, min_covR=50, min_freq=.05):
    '''
    Add information to SNP table. Update basesCounted and snpsCounted

    Also update clonality

    cryptic means that it's a SNP at lower mm levels, and then not again at higher levels

    6.10.18 - v0.4.1 - Notes

    What makes it into the Stable?
    *  There must be two bases that satisfy the following conditions:
        *  There are more than min_cov number of bases counted
        *  There are more than min_freq percent of reads at the variant base
        *  The number of variant bases is more than the null model for that coverage
    *  If those are satisfied twice, the base with the highest count is stored as the "varBase",
       and the the bases with the second highest is the "refBase"

    What is anySNP?
    *  For at least one mm position, there was a SNP called

    What is bases?
    *  A set of bases that are the varBase or refBase at least once

    What are the things that are being missed?
    *  The case where theres only one base that meets the chriteria above, but its
       different from the reference base
    *  The case where >2 variant bases are present

    What is SNP table used for?
    *  Looking for specific variant positions; performing dn/ds analysis
    *  Establishing linkage patterns

    What to do about SNPs that are fixed opposite of reference?
    *  Create a "allele_count" column (formerly morphia); 2 means there are 2 bases; 1 means theres 1
        (and its opposite of reference); 3/4 means there multiple, though only
        two will be logged
        *  Right now only 2+ are being shown; so it will be easy to filter others out
    *  Throw in a check for if there's only one bases that fits that chriteria,
       is it not the reference base?
    *  Add "refbase" back into the table to make it easier to compare SNP tables

    '''
    anySNP = False
    bases = set()
    ret_counts = np.zeros(4, dtype=int)
    for mm in sorted(list(MMcounts.keys())):
        counts = _mm_counts_to_counts(MMcounts, mm)
        snp, morphia = call_snv_site(counts, refBase, min_cov=min_cov, min_freq=min_freq) # Call SNP

        # Update clonality
        if mm not in clonT:
            clonT[mm] = np.full(mLen, fill_value=np.nan, dtype="float32")
        if sum(counts) >= min_cov:
            clonT[mm][pos] = calculate_clonality(counts)

        # Update rarefied clonality
        if mm not in clonTR:
            clonTR[mm] = np.full(mLen, fill_value=np.nan, dtype="float32")
        if sum(counts) >= min_covR:
            clonTR[mm][pos] = calculate_rarefied_clonality(counts, rarefied_coverage=min_covR)

        if not snp: # means base was not counted
            continue

        elif snp != -1: # means this is a SNP

           # calculate varBase
           counts_temp = list(counts)
           counts_temp[P2C[snp]] = 0
           varbase = C2P[list(counts_temp).index(sorted(counts_temp)[-1])] # this fixes the varbase = refbase error when there's a tie - alexcc 5/8/2019

           Stable['scaffold'].append(scaff)
           Stable['position'].append(pos)
           Stable['refBase'].append(refBase)
           for b, c in zip(['A', 'C', 'T', 'G'], counts):
               Stable[b].append(c)
           Stable['conBase'].append(snp)
           Stable['varBase'].append(varbase)
           Stable['mm'].append(mm)
           Stable['allele_count'].append(morphia)

           if morphia >= 2:
               anySNP = True
               bases.add(snp)
               bases.add(varbase)
               ret_counts = counts

           elif (morphia == 1) & (anySNP == True):
               p2c[pos] = True

           # if mm not in snpsCounted:
           #     snpsCounted[mm] = np.zeros(mLen, dtype=bool)
           # snpsCounted[mm][pos] = True

        # if it's now not a SNP, but it was in the past, mark it cryptic
        elif (snp == -1) & (anySNP == True):
            p2c[pos] = True

        # if mm not in basesCounted:
        #     basesCounted[mm] = np.zeros(mLen, dtype=bool)
        # basesCounted[mm][pos] = True # count everything that's not skipped

    return anySNP, bases, ret_counts

def _mm_counts_to_counts(MMcounts, maxMM=100):
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

def _mm_counts_to_counts_shrunk(MMcounts, maxMM=100, fill_zeros=False):
    '''
    Take mm counts and return just counts
    THIS IS WITH THE SHRUNK DATAFRAME

    If fill_zeros is an int, at the end make sure the array is of that lengths
    '''
    counts = pd.Series()
    for mm, count in [(mm, count) for mm, count in MMcounts.items() if mm <= maxMM]:
        counts = counts.add(count, fill_value=0)

    if fill_zeros:
        counts = counts.append(pd.Series(np.zeros(fill_zeros - len(counts))))

    return counts


P2C = {'A':0, 'C':1, 'T':2, 'G':3}
C2P = {0:'A', 1:'C', 2:'T', 3:'G'}
def call_snv_site(counts, refBase, min_cov=5, min_freq=0.05, model=None):
    '''
    Determines whether a site has a variant based on its nucleotide count frequencies.

    Returns one of the following and the "morphia" (number of bases present in the reads):
        Base if SNP
        -1 if not SNP
        None if unCounted

    A base is considered present in the reads if:
    *  There are more than min_cov number of bases counted
    *  There are more than min_freq percent of reads at the variant base
    *  The number of variant bases is more than the null model for that coverage

    The morphia is the number of bases present in the reads. Things are returned as SNPs if:
    *  Morphia >= 2
    *  Morphia == 1, and the base present is not the genome reference base
    *  Morphia == 0 (this is a special case meaning no bases are really present)
    '''
    # Get the null model
    if model: #alexcc - so you can call this function from outside of the file
        model_to_use = model
    else:
        model_to_use = null_model

    # Make sure you have the coverage
    total = sum(counts)
    if total < min_cov:
        return None, 0

    # Count how many bases are there
    i = 0
    for c in counts:
        if c >= model_to_use[total] and float(c) / total >= min_freq:
            i += 1

    # If there are 2, say that
    if i > 1:
        return C2P[np.argmax(counts)], i

    # If theres only 1, see if its the reference base
    elif (i == 1) & (C2P[np.argmax(counts)] != refBase):
        return C2P[np.argmax(counts)], i

    # That means this is a super polymorphic position with no dominant bases
    elif i == 0:
        return C2P[np.argmax(counts)], i

    # This means its just not a SNP; one dominant reference base
    else:
        return -1, i


def calculate_clonality(counts):
    '''
    Calculates the probability that two reads have the same allele at a position (per site nucleotide diversity)
    '''
    s = sum(counts)
    prob = (float(counts[0]) / s) * (float(counts[0]) / s) + (float(counts[1]) / s) * (float(counts[1]) / s) + (float(counts[2]) / s) * (float(counts[2]) / s) + (float(counts[3]) / s) * (float(counts[3]) / s)
    return prob

def calculate_rarefied_clonality(counts, rarefied_coverage=50):
    '''
    Calculates the probability that two reads have the same allele at a position (per site nucleotide diversity)
    '''
    # Figure out the probability distribution
    s = sum(counts)
    p = [counts[i]/s for i in [0,1,2,3]]

    # Make rarefied counts
    rcounts_list = np.random.choice([0,1,2,3], rarefied_coverage, p=p)
    unique, ncounts = np.unique(rcounts_list, return_counts=True)
    item2counts = dict(zip(unique, ncounts))
    rcounts = [item2counts[i] if i in item2counts else 0 for i in [0,1,2,3]]

    return calculate_clonality(rcounts)

def get_lowest_mm(clonT, mm):
    mms = [int(m) for m in list(clonT.keys()) if int(m) <= int(mm)]
    if len(mms) == 0:
        return None
    else:
        return max(mms)

def _get_basewise_clons2(clonT, MM, fill_zeros=False):
    p2c = {}
    mms = sorted([int(mm) for mm in list(clonT.keys()) if int(mm) <= int(MM)])
    for mm in mms:
        p2c.update(clonT[mm].to_dict())

    counts = list(p2c.values())

    if fill_zeros:
        counts = counts.append(pd.Series(np.zeros(fill_zeros - len(counts))))

    return counts

def _calc_snps(Odb, mm, min_freq=0.05):
    '''
    Calculate the number of reference SNPs, bi-allelic SNPs, multi-allelic SNPs (>2), total SNPs, consensus_SNPs, and poplation_SNPs
    '''
    if len(Odb) == 0:
        return [0, 0, 0, 0, 0, 0]

    db = Odb[Odb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['position'], keep='last')

    ref_snps = len(db[(db['allele_count'] == 1)])
    bi_snps = len(db[(db['allele_count'] == 2)])
    multi_snps = len(db[(db['allele_count'] > 2)])

    con_snps = len(db[(db['conBase'] != db['refBase']) & (db['allele_count'] > 0)])

    # One pool of population SNPs are those of morphia 1 where con != ref
    p1 = len(db[(db['refBase'] != db['conBase']) & (db['allele_count'] == 1)])

    # Another pool is when its biallelic but neither are the reference
    p2 = len(db[(db['refBase'] != db['conBase']) & (db['allele_count'] == 2) & (db['refBase'] != db['varBase'])])

    # Finally, and the hardest to detect, are morphia of three where
    p3 = 0
    pdb = db[(db['refBase'] != db['conBase']) & (db['allele_count'] == 3) & (db['refBase'] != db['varBase'])]
    for i, row in pdb.iterrows():
        # Below is trying to get the count of Ns; but there is no count of Ns
        if row['refBase'] not in ['A', 'C', 'T', 'G']:
            continue

        if not inStrain.readComparer.is_present(int(row[row['refBase']]), int(row['baseCoverage']), null_model, float(min_freq)):
            p3 += 1

    pop_snps = p1 + p2 + p3

    return [ref_snps, bi_snps, multi_snps, len(db), con_snps, pop_snps]

def make_coverage_table(covT, clonT, clonTR, lengt, scaff, SNPTable, min_freq=0.05, debug=False):
    '''
    Add information to the table
    Args:
        covT: list of coverage values
        lengt: length of scaffold
        scaff: name of scaffold
        Wdb: windows to process over (NOT YET SUPPORTED)
    '''
    table = defaultdict(list)
    #CLdb = _clonT_to_table(clonT)
    for mm in sorted(list(covT.keys())):

        covs = _mm_counts_to_counts_shrunk(covT, mm, fill_zeros=lengt)

        if len(covs) == 0:
            covs = pd.Series([0]*lengt)

        nonzeros = np.count_nonzero(covs)
        zeros = lengt - nonzeros

        # Get clonalities
        clons = _get_basewise_clons2(clonT, mm)
        Rclons = _get_basewise_clons2(clonTR, mm)

        counted_bases = len(clons)
        rarefied_bases = len(Rclons)
        ref_snps, bi_snps, multi_snps, counted_snps, con_snps, pop_snps = _calc_snps(SNPTable, mm, min_freq=min_freq)
        #counted_snps = _calc_counted_bases(snpsCounted, mm)

        assert len(covs) == lengt, [covs, lengt, mm]

        # fill in all coverage information
        table['scaffold'].append(scaff)
        table['length'].append(lengt)
        table['breadth'].append(nonzeros/lengt)
        table['coverage'].append(np.mean(covs))
        table['median_cov'].append(int(np.median(covs)))
        table['std_cov'].append(np.std(covs))
        table['bases_w_0_coverage'].append(zeros)

        if len(clons) > 0:
            mean_c = np.mean(clons)
            median_c = np.median(clons)
            table['mean_clonality'].append(mean_c)
            table['median_clonality'].append(median_c)
            table['mean_microdiversity'].append(1-mean_c)
            table['median_microdiversity'].append(1-median_c)
        else:
            table['mean_clonality'].append(np.nan)
            table['median_clonality'].append(np.nan)
            table['mean_microdiversity'].append(np.nan)
            table['median_microdiversity'].append(np.nan)

        if len(Rclons) > 0:
            mean_c = np.mean(Rclons)
            median_c = np.median(Rclons)
            table['rarefied_mean_microdiversity'].append(1-mean_c)
            table['rarefied_median_microdiversity'].append(1-median_c)
        else:
            table['rarefied_mean_microdiversity'].append(np.nan)
            table['rarefied_median_microdiversity'].append(np.nan)

        table['unmaskedBreadth'].append(len(clons) / lengt)
        table['rarefied_breadth'].append(rarefied_bases / lengt)
        table['expected_breadth'].append(estimate_breadth(table['coverage'][-1]))

        table['SNPs'].append(counted_snps)

        table['Reference_SNPs'].append(ref_snps)
        table['BiAllelic_SNPs'].append(bi_snps)
        table['MultiAllelic_SNPs'].append(multi_snps)

        table['consensus_SNPs'].append(con_snps)
        table['population_SNPs'].append(pop_snps)

        if counted_bases == 0:
            table['conANI'].append(0)
            table['popANI'].append(0)
        else:
            table['conANI'].append((counted_bases - con_snps)/ counted_bases)
            table['popANI'].append((counted_bases - pop_snps)/ counted_bases)

        table['mm'].append(mm)

    # if debug == True:
    #     covs = _mm_counts_to_counts(covT, max(list(covT.keys())))
    #     zero_pos = [i+1 for i,x in enumerate(covs) if x==0]
    #     print(zero_pos)
    #     print(len(zero_pos))

    return pd.DataFrame(table)

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
    # for mm, x in clonT.items():
    #     items += [(c, i) for i, c in enumerate(x) if c == c]
    #     mms += [mm]*(len(items) - len(mms))
    # db = pd.DataFrame({'clonality':[i[0] for i in items], 'mm':mms,
    #                     'position':[i[1] for i in items]})
    # db['position'] = db['position'].astype(int)
    #return db

# def _clonT_to_table(clonT):
#     '''
#     Much more RAM efficient way of doing this than before!
#     '''
#     items = []
#     mms = []
#     for mm, x in clonT.items():
#         items += [(c, i) for i, c in enumerate(x) if c == c]
#         mms += [mm]*(len(items) - len(mms))
#     db = pd.DataFrame({'clonality':[i[0] for i in items], 'mm':mms,
#                         'position':[i[1] for i in items]})
#     db['position'] = db['position'].astype(int)
#     return db

def _update_linked_reads(read_to_snvs, pileupcolumn, MMcounts, position, bases,
                        min_freq=0.05, scaffold=False):
    '''
    Find and update linked reads

    try:
        table[R2M[pileupread.alignment.query_name]]\
        [P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
    '''

    position = str(position)

    # Get reads to snvs
    for pileupread in pileupcolumn.pileups:
        read_name = pileupread.alignment.query_name
        if not pileupread.is_del and not pileupread.is_refskip and read_name in R2M:
            if type(R2M) == type({}):
                mm = R2M[read_name]
            else:
                mm = 0
            #local_min_snp = mm2min[mm]

            try:
                val = pileupread.alignment.query_sequence[pileupread.query_position]
                #val_counts = mm2counts[mm][P2C[val]]

                # search here for debug
                # if (scaffold == 'N5_271_010G1_scaffold_1') & ((position == '2130') | (position == '2139')):
                #     print('At position {0}, mm {1}'.format(position, mm))
                #     print("{0} - {1}".format(read_name, val))#, val_counts))

                # if value is not the consensus value
                if val in bases:
                #if val_counts >= local_min_snp and float(val_counts) / sum(mm2counts[mm]) >= min_freq:
                    # this read is of a variant position
                    read_to_snvs[mm][read_name].append(position + ":" + val)
            except KeyError: # This would be like an N or something not A/C/T/G
                pass

def _calc_counted_bases(basesCounted, maxMM):
    '''
    Return the number of bases counted at a particular mm level
    Returns the number of bases at this level and all levels below it
    '''
    counts = None
    for mm, count in [(mm, count) for mm, count in basesCounted.items() if mm <= maxMM]:
        if counts is None:
            counts = count
        else:
            counts = counts.add(count, fill_value=False)

    if counts is None:
        return 0

    else:
        return counts.astype('bool').sum()

def _make_snp_table(Stable):
    if Stable is not False:
        try:
            Sdb = pd.DataFrame(Stable)
            Sdb['scaffold'] = Sdb['scaffold'].astype('category')
            Sdb['conBase'] = Sdb['conBase'].astype('category')
        except KeyError:
            #logging.info("No SNPs detected!")
            Sdb = pd.DataFrame()
    else:
        Sdb = pd.DataFrame()

    return Sdb

def _make_snp_table2(Stable, scaffold, p2c):
    SNPTable = pd.DataFrame(Stable)
    if len(SNPTable) > 0:
        SNPTable['scaffold'] = scaffold

        # Add cryptic SNPs
        SNPTable['cryptic'] = SNPTable['position'].map(p2c)
        SNPTable['cryptic'] = SNPTable['cryptic'].fillna(False)

        # Calc base coverage
        SNPTable['baseCoverage'] = [sum([a,c,t,g]) for a,c,t,g in zip(SNPTable['A'],SNPTable['C'],SNPTable['T'],SNPTable['G'])]

    del Stable
    return SNPTable

def _parse_Sdb(sdb):
    '''
    Add some information to sdb
    '''
    if len(sdb) == 0:
        return sdb

    sdb['varFreq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['varBase'], sdb['baseCoverage'])]
    sdb['conFreq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['conBase'], sdb['baseCoverage'])]
    sdb['refFreq'] = [[a,c,t,g][['A','C','T','G'].index(v)]/s if v in ['A','C','T','G'] else np.nan for a,c,t,g,v,s in zip(\
                        sdb['A'], sdb['C'], sdb['T'], sdb['G'], sdb['refBase'], sdb['baseCoverage'])]

    return sdb



def gen_snv_profile(Sprofiles, **kwargs):
    '''
    Take a bunch of scaffold profiles and make an SNVprofile
    '''
    location = kwargs.get('output', False)
    store_everything = kwargs.get('store_everything', False)

    # Merge things
    scaffold_list = []
    bam_list = []
    scaffold2length = {}

    #raw_cov_dbs = []
    #raw_ani_dbs = []
    raw_snp_dbs = []
    raw_link_dbs = []

    cumu_scaff_dbs = []
    cumu_snv_dbs = []

    scaffold_2_mm_2_read_2_snvs = {}

    for Sprof in Sprofiles:

        scaffold_list.append(Sprof.scaffold)
        bam_list.append(Sprof.bam)
        scaffold2length[Sprof.scaffold] = Sprof.length

        #raw_cov_dbs.append(Sprof.raw_coverage_table)
        #raw_ani_dbs.append(Sprof.raw_ANI_table)
        raw_snp_dbs.append(Sprof.raw_snp_table)
        raw_link_dbs.append(Sprof.raw_linkage_table)

        cumu_scaff_dbs.append(Sprof.cumulative_scaffold_table)
        cumu_snv_dbs.append(Sprof.cumulative_snv_table)

        if hasattr(Sprof, 'mm_reads_to_snvs'):
            scaffold_2_mm_2_read_2_snvs[Sprof.scaffold] = Sprof.mm_reads_to_snvs

    # Make some dataframes
    raw_snp_table = pd.concat(raw_snp_dbs).reset_index(drop=True)
    if len(raw_snp_table) > 0:
        for col in ['A', 'C', 'G', 'T', 'mm', 'position']:
            raw_snp_table[col] = raw_snp_table[col].astype(int)
            raw_snp_table['scaffold'] = raw_snp_table['scaffold'].astype('category')
            raw_snp_table['conBase'] = raw_snp_table['conBase'].astype('category')

    # convert to numpy array
    #counts_table = np.array(counts_table)

    # Make the profile
    Sprofile = inStrain.SNVprofile.SNVprofile(location)

    # Add some things
    Sprofile.store('bam_loc', bam_list[0], 'value', 'Location of .bam file')
    Sprofile.store('scaffold_list', scaffold_list, 'list', '1d list of scaffolds, in same order as counts_table')
    Sprofile.store('scaffold2length', scaffold2length, 'dictionary', 'Dictionary of scaffold 2 length')
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

    # On the fly testing
    # if True:
    #     assert scaffold_2_mm_2_read_2_snvs == Sprofile.get('scaffold_2_mm_2_read_2_snvs')

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
