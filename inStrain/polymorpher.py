#!/usr/bin/env python

"""
Methods related to SNV pooling and extracting SNVs from .bam files
"""

import os
import sys
import json
import h5py
import time
import psutil
import pickle
import logging
import argparse
import traceback
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
import multiprocessing
import concurrent.futures
from concurrent import futures
from collections import defaultdict

from Bio import SeqIO

if __name__ != '__main__':
    from ._version import __version__

import inStrain.profile.fasta
import inStrain.profile.snv_utilities
import inStrain.profile.profile_utilities

import inStrain.SNVprofile
import inStrain.controller
import inStrain.logUtils


class PoolController(object):
    """
    Main controller of the polymorpher pooling commands
    """

    def __init__(self, SNPtables, names, name2Rdic, name2bam_loc, scaffold, **kwargs):
        """
        Main constructor; to be called from within compare

        *Right now this is set up for a single scaffold*
        """
        self.scaffold = scaffold

        # Edit the SNPtables to remove mm
        self.SNPtables = [x.sort_values('mm')\
            .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
            .sort_index().drop(columns=['mm']) if 'mm' in list(x.columns) else x for x in SNPtables]

        self.names = names
        self.name2Rdic = name2Rdic
        self.name2bam_loc = name2bam_loc

        self.verify_input()

    def verify_input(self):
        """
        Check that input is OK
        """
        # 1) Make sure you don't have SNPtables on mm level
        for t in self.SNPtables:
            assert(len(t) == len(t.drop_duplicates(subset=['scaffold', 'position'])))

    def main(self):
        """
        The main method
        """
        # Extract a list of SNV sites
        self.name2locs, self.all_SNVs = extract_SNV_positions(self.SNPtables, self.names)

        # Pull needed SNV sites from .bam files
        self.pull_SNVS_from_bams()

        # Run processing to create the pooled table
        self.make_pooled_SNV_table()

        # Run processing to make the SNV info table
        self.make_pooled_SNV_info_table()


    def pull_SNVS_from_bams(self):
        """
        Using attributes of this object, pull the needed SNVs from .bam files
        """
        name2position2counts = {}
        for name in self.names:
            name2position2counts[name] = extract_SNVS_from_bam(self.name2bam_loc[name], self.name2Rdic[name], self.name2locs[name], self.scaffold)
        self.name2position2counts = name2position2counts

    def make_pooled_SNV_table(self):
        """
        Using attributes of this object, create an SNV table
        """
        self.DDST = genreate_pooled_SNV_table(self.SNPtables, self.names, self.name2position2counts, self.all_SNVs)

    def make_pooled_SNV_info_table(self):
        """
        Using the pooled SNV table, create a kind of "summary" object
        """
        self.PST = generate_pooled_SNV_summary_table(self.DDST, self.SNPtables, self.names)


def extract_SNVS_from_bam(bam_loc, R2M, positions, scaffold, **kwargs):
    """
    From a .bam file and a list of SNV locations, extract and return counts

    MY GOD IS PYSAM A PAIN - the reason I have the "start" and "stop" so weird on that iterator is because without it, this fails
    """
    logging.debug(f"Extracting {len(positions)} SNVs from {scaffold} {bam_loc}")

    # Initilize .bam
    bam_init = pysam.AlignmentFile(bam_loc)

    # Initialize the object to return
    position2counts = {}

    # Iterate positions
    for pos in positions:
        pileupcolumn = None
        p = pos

        # Get the pileup for this position
        try:
            iter = bam_init.pileup(scaffold, truncate=True, max_depth=100000,
                                  stepper='nofilter', compute_baq=True,
                                  ignore_orphans=True, ignore_overlaps=True,
                                  min_base_quality=30, start=p - 1, stop=p + 1)
            for x in iter:
                if x.pos == pos:
                    pileupcolumn = x
                    break
            assert p == pileupcolumn.pos, pileupcolumn

        except ValueError:
            logging.error("scaffold {0} position {2} is not in the .bam file {1}!".format(scaffold, bam_loc, p))
            continue

        # Process the pileup column
        MMcounts = inStrain.profile.profile_utilities.get_base_counts_mm(pileupcolumn, R2M)


        counts = inStrain.profile.profile_utilities.mm_counts_to_counts(MMcounts)
        position2counts[pos] = counts

    return position2counts

def extract_SNV_positions(SNPtables, names):
    """
    From a list of tables and names, figure out how many SNVs need to be extracted from .bam files
    """
    name2locs = {}
    ALL_SNPS = set()

    # Loop once to get all locations
    for sdb in SNPtables:
        if len(sdb) == 0:
            continue

        assert len(sdb['scaffold'].unique()) == 1
        ALL_SNPS = ALL_SNPS.union(set(sdb['position']))

    # Loop once to get all needed locations
    for sdb, name in zip(SNPtables, names):
        if len(sdb) > 0:
            name2locs[name] = ALL_SNPS - set(sdb['position'])
        else:
            name2locs[name] = ALL_SNPS

    return name2locs, ALL_SNPS

def genreate_pooled_SNV_table(SNPtables, names, name2position2counts, all_SNVs):
    """
    Generate the "DDST" (data-dense SNV table) and "Info table" that can be used for further analysis
    """

    dbs = []
    for name, ori_table in zip(names, SNPtables):
        # Make a table for this sample from pulled SNVs
        position2counts = name2position2counts[name]
        db = pd.DataFrame.from_dict(position2counts, orient='index', columns=['A', 'C', 'T', 'G'])

        # Merge with existing table
        if len(ori_table) > 0:
            sdb = pd.concat([db, ori_table[['position', 'A', 'C', 'T', 'G']].set_index('position')]).sort_index()
        else:
            sdb = db.sort_index()

        # Verify
        assert set(sdb.index) == set(all_SNVs)
        dbs.append(sdb)

    DDST = pd.concat(dbs, keys=names)
    return DDST

def generate_pooled_SNV_summary_table(DDST, SNPtables, names):
    """
    Generate a table with information about all the SNVs in the DDST
    """
    if len(DDST) == 0:
        return pd.DataFrame()

    cdb = pd.concat(SNPtables)

    # Make a base table with information from the SNP tables
    Bdb = cdb[['position', 'ref_base']].drop_duplicates().set_index('position').sort_index()

    # Make a table of class counts
    class_options = ['DivergentSite', 'SNS', 'SNV', 'con_SNV', 'pop_SNV']
    ccdb = cdb.groupby('position')['class'].value_counts().to_frame().rename(
        columns={'class': 'count'}).reset_index().pivot('position', 'class', 'count').fillna(0).reset_index()
    for c in class_options:
        if c not in ccdb.columns:
            ccdb[c] = 0
    ccdb = ccdb[['position'] + class_options]
    ccdb = ccdb.astype({c: int for c in class_options})
    ccdb = ccdb.rename(columns={c: c + '_count' for c in class_options})
    ccdb = ccdb.set_index('position')

    # Make a table of var_base options
    vdb = cdb.groupby('position')['con_base'].unique().to_frame().rename(columns={'con_base':'sample_con_bases'})
    vdb['sample_con_bases'] = vdb['sample_con_bases'].astype(str)

    # Make a depth table summarizing depth for all samples
    Ddb = DDST.groupby(level=[1]).sum()
    Ddb['depth'] = Ddb['A'] + Ddb['C'] + Ddb['G'] + Ddb['T']


    # Make a table with sample detection numbers
    table = defaultdict(list)
    for SNV, db in DDST.groupby(level=[1]):
        table['position'].append(SNV)
        table['sample_detections'].append(len(db[(db['A'] > 0) | (db['C'] > 0) | (db['T'] > 0) | (db['G'] > 0)]))
        table['sample_5x_detections'].append(len(db[(db['A'] + db['C'] + db['T'] + db['G']) >= 5]))
    DEdb = pd.DataFrame(table).set_index('position')

    # Merge
    Mdb = pd.merge(Ddb, Bdb, left_index=True, right_index=True).join(DEdb).join(ccdb).join(vdb).astype(
                    {'A':int, 'C':int, 'T':int, 'G':int, 'depth':int,
                     'sample_detections':int}).sort_index()

    # Add more calculations
    Mdb['con_base'] = Mdb.apply(calc_con_base, axis=1)
    Mdb['var_base'] = Mdb.apply(calc_var_base, axis=1)

    return Mdb

P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}  # base -> position
C2P = {0:'A', 1:'C', 2:'T', 3:'G'} # position -> base
def calc_var_base(row):
    """
    Calcualte variant base from a row with the counts as columns
    """
    # calculate var_base
    counts_temp = [row['A'], row['C'], row['T'], row['G']]
    counts_temp[P2C[row['con_base']]] = 0
    var_base = C2P[list(counts_temp).index(
        sorted(counts_temp)[-1])]  # this fixes the var_base = ref_base error when there's a tie - alexcc 5/8/2019
    return var_base

def calc_con_base(row):
    """
    calculate the consensus base
    """
    counts = [row['A'], row['C'], row['T'], row['G']]
    return C2P[np.argmax(counts)]
