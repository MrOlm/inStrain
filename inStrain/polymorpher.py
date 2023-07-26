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

    def __init__(self, SC_objects, name2bam, **kwargs):
        """
        Main constructor; to be called from within compare_controller.run_auxillary_processing

        *Right now this is set up for a single scaffold*
        """
        self.name2bam_loc = name2bam
        self.SC_objects = SC_objects

        # Load SNPtables
        name2snpTable = {}
        name2scaffs = defaultdict(list)
        name2isp = {}
        for sc in self.SC_objects:
            for isp, name in zip(sc.profiles, sc.names):
                name2scaffs[name].append(sc.scaffold)
                if name in name2snpTable:
                    continue
                db = isp.get('cumulative_snv_table').rename(
                    columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                             'baseCoverage': 'position_coverage'})
                db = db[db['cryptic'] == False]
                db['scaffold'] = db['scaffold'].astype(str)
                if 'mm' in list(db.columns):
                    db = db.sort_values('mm').drop_duplicates(subset=['scaffold', 'position'], keep='last').sort_index().drop(columns=['mm'])
                name2snpTable[name] = db
                name2isp[name] = isp

        self.name2isp = name2isp
        self.name2snpTable = name2snpTable
        self.name2scaffs = name2scaffs
        self.verify_input()

    def verify_input(self):
        """
        Check that input is OK
        """
        # 1) Make sure you don't have SNPtables on mm level
        for name, t in self.name2snpTable.items():
            assert(len(t) == len(t.drop_duplicates(subset=['scaffold', 'position'])))

    def main(self):
        """
        The main method
        """
        # Extract a list of SNV sites
        self.name2scaff2locs, self.scaff2all_SNVs = extract_SNV_positions(self.name2snpTable, self.name2scaffs)

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
        inStrain.logUtils.log_checkpoint("Compare", "PullSNVsFromBAMs", "start")

        scaff2name2position2counts = {}
        num_to_run = estimate_timing(self.name2scaff2locs)
        pbar = tqdm(desc='Pulling SNVs from BAMs', total=num_to_run)
        for name, scaff2locs in self.name2scaff2locs.items():

            # Load .bam-level cache
            bam_loc = self.name2bam_loc[name]
            Rdic = self.name2isp[name].get("Rdic")

            for scaff, locs in scaff2locs.items():
                if scaff not in scaff2name2position2counts:
                    scaff2name2position2counts[scaff] = {}

                scaff2name2position2counts[scaff][name] = extract_SNVS_from_bam(bam_loc, Rdic, locs, scaff)
                pbar.update(1)

            # Delete .bam-level cache
            del Rdic

        self.scaff2name2position2counts = scaff2name2position2counts
        pbar.close()

        inStrain.logUtils.log_checkpoint("Compare", "PullSNVsFromBAMs", "end")

    def make_pooled_SNV_table(self):
        """
        Using attributes of this object, create an SNV table
        """
        self.DSTdb = genreate_pooled_SNV_table(self.name2snpTable, self.name2scaffs, self.scaff2name2position2counts, self.scaff2all_SNVs)

    def make_pooled_SNV_info_table(self):
        """
        Using the pooled SNV table, create a kind of "summary" object
        """
        self.PMdb = generate_pooled_SNV_summary_table(self.DSTdb, self.name2snpTable, self.name2scaffs)

"""
Below is the deprecated ".dev2" version of the PoolController
"""
# class PoolController(object):
#     """
#     Main controller of the polymorpher pooling commands
#     """
#
#     def __init__(self, SNPtables, names, name2Rdic, name2bam_loc, scaffold, **kwargs):
#         """
#         Main constructor; to be called from within compare
#
#         *Right now this is set up for a single scaffold*
#         """
#         self.scaffold = scaffold
#
#         # Edit the SNPtables to remove mm
#         self.SNPtables = [x.sort_values('mm')\
#             .drop_duplicates(subset=['scaffold', 'position'], keep='last')\
#             .sort_index().drop(columns=['mm']) if 'mm' in list(x.columns) else x for x in SNPtables]
#
#         self.names = names
#         self.name2Rdic = name2Rdic
#         self.name2bam_loc = name2bam_loc
#
#         self.verify_input()
#
#     def verify_input(self):
#         """
#         Check that input is OK
#         """
#         # 1) Make sure you don't have SNPtables on mm level
#         for t in self.SNPtables:
#             assert(len(t) == len(t.drop_duplicates(subset=['scaffold', 'position'])))
#
#     def main(self):
#         """
#         The main method
#         """
#         # Extract a list of SNV sites
#         self.name2locs, self.all_SNVs = extract_SNV_positions(self.SNPtables, self.names)
#
#         # Pull needed SNV sites from .bam files
#         self.pull_SNVS_from_bams()
#
#         # Run processing to create the pooled table
#         self.make_pooled_SNV_table()
#
#         # Run processing to make the SNV info table
#         self.make_pooled_SNV_info_table()
#
#
#     def pull_SNVS_from_bams(self):
#         """
#         Using attributes of this object, pull the needed SNVs from .bam files
#         """
#         name2position2counts = {}
#         for name in self.names:
#             name2position2counts[name] = extract_SNVS_from_bam(self.name2bam_loc[name], self.name2Rdic[name], self.name2locs[name], self.scaffold)
#         self.name2position2counts = name2position2counts
#
#     def make_pooled_SNV_table(self):
#         """
#         Using attributes of this object, create an SNV table
#         """
#         self.DDST = genreate_pooled_SNV_table(self.SNPtables, self.names, self.name2position2counts, self.all_SNVs)
#
#     def make_pooled_SNV_info_table(self):
#         """
#         Using the pooled SNV table, create a kind of "summary" object
#         """
#         self.PST = generate_pooled_SNV_summary_table(self.DDST, self.SNPtables, self.names)


# def extract_SNVS_from_bam(bam_loc, R2M, positions, scaffold, **kwargs):
#     """
#     From a .bam file and a list of SNV locations, extract and return counts
#
#     MY GOD IS PYSAM A PAIN - the reason I have the "start" and "stop" so weird on that iterator is because without it, this fails
#     """
#
#     logging.debug(f"Extracting {len(positions)} SNVs from {scaffold} {bam_loc}")
#
#     # Initilize .bam
#     bam_init = pysam.AlignmentFile(bam_loc)
#
#     # Initialize the object to return
#     position2counts = {}
#
#     # Iterate positions
#     for pos in positions:
#         pileupcolumn = None
#         p = pos
#
#         # Get the pileup for this position
#         try:
#             iter = bam_init.pileup(scaffold, truncate=True, max_depth=100000,
#                                   stepper='nofilter', compute_baq=True,
#                                   ignore_orphans=True, ignore_overlaps=True,
#                                   min_base_quality=30, start=max(p - 1, 0), stop=p + 1)
#             for x in iter:
#                 if x.pos == pos:
#                     pileupcolumn = x
#                     break
#
#             # This means there aren't any reads mapping here, but pysam doesn't return it for whatever reason
#             if pileupcolumn is None:
#                 broken = True
#             else:
#                 assert p == pileupcolumn.pos, pileupcolumn
#                 broken = False
#
#         # This is an error thrown by pysam
#         except ValueError:
#             logging.error("scaffold {0} position {2} is not in the .bam file {1}!".format(scaffold, bam_loc, p))
#             continue
#
#         if broken:
#             counts = np.zeros(4, dtype=int)
#             logging.debug("scaffold {0} position {2} had a funny pysam bug {1}. Should be OK, but just wanted to let you know".format(scaffold, bam_loc, p))
#
#         else:
#             # Process the pileup column
#             MMcounts = inStrain.profile.profile_utilities.get_base_counts_mm(pileupcolumn, R2M)
#             counts = inStrain.profile.profile_utilities.mm_counts_to_counts(MMcounts)
#
#         position2counts[pos] = counts
#
#     return position2counts


def extract_SNVS_from_bam(bam_loc, R2M, positions, scaffold, **kwargs):
    """
    From a .bam file and a list of SNV locations, extract and return counts

    Setting up a single iterator like this really improves speed
    """
    logging.debug(f"Extracting {len(positions)} SNVs from {scaffold} {bam_loc}")

    if len(positions) == 0:
        return {}

    # Initilize .bam
    bam_init = pysam.AlignmentFile(bam_loc)

    # Initiize iterator
    biter = bam_init.pileup(scaffold, truncate=True, max_depth=100000,
                            stepper='nofilter', compute_baq=True,
                            ignore_orphans=True, ignore_overlaps=True,
                            min_base_quality=30, start=max(min(positions) - 1, 0), stop=max(positions) + 1)

    # Initialize the object to return
    position2counts = {}

    # Iterate positions
    pset = set(positions)
    for pilecol in biter:

        # You have it
        if pilecol.pos in pset:
            position2counts[pilecol.pos] = get_pooling_counts(pilecol, R2M)

    # Fill in blanks
    for p in pset - set(position2counts.keys()):
        position2counts[p] = np.zeros(4, dtype=int)

    return position2counts


def get_pooling_counts(pileupcolumn, R2M):
    MMcounts = inStrain.profile.profile_utilities.get_base_counts_mm(pileupcolumn, R2M)
    counts = inStrain.profile.profile_utilities.mm_counts_to_counts(MMcounts)
    return counts

def extract_SNV_positions(name2SNPtables, name2scaffs):
    """
    From a list of tables and names, figure out how many SNVs need to be extracted from .bam files
    """
    # Loop once to get all locations
    scaff2all = {}
    for name, ssdb in name2SNPtables.items():
        if len(ssdb) == 0:
            continue

        for scaff, sdb in ssdb.groupby('scaffold'):
            if scaff not in name2scaffs[name]:
                continue

            if scaff in scaff2all:
                scaff2all[scaff] = scaff2all[scaff].union(set(sdb['position']))
            else:
                scaff2all[scaff] = set(sdb['position'])

    # Loop again to get all needed locations
    name2scaff2locs = {}
    for name, ssdb in name2SNPtables.items():

        if len(ssdb) == 0:
            continue

        # Do an initial loop with "groupby"
        got_scaffs = []
        for scaff, sdb in ssdb.groupby('scaffold'):
        #for scaff in name2scaffs[name]:
            #sdb = ssdb[ssdb['scaffold'] == scaff]

            if scaff not in name2scaffs[name]:
                continue

            if name not in name2scaff2locs:
                name2scaff2locs[name] = {}

            if len(sdb) > 0:
                name2scaff2locs[name][scaff] = scaff2all[scaff] - set(sdb['position'])
            elif scaff in scaff2all:
                name2scaff2locs[name][scaff] = scaff2all[scaff]
            else: # This means you're comparing the scaffold, but it has no SNVs in any sample
                pass
            got_scaffs.append(scaff)

        # Do a second loop to catch stragglers
        for scaff in set(name2scaffs[name]) - set(got_scaffs):

            if name not in name2scaff2locs:
                name2scaff2locs[name] = {}

            if scaff not in scaff2all:
                continue

            name2scaff2locs[name][scaff] = scaff2all[scaff]

    return name2scaff2locs, scaff2all
    #
    # name2locs = {}
    # ALL_SNPS = set()
    #
    # # Loop once to get all locations
    # for sdb in SNPtables:
    #     if len(sdb) == 0:
    #         continue
    #
    #     assert len(sdb['scaffold'].unique()) == 1
    #     ALL_SNPS = ALL_SNPS.union(set(sdb['position']))
    #
    # # Loop once to get all needed locations
    # for sdb, name in zip(SNPtables, names):
    #     if len(sdb) > 0:
    #         name2locs[name] = ALL_SNPS - set(sdb['position'])
    #     else:
    #         name2locs[name] = ALL_SNPS
    #
    # return name2locs, ALL_SNPS

def genreate_pooled_SNV_table(name2snpTable, name2scaffs, scaff2name2position2counts, scaff2all_SNVs):
    #def genreate_pooled_SNV_table(SNPtables, names, name2position2counts, all_SNVs):
    """
    Generate the "DDST" (data-dense SNV table) and "Info table" that can be used for further analysis
    """
    DDSTs = []
    order = []
    for scaff, name2position2counts in scaff2name2position2counts.items():
        dbs = []
        names = []
        order.append(scaff)
        for name, ori_table in name2snpTable.items():
            if scaff not in name2scaffs[name]:
                continue

            # Make a table for this sample from pulled SNVs
            if name in name2position2counts:
                position2counts = name2position2counts[name]
                db = pd.DataFrame.from_dict(position2counts, orient='index', columns=['A', 'C', 'T', 'G'])
            else:
                msg = f"{name} had no SNVs pulled for {scaff}. Has {len(ori_table[ori_table['scaffold'] == scaff])} ori, of {len(set(scaff2all_SNVs[scaff]))} expected"
                logging.debug(msg)
                db = pd.DataFrame()

            ori_table = ori_table[ori_table['scaffold'] == scaff]

            # Merge with existing table
            if len(ori_table) > 0:
                sdb = pd.concat([db, ori_table[['position', 'A', 'C', 'T', 'G']].set_index('position')]).sort_index()
            else:
                sdb = db.sort_index()

            # Verify
            assert set(sdb.index) == set(scaff2all_SNVs[scaff])
            dbs.append(sdb)
            names.append(name)

        DDST = pd.concat(dbs, keys=names)
        DDSTs.append(DDST)

    # Merge DDSTs into one table
    from pandas.api.types import CategoricalDtype
    scaff_cat = CategoricalDtype(order)

    # Add this category to all the dataframes
    for d, s in zip(DDSTs, order):
        d['scaffold'] = s
        d['scaffold'] = d['scaffold'].astype(scaff_cat)

    DSTdb = pd.concat(DDSTs)

    return DSTdb


    # dbs = []
    # for name, ori_table in zip(names, SNPtables):
    #     # Make a table for this sample from pulled SNVs
    #     position2counts = name2position2counts[name]
    #     db = pd.DataFrame.from_dict(position2counts, orient='index', columns=['A', 'C', 'T', 'G'])
    #
    #     # Merge with existing table
    #     if len(ori_table) > 0:
    #         sdb = pd.concat([db, ori_table[['position', 'A', 'C', 'T', 'G']].set_index('position')]).sort_index()
    #     else:
    #         sdb = db.sort_index()
    #
    #     # Verify
    #     assert set(sdb.index) == set(all_SNVs)
    #     dbs.append(sdb)
    #
    # DDST = pd.concat(dbs, keys=names)
    # return DDST


def generate_pooled_SNV_summary_table(DSTdb, name2snpTable, name2scaffs):
    """
    Generate a table with information about all the SNVs in the DDST
    """
    # def generate_pooled_SNV_summary_table(DDST, SNPtables, names):
    if len(DSTdb) == 0:
        return pd.DataFrame()

    Mdbs = []
    order = []
    pbar = tqdm(desc='Generating pooled SNV summary table', total=len(pd.unique(DSTdb['scaffold'])))
    for scaff, DDST in DSTdb.groupby('scaffold'):
        order.append(scaff)

        # Get the base table
        cdb = pd.concat([snptable[snptable['scaffold'] == scaff] for name, snptable in name2snpTable.items() if scaff in name2scaffs[name]])

        # Make a base table with information from the SNP tables
        Bdb = cdb[['position', 'ref_base']].drop_duplicates().set_index('position').sort_index()

        # Make a table of class counts
        class_options = ['DivergentSite', 'SNS', 'SNV', 'con_SNV', 'pop_SNV']
        ccdb = cdb.groupby('position')['class'].value_counts().to_frame().rename(
            columns={'class': 'count'}).reset_index().pivot(index='position', columns='class', values='count').fillna(0).reset_index()
        for c in class_options:
            if c not in ccdb.columns:
                ccdb[c] = 0
        ccdb = ccdb[['position'] + class_options]
        ccdb = ccdb.astype({c: int for c in class_options})
        ccdb = ccdb.rename(columns={c: c + '_count' for c in class_options})
        ccdb = ccdb.set_index('position')

        # Make a table of var_base options
        vdb = cdb.groupby('position')['con_base'].unique().to_frame().rename(columns={'con_base': 'sample_con_bases'})
        vdb['sample_con_bases'] = vdb['sample_con_bases'].astype(str)

        # Make a depth table summarizing depth for all samples
        Ddb_list = []
        index_list = []
        for myindex, indexDF in DDST.groupby(level=[1]):
            Ddb_list.append(indexDF[['A','C','T','G']].sum(axis=0))
            index_list.append(myindex)
        Ddb = pd.DataFrame(Ddb_list)
        Ddb = Ddb.set_index(pd.Index(index_list))
        Ddb['scaffold'] = scaff
        Ddb['depth'] = Ddb['A'] + Ddb['C'] + Ddb['G'] + Ddb['T']

        # Make a table with sample detection numbers
        xdb = DDST[(DDST['A'] + DDST['C'] + DDST['T'] + DDST['G']) >= 5].groupby(level=[1])['A'].agg('count').to_frame() \
            .rename(columns={'A': 'sample_5x_detections'})
        xdb2 = DDST[(DDST['A'] > 0) | (DDST['C'] > 0) | (DDST['T'] > 0) | (DDST['G'] > 0)].groupby(level=[1])['A'].agg(
            'count').to_frame() \
            .rename(columns={'A': 'sample_detections'})
        DEdb = pd.merge(xdb, xdb2, left_index=True, right_index=True)

        # Merge
        Mdb = pd.merge(Ddb, Bdb, left_index=True, right_index=True).join(DEdb).join(ccdb).join(vdb).astype(
            {'A': int, 'C': int, 'T': int, 'G': int, 'depth': int,
             'sample_detections': int}).sort_index()

        # Add more calculations
        Mdb['con_base'] = Mdb.apply(calc_con_base, axis=1)
        Mdb['var_base'] = Mdb.apply(calc_var_base, axis=1)

        Mdbs.append(Mdb)
        pbar.update(1)
    pbar.close()
    # Parse
    from pandas.api.types import CategoricalDtype
    scaff_cat = CategoricalDtype(order)

    # Add this category to all the dataframes
    for d, s in zip(Mdbs, order):
        d['scaffold'] = s
        d['scaffold'] = d['scaffold'].astype(scaff_cat)

    PMdb = pd.concat(Mdbs).astype(
        {'A': int, 'C': int, 'T': int, 'G': int, 'depth': int,
         'sample_detections': int, 'DivergentSite_count': int, 'SNS_count': int, 'SNV_count': int,
         'con_SNV_count': int, 'pop_SNV_count': int, 'sample_5x_detections': int})
    return PMdb

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

def estimate_timing(name2scaff2locs):
    scaffs = 0
    SNVs = 0
    for name, scaff2locs in name2scaff2locs.items():
        for scaff, locs in scaff2locs.items():
            scaffs += 1
            SNVs += len(locs)

    # Convert to minutes
    time = scaffs / 60

    logging.info(f"Pulling {SNVs} SNVs from {scaffs} scaffolds. This should take ~ {time:.1f} min in total")

    return scaffs
