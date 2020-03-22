#!/usr/bin/env python

import os
import csv
import sys
import glob
import logging
import argparse
import traceback

import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import concurrent.futures
from concurrent import futures
from inStrain import SNVprofile
from collections import defaultdict

import inStrain.SNVprofile
import inStrain.controller
import inStrain.profileUtilities
import inStrain.quickProfile

# VERSION 0.1.1
# 10/12/19
# Made mm work better

# VERSION 0.1.0
# 10/12/19
# Lots of little changes to make this work as a script

# VERSION 0.0.9
# 9/27/19
# RC now is mm-afied

# VERSION 0.0.8
# 9/23/19
# Fix a bug with read comparer loading

# VERSION 0.0.7
# 8/28/19
# Make read comparer compatible with inStrain v0.6

# VERSION 0.0.6
# 8/9/19
# Add read comparer

# VERSION 0.0.5
# 6/5/19
# Make the mm-afied scaffold table

# VERSION 0.0.4
# 6/4/19
# Fixed a bug with read reports

# VERSION 0.0.3
# 4/30/19
# Add expected breadth

# VERSION 0.0.2
# 4/29/19
# Many edits

# VERSION 0.0.1
# 4/29/19
# Created from https://biotite.berkeley.edu/j/user/mattolm/notebooks/NIH_Infants/Twins/TwinsStudy_4_typeStrainAudit_2_compareMadness_2.ipynb

def genomeWideFromIS(IS, thing, **kwargs):
    stb = IS.get('scaffold2bin')
    b2l = IS.get('bin2length')
    mm_level = kwargs.get('mm_level', False)

    if stb == None:
        logging.error("This IS object does not have an .stb file; cant use it")
        assert False

    if thing == 'scaffold_info':
        db = IS.get('cumulative_scaffold_table')
        gdb = _add_stb(db, stb)
        assert _validate_b2l(gdb, b2l)
        return _genome_wide_si_2(gdb, stb, b2l, **kwargs)

    if thing == 'read_report':
        db = IS.get('read_report')
        db = db[db['scaffold'] != 'all_scaffolds']
        gdb = _add_stb(db, stb)
        return _genome_wide_rr(gdb, stb, **kwargs)

    if thing == 'read_compaer':
        db = IS.get('comparisonsTable')
        gdb = _add_stb(db, stb)
        return _genome_wide_readComparer(gdb, stb, b2l, **kwargs)

    else:
        logging.error('{0} is not possible to get from genomeWideIS'.format(thing))
        assert False

def makeGenomeWide(thing, db, stb, b2l=None, **kwargs):
    '''
    Convert scaffold-level things to be genome-level things
    '''
    assert type(stb) == type({})

    # # Load it if you have to
    # if thing == 'mm_genome_info':
    #     db = _get_mm_db(db)

    # Validate the stb
    gdb = _add_stb(db, stb)

    # Validate the b2l
    if b2l != None:
        if not _validate_b2l(gdb, b2l):
            return

#     if thing == 'mm_genome_info':
#         return _genome_wide_mm_genomeInfo(gdb, stb, b2l, **kwargs)

    if thing == 'read_report':
        return _genome_wide_rr(gdb, stb, **kwargs)

    elif thing == 'scaffold_info':
        return _genome_wide_si_2(gdb, stb, b2l, **kwargs)

    elif thing == 'readComparer':
        return _genome_wide_readComparer(gdb, stb, b2l, **kwargs)

def _add_stb(db, stb):
    gdb = db.copy()

    if len(gdb) == 0:
        logging.error('Error- no scaffolds detected.')
        return

    gdb['genome'] = gdb['scaffold'].map(stb)

    if len(gdb['genome'].dropna().unique()) == 0:
        logging.error('Error- no genomes detected. Example: stb has scaffold {0}, database has scaffold {1}'.format(
                list(stb.keys())[0], db['scaffold'].tolist()[0]))
        return

    if len(gdb[gdb['genome'].isna()]) > 0:
        logging.info("{0:.2f}% of scaffolds have a genome".format((len(gdb[~gdb['genome'].isna()])/\
                                                      len(gdb))*100))
    return gdb

def _validate_b2l(gdb, b2l):
    genomes = list(gdb['genome'].unique())
    missing = [g for g in genomes if g not in b2l]

    if len(missing) > 1:
        logging.error("{0} of {1} genomes are missing from b2l".format(len(missing), len(genomes)))
        return False

    return True

def _genome_wide_si_2(gdb, stb, b2l, **kwargs):
    '''
    This version can handle mm level
    '''
    mm_level = kwargs.get('mm_level', False)

    if not mm_level:
        gdb = gdb.sort_values('mm').drop_duplicates(
                subset=['scaffold'], keep='last')\
                .sort_values('scaffold')
        gdb['mm'] = 0

    table = defaultdict(list)
    for mm in sorted(list(gdb['mm'].unique())):
        Odb = gdb[gdb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['scaffold'], keep='last')
        for genome, df in Odb.groupby('genome'):
            cols = list(df.columns)

            table['mm'].append(mm)
            table['genome'].append(genome)

            # Scaffolds
            table['detected_scaffolds'].append(len(df))
            table['true_scaffolds'].append(len([True for s, b in stb.items() if b == genome]))
            table['true_length'].append(int(b2l[genome]))

            # The summing columns
            for col in ['SNPs', 'Referece_SNPs', 'BiAllelic_SNPs', 'MultiAllelic_SNPs', 'consensus_SNPs', 'population_SNPs']:
                if col in cols:
                    table[col].append(df[col].fillna(0).sum())

            # Weighted average (over total length)
            for col in ['breadth', 'coverage', 'std_cov']:
                table[col].append(sum([x * y for x, y in zip(df[col].fillna(0), df['length'])]) / b2l[genome])

            # Weighted average (over detected scaffold length)
            df['considered_length'] = [x*y for x,y in zip(df['unmaskedBreadth'], df['length'])]
            considered_leng = float(df['considered_length'].sum())

            # To maintain backwards compatibility
            for col in ['mean_clonality']:
                if col not in df.columns:
                    continue
                if considered_leng != 0:
                    table[col].append(sum(x * y for x, y in zip(df[col].fillna(0), df['considered_length'])) / considered_leng)
                else:
                    table[col].append(np.nan)

            for col in ['mean_microdiverstiy', 'rarefied_mean_microdiversity']:
                if col not in df.columns:
                    continue
                if considered_leng != 0:
                    table[col].append(sum(x * y for x, y in zip(df[col].fillna(0), df['considered_length'])) / considered_leng)
                else:
                    table[col].append(np.nan)

            # ANI
            if 'consensus_SNPs' in cols:
                if considered_leng != 0:
                    table['conANI'].append((considered_leng - df['consensus_SNPs'].sum()) / considered_leng)
                    table['popANI'].append((considered_leng - df['population_SNPs'].sum()) / considered_leng)
                else:
                    table['conANI'].append(0)
                    table['popANI'].append(0)
            else:
                if considered_leng != 0:
                    table['ANI'].append((considered_leng - df['SNPs'].sum()) / considered_leng)
                else:
                    table['ANI'].append(0)
            #table['detected_scaff_length'].append(df['length'].sum())

            table['unmaskedBreadth'].append(considered_leng / b2l[genome])
            table['expected_breadth'].append(estimate_breadth(table['coverage'][-1]))

    db = pd.DataFrame(table)

    # Add back microdiversity
    if (('mean_microdiversity' not in df.columns) & ('mean_clonality' in df.columns)):
        db['mean_microdiversity'] = 1 - db['mean_clonality']

    if not mm_level:
        del db['mm']

    return db


# def _get_mm_db(obj):
#     '''
#     From an inStrain object (currently a pickle), get the raw scaffold table
#     '''
#     s = inStrain.SNVprofile()
#     s.load(obj)
#     db = s.cumulative_scaffold_table
#
#     return db

def _genome_wide_rr(gdb, stb, **kwrags):
    '''
    THIS REALLY SHOULD USE SCAFFOLD LENGTH FOR A WEIGHTED AVERAGE!
    '''
    table = defaultdict(list)
    for genome, df in gdb.groupby('genome'):
        table['genome'].append(genome)
        table['scaffolds'].append(len(df))
        for col in [c for c in list(df.columns) if c not in ['scaffold', 'genome']]:
            if col.startswith('pass') or col.startswith('unfiltered_') or col.startswith('filtered'):
                table[col].append(df[col].sum())
            else:
                table[col].append(df[col].mean())

    return pd.DataFrame(table)

# def _genome_wide_si(gdb, stb, b2l, **kwargs):
#     table = defaultdict(list)
#     for genome, df in gdb.groupby('genome'):
#         table['genome'].append(genome)

#         # Scaffolds
#         table['detected_scaffolds'].append(len(df))
#         table['true_scaffolds'].append(len([True for s, b in stb.items() if b == genome]))
#         table['true_length'].append(int(b2l[genome]))

#         # The summing columns
#         for col in ['SNPs']:
#             table[col].append(df[col].sum())

#         # Weighted average (over total length)
#         for col in ['breadth', 'coverage', 'std_cov']:
#             table[col].append(sum(x * y for x, y in zip(df[col], df['length'])) / b2l[genome])

#         # Weighted average (over detected scaffold length)
#         df['considered_length'] = [x*y for x,y in zip(df['unmaskedBreadth'], df['length'])]
#         considered_leng = df['considered_length'].sum()
#         for col in ['mean_clonality']:
#             table[col].append(sum(x * y for x, y in zip(df[col], df['considered_length'])) / considered_leng)

#         # ANI
#         if considered_leng != 0:
#             table['ANI'].append((considered_leng - df['SNPs'].sum()) / considered_leng)
#         else:
#             table['ANI'].append(0)
#         #table['detected_scaff_length'].append(df['length'].sum())

#         table['unmaskedBreadth'].append(considered_leng / b2l[genome])
#         table['expected_breadth'].append(estimate_breadth(table['coverage'][-1]))

#     return pd.DataFrame(table)

# def _genome_wide_mm_genomeInfo(gdb, stb, b2l, **kwargs):
#     '''
#     JupyterNotebooks/Infant_Eukaryotes/Euks_13_mappingListGamma_bams_2_load_v2_2.ipynb
#     '''
#     gdb['genome'] = gdb['scaffold'].map(stb)
#     gdb['considered_length'] = [x*y for x,y in zip(gdb['unmaskedBreadth'], gdb['length'])]

#     # Add blanks for scaffolds that aren't there
#     #bdb = pd.DataFrame({s:})

#     table = defaultdict(list)
#     for genome, mm, db in interate_sdb_mm(gdb, s2l=None, stb=stb):
#         table['genome'].append(genome)
#         table['mm'].append(mm)

#         # Sum columns
#         for col in ['bases_w_0_coverage', 'length']:
#             table[col].append(db[col].sum())

#         # Max columns
#         for col in ['SNPs']:
#             table[col].append(db[col].max())

#         # Weighted average columns
#         for col in ['breadth', 'coverage', 'std_cov', 'unmaskedBreadth']:
#             table[col].append(sum(x * y for x, y in zip(db[col], db['length'])) / sum(db['length']))

#         # Special weighted average
#         db['considered_length'] = [x*y for x,y in zip(db['unmaskedBreadth'], db['length'])]
#         considered_leng = db['considered_length'].sum()

#         if considered_leng != 0:
#             table['ANI'].append((considered_leng - db['SNPs'].sum()) / considered_leng)
#         else:
#             table['ANI'].append(0)

#         # Special
# #         table['max_cov'].append(db['max_cov'].max())
# #         table['min_cov'].append(db['min_cov'].min())

#     return pd.DataFrame(table)

def _genome_wide_readComparer(gdb, stb, b2l, **kwargs):
    '''
    Now this does work on the genome level
    '''
    mm_level = kwargs.get('mm_level', False)

    if not mm_level:
        gdb = gdb.sort_values('mm').drop_duplicates(
                subset=['scaffold', 'name1', 'name2'], keep='last')\
                .sort_values('scaffold')
        gdb['mm'] = 0
    table = defaultdict(list)

    # Backwards compatibility
    #gdb = gdb.rename(columns={'ANI':'popANI'})

    for mm in sorted(list(gdb['mm'].unique())):
        Odb = gdb[gdb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['scaffold', 'name1', 'name2'], keep='last')
        for things, db in Odb.groupby(['genome', 'name1', 'name2']):
#         for genome, df in Odb.groupby('genome'):

#     for things, db in gdb.groupby(['genome', 'name1', 'name2', 'mm']):
            genome, name1, name2 = things

            table['genome'].append(genome)
            table['name1'].append(name1)
            table['name2'].append(name2)
            table['mm'].append(mm)

            # Weighted average columns
            tcb = db['compared_bases_count'].sum()
            for col in ['coverage_overlap']:
                if tcb == 0:
                    table[col].append(np.nan)
                else:
                    table[col].append(sum([x * y for x, y in zip(db[col], db['compared_bases_count'])]) / tcb)

            # Sum columns
            for col in ['compared_bases_count', 'consensus_SNPs', 'population_SNPs']:
                if col in list(db.columns):
                    table[col].append(db[col].sum())

            # Special ANI column
            for col in ['ANI', 'popANI', 'conANI']:
                if col in list(db.columns):
                    if tcb == 0:
                        table[col].append(np.nan)
                    else:
                        table[col].append(sum(a * c  if a == a else 0 for a, c in
                                            zip(db[col], db['compared_bases_count'])) / tcb)

            # Percent compared
            if b2l != None:
                table['percent_compared'].append(tcb / b2l[genome])

    db = pd.DataFrame(table)

    if not mm_level:
        del db['mm']

    return db

def _backfill_blanks(db, s2l):
    scaffs = list(set(s2l.keys()) - set(db['scaffold'].unique()))
    bdb = pd.DataFrame({'scaffold':scaffs, 'length':[s2l[s] for s in scaffs]})

    # make some adjustments
    bdb['bases_w_0_coverage'] = bdb['length']

    # append
    db = db.append(bdb, sort=True)

    # fill 0
    return db.fillna(0)

def interate_sdb_mm(sdb, on='genome', s2l=None, stb=None):
    '''
    For the dataframe, iterate through each mm and a dataframe with ALL scaffolds at that mm level
    (including blanks)
    '''
    for g, db in sdb.groupby(on):
        if s2l == None:
            gs2l = db.set_index('scaffold')['length'].to_dict()
        else:
            gs2l = {s:s2l[s] for s in [x for x, b in stb.items() if b == g]}
        mms = sorted(db['mm'].unique())
        for mm in mms:
            # get all the ones that you can
            dd = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset='scaffold', keep='last')
            #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))

            # backfill with blanks
            dd = _backfill_blanks(dd, gs2l)
            #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))

            yield g, mm, dd

def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage

    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    import numpy as np
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1

class Controller():

    def main(self, args):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        args = self.validate_input(args)
        vargs = vars(args)
        mm_level = vargs.get('mm_level', False)
        IS = vargs.pop('IS')

        # Read the scaffold to bin file
        logging.debug('Loading scaffold to bin')
        stb = load_scaff2bin(args.stb, IS)

        # Make bin to length
        try:
            stl = IS.get('scaffold2length')
            b2l = {}
            for scaffold, bin in stb.items():
                if bin not in b2l:
                    b2l[bin] = 0
                if scaffold in stl:
                    b2l[bin] += stl[scaffold]
        except:
            logging.error('Could not make bin to length')
            b2l = None

        # Figure out output base
        out_base = IS.get_location('output') + os.path.basename(IS.get('location')) + '_'

        # Store the scaffold to bin and bin to length file
        IS.store('scaffold2bin', stb, 'dictionary', 'Dictionary of scaffold 2 bin')
        IS.store('bin2length', b2l, 'dictionary', 'Dictionary of bin 2 total length')

        # Do genome info
        try:
            gdb = genomeWideFromIS(IS, 'scaffold_info', mm_level=mm_level)
            gdb.to_csv(out_base + 'genomeWide_scaffold_info.tsv', index=False, sep='\t')
        except:
            logging.debug("GenomeWide scaffold info failed ", exc_info=1)
            #traceback.print_exc()

        # Do read report
        try:
            gdb = genomeWideFromIS(IS, 'read_report', mm_level=mm_level)
            gdb.to_csv(out_base + 'genomeWide_read_report.tsv', index=False, sep='\t')
        except:
            logging.debug("GenomeWide read report failed ", exc_info=1)
            #traceback.print_exc()

        # Do read comparer
        try:
            gdb = genomeWideFromIS(IS, 'read_compaer', mm_level=mm_level)
            gdb.to_csv(out_base + 'genomeWide_compare.tsv', index=False, sep='\t')
        except:
            logging.debug("GenomeWide compare failed ", exc_info=1)
            # traceback.print_exc()

    def validate_input(self, args):
        '''
        Validate and mess with the arguments a bit
        '''
        # Make sure the IS object is OK
        assert os.path.exists(args.IS)
        args.IS = inStrain.SNVprofile.SNVprofile(args.IS)

        # Set up the logger
        log_loc = args.IS.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)

        return args

def _genome():
    return 'genome'

def load_scaff2bin(input_stb, IS):
    '''
    From the input, load an .stb. The input can be a lot of things, though
    '''
    # Check if there is nothing there
    if input_stb == []:
        s2l = IS.get('scaffold2length')
        stb = {}
        for scaffold in list(s2l.keys()):
            stb[scaffold] = 'all_scaffolds'
        logging.info('Scaffold to bin will consider all scaffolds the same genome')
    else:
        stb = None

    # Check if this is a .fasta file
    if stb == None:
        try:
            stb = gen_stb(input_stb)
            logging.info('Scaffold to bin was made using .fasta files')
        except:
            stb = None
            #traceback.print_exc()

    # Check if this is a regular stb file
    if (stb == None) & (len(input_stb) == 1):
        try:
            stb = inStrain.quickProfile.parse_stb(input_stb[0])
            logging.info('Scaffold to bin was made using .stb file')
        except:
            stb = None

    if stb == None:
        logging.error('Could not load the scaffold to bin file!')
        assert False

    return stb

def gen_stb(fastas):
    stb = {}
    for fasta in fastas:
        bin = os.path.basename(fasta)
        for seq_record in SeqIO.parse(fasta, "fasta"):
            id = str(seq_record.id).strip()
            stb[id] = bin
    if len(stb.keys()) == 0:
        raise Exception
    return stb
