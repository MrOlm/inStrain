#!/usr/bin/env python

import os
import csv
import sys
import glob
import scipy
import logging
import warnings
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
import inStrain.profile.profile_utilities
import inStrain.quickProfile
import inStrain.irep_utilities
import inStrain.logUtils

class Controller():

    def main(self, args):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        IS, vargs = self.validate_input(args)

        # Load scaffold to bin and bin to length; put in IS
        logging.debug('Loading scaffold to bin')
        self.prepare_genome_wide(IS, vargs)

        # Make the genome-wide information table for a profile object
        object_type = IS.get('object_type')
        if object_type is None:
            if IS.get('comparisonsTable') is None:
                object_type = 'profile'
            else:
                object_type = 'compare'
            # logging.error("Theres no object type! This is bad! I'll assume profile")
            # print(IS)
            # object_type = 'profile'

        if object_type == 'profile':
            # Generate and store the genome level info
            GIdb = genomeLevel_from_IS(IS, **vargs)
            IS.store('genome_level_info', GIdb, 'pandas', 'Table of genome-level information')

            # Store the output
            IS.generate('genome_info', **vargs)
            IS.generate('SNVs')

        else:
            self.prepare_genome_wide(IS, vargs)
            gdb = genomeWideFromIS(IS, 'read_compaer', **vargs)
            out_base = IS.get_location('output') + os.path.basename(IS.get('location')) + '_'
            gdb.to_csv(out_base + 'genomeWide_compare.tsv', index=False, sep='\t')
            # out_base = IS.get_location('output') + os.path.basename(IS.get('location')) + '_'
            # GIdb.to_csv(out_base + 'genome_info.tsv', index=False, sep='\t')

            # Store the SNV output
            # IS.generate('SNVs')

        # Make the genome-wide information table for a compare object
        # elif object_type == 'compare':
        #     pass
        #
        # # Figure out output base
        # out_base = IS.get_location('output') + os.path.basename(IS.get('location')) + '_'
        #
        #
        # # Do genome info
        # try:
        #     gdb = genomeWideFromIS(IS, 'scaffold_info', mm_level=mm_level)
        #     gdb.to_csv(out_base + 'genomeWide_scaffold_info.tsv', index=False, sep='\t')
        # except:
        #     logging.debug("GenomeWide scaffold info failed ", exc_info=1)
        #     #traceback.print_exc()
        #
        # # Do read report
        # try:
        #     gdb = genomeWideFromIS(IS, 'mapping_info', mm_level=mm_level)
        #     gdb.to_csv(out_base + 'genomeWide_mapping_info.tsv', index=False, sep='\t')
        # except:
        #     logging.debug("GenomeWide read report failed ", exc_info=1)
        #     #traceback.print_exc()
        #
        # # Do read comparer
        # try:
        #     gdb = genomeWideFromIS(IS, 'read_compaer', mm_level=mm_level)
        #     gdb.to_csv(out_base + 'genomeWide_compare.tsv', index=False, sep='\t')
        # except:
        #     logging.debug("GenomeWide compare failed ", exc_info=1)
        #     # traceback.print_exc()

    def validate_input(self, args):
        '''
        Validate and mess with the arguments a bit
        '''
        # Make sure the IS object is OK
        assert os.path.exists(args.IS)
        IS = inStrain.SNVprofile.SNVprofile(args.IS)

        # Set up the logger
        log_loc = IS.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)

        # Convert argument type
        vargs = vars(args)
        vargs.pop('IS')

        return IS, vargs

    def prepare_genome_wide(self, IS, vargs):
        '''
        Load and store the scaffold to bin file and bin to length file
        '''
        stb = load_scaff2bin(vargs.get('stb'), IS)

        # Make bin to length
        stl = IS.get('scaffold2length')
        b2l = {}
        for scaffold, bin in stb.items():
            if bin not in b2l:
                b2l[bin] = 0

            if scaffold in stl:
                b2l[bin] += stl[scaffold]
            else:
                logging.error("FAILURE StbError {0} {1} no_length will not be considered as part of the genome".format(scaffold, bin))

        # Store these
        IS.store('scaffold2bin', stb, 'dictionary', 'Dictionary of scaffold 2 bin')
        IS.store('bin2length', b2l, 'dictionary', 'Dictionary of bin 2 total length')

def genomeLevel_from_IS(IS, **kwargs):
    '''
    Calculate genome-wide metrics based on scaffold-level info.

    IS object must already have:
        scaffold2bin
        bin2length

    Arguments:
        IS = InStrain Profile object; must already be initialized

    kwargs:
        mm-level = If True, calculate all metrics on the mm-level
    '''
    inStrain.logUtils.log_checkpoint("GenomeLevel", "genomeLevel_from_IS", "start")

    # Load necessary stuff from IS file
    stb = IS.get('scaffold2bin')
    b2l = IS.get('bin2length')
    s2l = IS.get('scaffold2length')

    # Calculate averaing and summing metrics from the scaffold table
    db = IS.get('cumulative_scaffold_table')
    gdb = _add_stb(db, stb)

    # Handle mm level
    skip_mm_level = kwargs.get('skip_mm_profiling', False)
    if skip_mm_level:
        gdb = gdb.sort_values('mm').drop_duplicates(
                subset=['scaffold'], keep='last')\
                .sort_values('scaffold')
        gdb['mm'] = 1000

    # Calculate genome level scaffold info
    inStrain.logUtils.log_checkpoint("GenomeLevel", "scaffold_info", "start")

    GSI_db = _genomeLevel_scaffold_info_v3(gdb, stb, b2l, **kwargs)

    inStrain.logUtils.log_checkpoint("GenomeLevel", "scaffold_info", "end")

    # Prepare to calculate coverage distribution metrics
    table = defaultdict(list)
    bin2scaffolds = calc_bin2scaffols(stb)
    scaff2sequence = _load_scaff2sequence(IS, **kwargs)

    relevant_genomes = calc_relevant_genomes(GSI_db, IS, **kwargs)
    relevant_scaffolds = set()
    for g in relevant_genomes:
        for s in bin2scaffolds[g]:
            relevant_scaffolds.add(s)
    covT = IS.get('covT', scaffolds=relevant_scaffolds)

    # Figure out total number of mms
    if skip_mm_level:
        mms = [1000]
    else:
        mms = calc_mms(covT)

    # Calculate expensive coverage distribution metrics
    inStrain.logUtils.log_checkpoint("GenomeLevel", "coverage_info", "start")

    EG_db = genomeLevel_coverage_info(covT, bin2scaffolds, relevant_genomes,
                                        s2l, scaff2sequence, mms, **kwargs)

    inStrain.logUtils.log_checkpoint("GenomeLevel", "coverage_info", "end")

    # Calculate genome-level read mapping
    rdb = IS.get('mapping_info')
    rdb = rdb[rdb['scaffold'] != 'all_scaffolds']
    rdb = _add_stb(rdb, stb)

    inStrain.logUtils.log_checkpoint("GenomeLevel", "mapping_info", "start")

    rdb = _genome_wide_rr(rdb, stb, **kwargs)
    rdb = rdb.rename(columns={'reads_filtered_pairs':'filtered_read_pair_count'})
    if 'reads_pass_pairing_filter' in rdb.columns:
        del rdb['reads_pass_pairing_filter']

    inStrain.logUtils.log_checkpoint("GenomeLevel", "mapping_info", "end")

    # Merge
    mdb = pd.merge(GSI_db, EG_db, on=['genome', 'mm'], how='outer')
    mdb = pd.merge(mdb, rdb, on=['genome'], how='left')

    # Calculate genome-level linkage metrics
    ldb = IS.get('raw_linkage_table')
    if len(ldb) > 0:
        if skip_mm_level:
            ldb = ldb.sort_values('mm').drop_duplicates(
                    subset=['scaffold', 'position_A', 'position_B'], keep='last')\
                    .sort_values(['scaffold', 'position_A', 'position_B'])
            gdb['mm'] = 1000

        ldb = _add_stb(ldb, stb)
        if (ldb is not None) and (len(ldb) > 0):

            inStrain.logUtils.log_checkpoint("GenomeLevel", "linkage", "start")

            ldb = _genome_wide_linkage(ldb, stb, mms, **kwargs)

            inStrain.logUtils.log_checkpoint("GenomeLevel", "linkage", "end")


            mdb = pd.merge(mdb, ldb, on=['genome', 'mm'], how='left')


    if skip_mm_level:
        del mdb['mm']

    inStrain.logUtils.log_checkpoint("GenomeLevel", "genomeLevel_from_IS", "end")

    return mdb

def _load_scaff2sequence(IS, **kwargs):
    '''
    Try and get scaff2sequence to GC-correct iRep. Return None if you cant
    '''
    if kwargs.get('fasta', None) is not None:
        fasta_loc = kwargs.get('fasta')
    else:
        fasta_loc = IS.get('fasta_loc')

    if fasta_loc is None:
        logging.error("Do not have a .fasta file; will not GC correct iRep")
        scaff2sequence = None
    else:
        try:
            scaff2sequence = SeqIO.to_dict(SeqIO.parse(fasta_loc, "fasta"))
        except:
            logging.error("Could not load .fasta file {0}; will not GC correct iRep".format(fasta_loc))
            scaff2sequence = None
    return scaff2sequence

def calc_relevant_genomes(GSI_db, IS, **kwargs):
    '''
    Depending on kwrags, determine which scaffolds should have expensive profiling done
    '''
    return set(GSI_db['genome'].tolist())

def genomeLevel_coverage_info(covT, bin2scaffolds, relevant_genomes, s2l,
                                scaff2sequence, mms, **kwargs):
    '''
    THIS IS ALL CALCULATED WITH MASKED EDGES OF SCAFFOLDS
    '''

    dbs = []

    for genome, scaffolds in bin2scaffolds.items():
        if genome not in relevant_genomes:
            continue

        table = defaultdict(list)

        # Sort the scaffolds
        scaffolds = scaffolds.intersection(set(s2l.keys())) # THIS SHOULDN'T HAVE TO BE HERE!
        scaffolds = sorted(scaffolds, key=s2l.get, reverse=True)

        # Try and get the GC_correction
        gc_windows = None
        if scaff2sequence is not None:
            try:
                gc_windows = inStrain.irep_utilities.generate_gc_windows(scaffolds, scaff2sequence, mask_edges=100)
            except:
                pass

        # Set up default iRep
        iRep = np.nan
        iRep_accessory = {'iRep_GC_corrected':np.nan}

        # Iterate for this genome
        for mm in mms:

            # Make a coverage array for this genome in it's entirety
            covs, scaff2genome_index = generate_genome_coverage_array(covT, s2l, order=scaffolds, maxMM=mm, mask_edges=100)

            # Calculate iRep only on mm == 1, or if not on mm level
            if ((mm == 1) or (mm == 1000)):
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        iRep, iRep_accessory = inStrain.irep_utilities.calculate_iRep_from_coverage_array(covs, len(scaffolds), gc_windows)
                    logging.debug("iRep {0} {1} {2}".format(genome, mm, iRep_accessory))
                except:
                    logging.warning("FAILURE iRepError {0} {1}".format(genome, mm))
                    iRep = np.nan
                    iRep_accessory = {'iRep_GC_corrected':np.nan}

            # Calculate medians and other variants
            if len(covs) == 0:
                covs = pd.Series([0])

            table['mm'].append(mm)
            table['genome'].append(genome)

            # Dont print annoying runtime warnings here
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                table['coverage_median'].append(int(np.median(covs)))
                table['coverage_SEM'].append(scipy.stats.sem(covs))
                table['coverage_std'].append(np.std(covs))

        gdb = pd.DataFrame(table)
        gdb['iRep'] = iRep
        gdb.loc[:,'iRep_GC_corrected'] = iRep_accessory['iRep_GC_corrected']
        dbs.append(gdb)

    adb = pd.concat(dbs).reset_index(drop=True)
    return adb

def genomeWideFromIS(IS, thing, **kwargs):
    '''
    This is used by plotting utilities
    '''
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

    if thing == 'mapping_info':
        db = IS.get('mapping_info')
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

    if thing == 'mapping_info':
        return _genome_wide_rr(gdb, stb, **kwargs)

    elif thing == 'scaffold_info':
        return _genome_wide_si_2(gdb, stb, b2l, **kwargs)

    elif thing == 'readComparer':
        return _genome_wide_readComparer(gdb, stb, b2l, **kwargs)

def _add_stb(db, stb, verbose=True):
    gdb = db.copy()

    if len(gdb) == 0:
        logging.error('Error- no scaffolds detected.')
        return

    gdb.loc[:,'genome'] = gdb['scaffold'].map(stb)

    if len(gdb['genome'].dropna().unique()) == 0:
        logging.error('Error- no genomes detected. Example: stb has scaffold {0}, database has scaffold {1}'.format(
                list(stb.keys())[0], db['scaffold'].tolist()[0]))
        return

    if len(gdb[gdb['genome'].isna()]) > 0:
        if verbose:
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
    skip_mm_level = kwargs.get('skip_mm_level', False)

    if skip_mm_level:
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
            table['length'].append(int(b2l[genome]))

            # The summing columns
            for col in ['SNPs', 'Referece_SNPs', 'BiAllelic_SNPs', 'MultiAllelic_SNPs', 'consensus_SNPs', 'population_SNPs']:
                if col in cols:
                    table[col].append(df[col].fillna(0).sum())

            # Weighted average (over total length)
            for col in ['breadth', 'coverage', 'coverage_std']:
                table[col].append(sum([x * y for x, y in zip(df[col].fillna(0), df['length'])]) / b2l[genome])

            # Weighted average (over detected scaffold length)
            df.loc[:,'considered_length'] = [x*y for x,y in zip(df['breadth_minCov'], df['length'])]
            considered_leng = float(df['considered_length'].sum())

            # To maintain backwards compatibility
            for col in ['mean_clonality']:
                if col not in df.columns:
                    continue
                if considered_leng != 0:
                    table[col].append(sum(x * y for x, y in zip(df[col].fillna(0), df['considered_length'])) / considered_leng)
                else:
                    table[col].append(np.nan)

            for col in ['mean_microdiverstiy', 'nucl_diversity_rarefied']:
                if col not in df.columns:
                    continue
                if considered_leng != 0:
                    table[col].append(sum(x * y for x, y in zip(df[col].fillna(0), df['considered_length'])) / considered_leng)
                else:
                    table[col].append(np.nan)

            # ANI
            if 'consensus_SNPs' in cols:
                if considered_leng != 0:
                    table['conANI_reference'].append((considered_leng - df['consensus_SNPs'].sum()) / considered_leng)
                    table['popANI_reference'].append((considered_leng - df['population_SNPs'].sum()) / considered_leng)
                else:
                    table['conANI_reference'].append(0)
                    table['popANI_reference'].append(0)
            else:
                if considered_leng != 0:
                    table['ANI'].append((considered_leng - df['SNPs'].sum()) / considered_leng)
                else:
                    table['ANI'].append(0)
            #table['detected_scaff_length'].append(df['length'].sum())

            table['breadth_minCov'].append(considered_leng / b2l[genome])
            table['breadth_expected'].append(estimate_breadth(table['coverage'][-1]))

    db = pd.DataFrame(table)

    # Add back microdiversity
    if (('nucl_diversity' not in df.columns) & ('mean_clonality' in df.columns)):
        db['nucl_diversity'] = 1 - db['mean_clonality']

    if skip_mm_level:
        del db['mm']

    return db

def _genomeLevel_scaffold_info_v3(gdb, s2b, b2l, **kwargs):
    '''
    Private method for calculating genome-level scaffold info
    '''
    stb = s2b # This is because you also have the kwarg
    table = defaultdict(list)
    for mm in sorted(list(gdb['mm'].unique())):
        # Get scaffold information at this mm
        Odb = gdb[gdb['mm'] <= mm].sort_values('mm')\
                .drop_duplicates(subset=['scaffold'], keep='last')

        for genome, df in Odb.groupby('genome'):
            cols = list(df.columns)

            table['mm'].append(mm)
            table['genome'].append(genome)

            # Scaffolds
            table['detected_scaffolds'].append(len(df))
            table['true_scaffolds'].append(len([True for s, b in stb.items()
                                                                if b == genome]))
            table['length'].append(int(b2l[genome]))

            # The summing columns
            for col in ['SNS_count', 'SNV_count', 'divergent_site_count',
                    'consensus_divergent_sites', 'population_divergent_sites']:
                if col in cols:
                    table[col].append(df[col].fillna(0).sum())

            # Weighted average (over total length)
            for col in ['breadth', 'coverage']:
                table[col].append(sum([x * y for x, y in zip(df[col].fillna(0), df['length'])]) / b2l[genome])

            # Weighted average (over detected scaffold length)
            df['considered_length'] = [x*y for x,y in zip(
                                        df['breadth_minCov'], df['length'])]
            considered_leng = float(df['considered_length'].sum())

            for col in ['nucl_diversity', 'nucl_diversity_rarefied']:
                if col not in df.columns:
                    continue
                if considered_leng != 0:
                    table[col].append(sum(x * y for x, y in zip(df[col].fillna(0), df['considered_length'])) / considered_leng)
                else:
                    table[col].append(np.nan)

            # ANI
            if 'consensus_divergent_sites' in cols:
                if considered_leng != 0:
                    table['conANI_reference'].append((considered_leng - df['consensus_divergent_sites'].sum()) / considered_leng)
                    table['popANI_reference'].append((considered_leng - df['population_divergent_sites'].sum()) / considered_leng)
                else:
                    table['conANI_reference'].append(0)
                    table['popANI_reference'].append(0)

            # Special
            table['breadth_minCov'].append(considered_leng / b2l[genome])
            table['breadth_expected'].append(estimate_breadth(table['coverage'][-1]))

    db = pd.DataFrame(table)
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

def _genome_wide_rr(gdb, s2b, **kwrags):
    '''
    THIS REALLY SHOULD USE SCAFFOLD LENGTH FOR A WEIGHTED AVERAGE!
    '''
    stb = s2b
    table = defaultdict(list)
    for genome, df in gdb.groupby('genome'):
        table['genome'].append(genome)
        #table['scaffolds'].append(len(df))
        for col in [c for c in list(df.columns) if c not in ['scaffold', 'genome']]:
            if len(df[col].dropna()) == 0:
                table['reads_' + col].append(np.nan)
            elif col.startswith('pass') or col.startswith('unfiltered_') or col.startswith('filtered'):
                table['reads_' + col].append(df[col].sum())
            else:
                table['reads_' + col].append(df[col].mean())

    return pd.DataFrame(table)

def _genome_wide_linkage(ldb, s2b, mms, **kwrags):
    '''
    From the raw_linkage_table, calculate average d' and r2 for each genome
    '''
    stb = s2b
    table = defaultdict(list)
    for mm in mms:
        # Get scaffold information at this mm
        Odb = ldb[ldb['mm'] <= mm].sort_values('mm')\
                .drop_duplicates(subset=['scaffold', 'position_A', 'position_B'],
                keep='last')

        if len(Odb) == 0:
            continue

        for genome, df in Odb.groupby('genome'):
            table['genome'].append(genome)
            table['mm'].append(mm)
            table['r2_mean'].append(df['r2'].mean())
            table['d_prime_mean'].append(df['d_prime'].mean())
            table['SNV_distance_mean'].append(df['distance'].mean())
            table['linked_SNV_count'].append(len(df))

    return pd.DataFrame(table)

# def _genome_wide_si(gdb, stb, b2l, **kwargs):
#     table = defaultdict(list)
#     for genome, df in gdb.groupby('genome'):
#         table['genome'].append(genome)

#         # Scaffolds
#         table['detected_scaffolds'].append(len(df))
#         table['true_scaffolds'].append(len([True for s, b in stb.items() if b == genome]))
#         table['length'].append(int(b2l[genome]))

#         # The summing columns
#         for col in ['SNPs']:
#             table[col].append(df[col].sum())

#         # Weighted average (over total length)
#         for col in ['breadth', 'coverage', 'coverage_std']:
#             table[col].append(sum(x * y for x, y in zip(df[col], df['length'])) / b2l[genome])

#         # Weighted average (over detected scaffold length)
#         df['considered_length'] = [x*y for x,y in zip(df['breadth_minCov'], df['length'])]
#         considered_leng = df['considered_length'].sum()
#         for col in ['mean_clonality']:
#             table[col].append(sum(x * y for x, y in zip(df[col], df['considered_length'])) / considered_leng)

#         # ANI
#         if considered_leng != 0:
#             table['ANI'].append((considered_leng - df['SNPs'].sum()) / considered_leng)
#         else:
#             table['ANI'].append(0)
#         #table['detected_scaff_length'].append(df['length'].sum())

#         table['breadth_minCov'].append(considered_leng / b2l[genome])
#         table['breadth_expected'].append(estimate_breadth(table['coverage'][-1]))

#     return pd.DataFrame(table)

# def _genome_wide_mm_genomeInfo(gdb, stb, b2l, **kwargs):
#     '''
#     JupyterNotebooks/Infant_Eukaryotes/Euks_13_mappingListGamma_bams_2_load_v2_2.ipynb
#     '''
#     gdb['genome'] = gdb['scaffold'].map(stb)
#     gdb['considered_length'] = [x*y for x,y in zip(gdb['breadth_minCov'], gdb['length'])]

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
#         for col in ['breadth', 'coverage', 'coverage_std', 'breadth_minCov']:
#             table[col].append(sum(x * y for x, y in zip(db[col], db['length'])) / sum(db['length']))

#         # Special weighted average
#         db['considered_length'] = [x*y for x,y in zip(db['breadth_minCov'], db['length'])]
#         considered_leng = db['considered_length'].sum()

#         if considered_leng != 0:
#             table['ANI'].append((considered_leng - db['SNPs'].sum()) / considered_leng)
#         else:
#             table['ANI'].append(0)

#         # Special
# #         table['max_cov'].append(db['max_cov'].max())
# #         table['min_cov'].append(db['min_cov'].min())

#     return pd.DataFrame(table)

def _genome_wide_readComparer(gdb, s2b, b2l, **kwargs):
    '''
    Now this does work on the genome level
    '''
    mm_level = kwargs.get('mm_level', False)
    stb = s2b

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

# def _backfill_blanks(db, s2l):
#     scaffs = list(set(s2l.keys()) - set(db['scaffold'].unique()))
#     bdb = pd.DataFrame({'scaffold':scaffs, 'length':[s2l[s] for s in scaffs]})
#
#     # make some adjustments
#     bdb['bases_w_0_coverage'] = bdb['length']
#
#     # append
#     db = db.append(bdb, sort=True)
#
#     # fill 0
#     return db.fillna(0)

# def interate_sdb_mm(sdb, on='genome', s2l=None, stb=None):
#     '''
#     For the dataframe, iterate through each mm and a dataframe with ALL scaffolds at that mm level
#     (including blanks)
#     '''
#     for g, db in sdb.groupby(on):
#         if s2l == None:
#             gs2l = db.set_index('scaffold')['length'].to_dict()
#         else:
#             gs2l = {s:s2l[s] for s in [x for x, b in stb.items() if b == g]}
#         mms = sorted(db['mm'].unique())
#         for mm in mms:
#             # get all the ones that you can
#             dd = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset='scaffold', keep='last')
#             #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))
#
#             # backfill with blanks
#             dd = _backfill_blanks(dd, gs2l)
#             #print("mm={0}; len={1}; tl={2}".format(mm, len(dd), len(s2l.keys())))
#
#             yield g, mm, dd

def estimate_breadth(coverage):
    '''
    Estimate breadth based on coverage

    Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000
    '''
    import numpy as np
    return (-1) * np.exp(-1 * ((0.883) * coverage)) + 1



def _genome():
    return 'genome'

def load_scaff2bin(input_stb, IS=None):
    '''
    From the input, load an .stb. The input can be a lot of things, though
    '''

    # Check if there is nothing there
    if IS is not None:
        if (input_stb == []) | (input_stb is None):
            s2l = IS.get('scaffold2length')
            stb = {}
            for scaffold in list(s2l.keys()):
                stb[scaffold] = 'all_scaffolds'
            logging.info('Scaffold to bin will consider all scaffolds the same genome')
            return stb

    # Check if this is a .fasta file
    try:
        stb = gen_stb(input_stb)
        logging.info('Scaffold to bin was made using .fasta files')
        return stb
    except:
        pass

    # Check if this is a regular stb file
    if (len(input_stb) == 1):
        try:
            stb = parse_stb(input_stb[0])
            logging.info('Scaffold to bin was made using .stb file')
            return stb
        except:
            pass

    # Check if you didn't have an input
    if (input_stb == []):
        return {}

    # Fail
    logging.error('Could not load the scaffold to bin file!')
    assert False

def parse_stb(stb_loc):
    stb = {}
    with open(stb_loc, "r") as ins:
        for line in ins:
            linewords = line.strip().split('\t')
            scaffold,b = linewords[:2]
            scaffold = scaffold.strip()
            b = b.strip()
            if scaffold not in stb:
                stb[scaffold] = []
            stb[scaffold] = b
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

def calc_bin2scaffols(stb):
    b2s = {}
    for scaffold, bin in stb.items():
        if bin not in b2s:
            b2s[bin] = set()
        b2s[bin].add(scaffold)
    return b2s

def calc_mms(covT):
    '''
    Return a list of mms from scaffold to mm to cov
    '''
    mms = set()
    for scaff, covt in covT.items():
        mms = mms.union(covt.keys())
    return sorted(list(mms))

def generate_genome_coverage_array(covT, s2l, order=None, maxMM=100, mask_edges=0):
    '''
    Geneate a pandas series of coverage values indexed by genome position.

    Arguments:
        covT = The covT object (IS.get('covT')). In this case should have scaffold as well. (scaffold -> mm -> coverage array)
        order = Order of scaffolds to use. Otherwise will just use the sorted order of s2l keys
        s2l = Scaffold 2 length
        maxMM = The maximum number of mms to make it into the final coverage array
        mask_edges = Remove this much from the edge of each scaffold

    Returns:
        cov = Genome-wide array of coverage values
        scaff2index_addition = Dictionary with the property: scaff2index_addition[scaffold] = scalar to sum to scaffold index to get genome_index
    '''
    arrs = []
    tally = 0
    scaff2index_addition = {}

    if order is None:
        order = sorted(list(s2l.keys()), key=s2l.get, reverse=True)

    for scaff in order:
        slen = s2l[scaff]

        if scaff in covT:
            cov = inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(covT[scaff], maxMM=maxMM, fill_zeros=slen)
        else:
            cov = pd.Series(index=np.arange(slen))


        if mask_edges:
            if len(cov) >= (mask_edges * 2):
                cov = cov[mask_edges:len(cov)-mask_edges]
                slen = slen - (mask_edges * 2)
            else: # really short scaffold
                cov = pd.Series()
                slen = 0

        arrs.append(cov)

        # Set the index addition
        scaff2index_addition[scaff] = tally

        # Update the tally
        tally += slen

    cov = pd.concat(arrs, verify_integrity=False, ignore_index=True).reset_index(drop=True).fillna(0)

    return cov, scaff2index_addition
