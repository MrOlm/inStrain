#!/usr/bin/env python

import os
import csv
import sys
import time
import glob
import logging
import warnings
import argparse
import traceback
import multiprocessing

import Bio
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import concurrent.futures
from concurrent import futures
from inStrain import SNVprofile
from collections import defaultdict

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.codonalign.codonalphabet import default_codon_table

import inStrain.SNVprofile
import inStrain.controller
#import inStrain.profileUtilities
import inStrain.logUtils

class Controller():
    '''
    The command line access point to the program
    '''
    def main(self, args):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        args = self.validate_input(args)

        vargs = vars(args)
        IS = vargs.pop('IS')
        GF = vargs.pop('gene_file')

        # Read the genes file
        logging.debug('Loading genes')
        GdbP, gene2sequence = parse_genes(GF, **vargs)

        # Calculate all your parallelized gene-level stuff
        name2result = calculate_gene_metrics(IS, GdbP, gene2sequence, **vargs)

        # Store information
        IS.store('genes_fileloc', GF, 'value', 'Location of genes file that was used to call genes')
        IS.store('genes_table', GdbP, 'pandas', 'Location of genes in the associated genes_file')
        IS.store('genes_coverage', name2result['coverage'], 'pandas', 'Coverage of individual genes')
        IS.store('genes_clonality', name2result['clonality'], 'pandas', 'Clonality of individual genes')
        IS.store('genes_SNP_count', name2result['SNP_density'], 'pandas', 'SNP density and counts of individual genes')
        IS.store('SNP_mutation_types', name2result['SNP_mutation_types'], 'pandas', 'The mutation types of SNPs')

        if vargs.get('store_everything', False):
            IS.store('gene2sequence', gene2sequence, 'pickle', 'Dicitonary of gene -> nucleotide sequence')

        # Store the output
        IS.generate('gene_info')
        IS.generate("SNVs")

    def validate_input(self, args):
        '''
        Validate and mess with the arguments a bit
        '''
        # Make sure the IS object is OK and load it
        assert os.path.exists(args.IS)
        args.IS = inStrain.SNVprofile.SNVprofile(args.IS)

        # Set up the logger
        log_loc = args.IS.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)

        return args

def gene_profile_worker(gene_cmd_queue, gene_result_queue, single_thread=False):
    '''
    Worker to profile splits
    '''
    while True:
        # Get command
        if not single_thread:
            cmds = gene_cmd_queue.get(True)
        else:
            try:
                cmds = gene_cmd_queue.get(timeout=5)
            except:
                return

        # Process cmd
        GPs = profile_genes_wrapper(cmds)
        gene_result_queue.put(GPs)

def profile_genes_wrapper(cmds):
    '''
    Take a group of commands and run geneprofile
    '''
    results = []
    for cmd in cmds:
        try:
            results.append(profile_genes(cmd.scaffold, **cmd.arguments))
        except Exception as e:
            print(e)
            traceback.print_exc()
            logging.error("FAILURE GeneException {0}".format(str(cmd.scaffold)))
            results.append(None)
    return results

def calculate_gene_metrics(IS, GdbP, gene2sequenceP, **kwargs):
    '''
    Calculate the metrics of all genes on a parallelized scaffold-level basis

    IS = Initialized inStrain.SNVprofile
    GdbP = List of gene locations
    gene2sequenceP = Dicitonary of gene -> nucleotide sequence
    '''
    inStrain.logUtils.log_checkpoint("GeneProfile", "calculate_gene_metrics", "start")

    # Get key word arguments for the wrapper
    p = int(kwargs.get('processes', 6))

    # Make a list of scaffolds to profile the genes of
    scaffolds_with_genes = set(GdbP['scaffold'].unique())
    scaffolds_in_IS = set(IS._get_covt_keys())
    scaffolds_to_profile = scaffolds_with_genes.intersection(scaffolds_in_IS)
    logging.info("{0} scaffolds with genes in the input; {1} scaffolds in the IS, {2} to compare".format(
            len(scaffolds_with_genes), len(scaffolds_in_IS), len(scaffolds_to_profile)))

    # Calculate scaffold -> number of genes to profile
    s2g = GdbP['scaffold'].value_counts().to_dict()
    kwargs['s2g'] = s2g

    # Make global objects for the profiling
    inStrain.logUtils.log_checkpoint("GeneProfile", "make_globals", "start")
    global CumulativeSNVtable
    CumulativeSNVtable = IS.get('cumulative_snv_table')
    if len(CumulativeSNVtable) > 0:
        CumulativeSNVtable = CumulativeSNVtable.sort_values('mm')
    else:
        CumulativeSNVtable = pd.DataFrame(columns=['scaffold'])

    global covTs
    covTs = IS.get('covT', scaffolds=scaffolds_to_profile)

    global clonTs
    clonTs = IS.get('clonT', scaffolds=scaffolds_to_profile)

    global gene2sequence
    gene2sequence =  gene2sequenceP

    global Gdb
    Gdb = GdbP
    inStrain.logUtils.log_checkpoint("GeneProfile", "make_globals", "end")

    # Generate commands and queue them
    logging.debug('Creating commands')
    cmd_groups = [x for x in iterate_commands(scaffolds_to_profile, Gdb, kwargs)]
    logging.debug('There are {0} cmd groups'.format(len(cmd_groups)))

    inStrain.logUtils.log_checkpoint("GeneProfile", "create_queue", "start")
    gene_cmd_queue = multiprocessing.Queue()
    gene_result_queue = multiprocessing.Queue()
    GeneProfiles = []

    for cmd_group in cmd_groups:
        gene_cmd_queue.put(cmd_group)

    inStrain.logUtils.log_checkpoint("GeneProfile", "create_queue", "end")

    if p > 1:
        logging.debug('Establishing processes')
        processes = []
        for i in range(0, p):
            processes.append(multiprocessing.Process(target=gene_profile_worker, args=(gene_cmd_queue, gene_result_queue)))
        for proc in processes:
            proc.start()

        # Set up progress bar
        pbar = tqdm(desc='Profiling genes: ', total=len(cmd_groups))

        # Get the results
        recieved_profiles = 0
        while recieved_profiles < len(cmd_groups):
            GPs = gene_result_queue.get()
            recieved_profiles += 1
            pbar.update(1)
            for GP in GPs:
                if GP is not None:
                    logging.debug(GP[4])
                    GeneProfiles.append(GP)

        # Close multi-processing
        for proc in processes:
            proc.terminate()

        # Close progress bar
        pbar.close()

    else:
        gene_profile_worker(gene_cmd_queue, gene_result_queue, single_thread=True)
        logging.info("Done profiling genes")

        # Get the genes
        recieved_profiles = 0
        while recieved_profiles < len(cmd_groups):
            logging.debug('going to grab at {0}'.format(recieved_profiles))
            GPs = gene_result_queue.get(timeout=5)
            logging.debug('did a grab at {0}'.format(recieved_profiles))
            recieved_profiles += 1
            for GP in GPs:
                if GP is not None:
                    logging.debug(GP[4])
                    GeneProfiles.append(GP)

    inStrain.logUtils.log_checkpoint("GeneProfile", "return_results", "start")
    name2result = {}
    for i, name in enumerate(['coverage', 'clonality', 'SNP_density', 'SNP_mutation_types']):
        name2result[name] = pd.concat([G[i] for G in GeneProfiles])
    inStrain.logUtils.log_checkpoint("GeneProfile", "return_results", "end")

    inStrain.logUtils.log_checkpoint("GeneProfile", "calculate_gene_metrics", "end")
    return name2result

def profile_genes(scaffold, **kwargs):
    '''
    This is the money that gets multiprocessed

    Relies on having a global "Gdb", "gene2sequence", "CumulativeSNVtable", "covTs", and "clonTs"

    * Calculate the clonality, coverage, linkage, and SNV_density for each gene
    * Determine whether each SNP is synynomous or nonsynonymous
    '''
    # Log
    pid = os.getpid()
    log_message = "\nSpecialPoint_genes {0} PID {1} whole start {2}".format(scaffold, pid, time.time())

    # For testing purposes
    if ((scaffold == 'FailureScaffoldHeaderTesting')):
        assert False

    # Get the list of genes for this scaffold
    gdb = Gdb[Gdb['scaffold'] == scaffold]

    # Calculate gene-level coverage
    log_message += "\nSpecialPoint_genes {0} PID {1} coverage start {2}".format(scaffold, pid, time.time())
    if scaffold not in covTs:
        logging.info("{0} isnt in covT!".format(scaffold))
        cdb = pd.DataFrame()
    else:
        covT = covTs[scaffold]
        cdb = calc_gene_coverage(gdb, covT)
        del covT
    log_message += "\nSpecialPoint_genes {0} PID {1} coverage end {2}".format(scaffold, pid, time.time())

    # Calculate gene-level clonality
    log_message += "\nSpecialPoint_genes {0} PID {1} clonality start {2}".format(scaffold, pid, time.time())
    if scaffold not in clonTs:
        logging.info("{0} isnt in clovT!".format(scaffold))
        cldb = pd.DataFrame()
    else:
        clonT = clonTs[scaffold]
        cldb = calc_gene_clonality(gdb, clonT)
        del clonT
    log_message += "\nSpecialPoint_genes {0} PID {1} clonality end {2}".format(scaffold, pid, time.time())

    # Determine whether SNPs are synonmous or non-synonmous
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_character start {2}".format(scaffold, pid, time.time())
    Ldb = CumulativeSNVtable[CumulativeSNVtable['scaffold'] == scaffold]
    if len(Ldb) == 0:
        sdb = pd.DataFrame()
    else:
        sdb = Characterize_SNPs_wrapper(Ldb, gdb, gene2sequence)
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_character end {2}".format(scaffold, pid, time.time())

    # Calculate gene-level SNP counts
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts start {2}".format(scaffold, pid, time.time())
    if len(Ldb) == 0:
        ldb = pd.DataFrame()
        sublog = ''
    else:
        #ldb = calc_gene_snp_density(gdb, Ldb)
        ldb, sublog = calc_gene_snp_counts(gdb, Ldb, sdb, gene2sequence, scaffold=scaffold)
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts end {2}".format(scaffold, pid, time.time())
    log_message += sublog

    log_message += "\nSpecialPoint_genes {0} PID {1} whole end {2}".format(scaffold, pid, time.time())

    results = (cdb, cldb, ldb, sdb, log_message)

    return results

def profile_genes_from_profile(scaffold, gdb, covT, clonT, Ldb, gene2sequence):
    '''
    Call profile genes from elsewhere

    Arguments:
        scaffold = name of scaffold
        gdb = gene_datatable
        covT = covT for this scaffold
        clonT = clonT for this scaffold
        Ldb = cumulative SNP table for this scaffold

    * Calculate the clonality, coverage, linkage, and SNV_density for each gene
    * Determine whether each SNP is synynomous or nonsynonymous
    '''
    # Log
    log = inStrain.logUtils.get_worker_log('ProfileGenes', scaffold, 'start')

    # For testing purposes
    if ((scaffold == 'FailureScaffoldHeaderTesting')):
        assert False

    # Calculate gene-level coverage
    #log_message += "\nSpecialPoint_genes {0} PID {1} coverage start {2}".format(scaffold, pid, time.time())
    cdb = calc_gene_coverage(gdb, covT)
    #log_message += "\nSpecialPoint_genes {0} PID {1} coverage end {2}".format(scaffold, pid, time.time())

    # Calculate gene-level clonality
    # log_message += "\nSpecialPoint_genes {0} PID {1} clonality start {2}".format(scaffold, pid, time.time())
    cldb = calc_gene_clonality(gdb, clonT)
    # log_message += "\nSpecialPoint_genes {0} PID {1} clonality end {2}".format(scaffold, pid, time.time())

    # Determine whether SNPs are synonmous or non-synonmous
    #log_message += "\nSpecialPoint_genes {0} PID {1} SNP_character start {2}".format(scaffold, pid, time.time())
    sdb = Characterize_SNPs_wrapper(Ldb, gdb, gene2sequence)
    #log_message += "\nSpecialPoint_genes {0} PID {1} SNP_character end {2}".format(scaffold, pid, time.time())

    # Calculate gene-level SNP counts
    #log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts start {2}".format(scaffold, pid, time.time())
    ldb, sublog = calc_gene_snp_counts(gdb, Ldb, sdb, gene2sequence, scaffold=scaffold)
    #log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts end {2}".format(scaffold, pid, time.time())
    #log_message += sublog

    log += inStrain.logUtils.get_worker_log('ProfileGenes', scaffold, 'end')

    results = (cdb, cldb, ldb, sdb, log)

    return results

def calc_gene_coverage(gdb, covT):
    '''
    Gene-level and mm-level coverage
    '''
    table = defaultdict(list)

    for mm, cov in iterate_covT_mms(covT):
        if len(cov) == 0:
            continue

        for i, row in gdb.iterrows():
            gcov = cov.loc[int(row['start']):int(row['end'])]
            gLen = abs(row['end'] - row['start']) + 1

            table['gene'].append(row['gene'])
            table['coverage'].append(gcov.sum() / gLen)
            table['breadth'].append(len(gcov) / gLen)
            table['mm'].append(mm)

    return pd.DataFrame(table)

def iterate_clonT_mms(clonT):
    p2c = {}
    mms = sorted([int(mm) for mm in list(clonT.keys())])
    for mm in mms:
        for pos, val in clonT[mm].items():
            p2c[pos] = val

        inds = []
        vals = []
        for ind in sorted(p2c.keys()):
            inds.append(ind)
            vals.append(p2c[ind])

        yield mm, pd.Series(data = vals, index = np.array(inds).astype('int'))

def iterate_covT_mms(clonT):
    counts = pd.Series()
    mms = sorted([int(mm) for mm in list(clonT.keys())])
    for mm in mms:
        count = clonT[mm]
        counts = counts.add(count, fill_value=0)
        yield mm, counts

def calc_gene_clonality(gdb, clonT):
    '''
    Gene-level and mm-level clonality
    '''
    table = defaultdict(list)

    for mm, cov in iterate_clonT_mms(clonT):
        if len(cov) == 0:
            continue

        for i, row in gdb.iterrows():
            gcov = cov.loc[int(row['start']):int(row['end'])]
            gLen = abs(row['end'] - row['start']) + 1

            table['gene'].append(row['gene'])

            try:
                microdiversity = 1 - gcov.mean()
            except :
                microdiversity = np.nan

            #table['clonality'].append(gcov.mean())
            table['nucl_diversity'].append(microdiversity)
            table['breadth_minCov'].append(len(gcov) / gLen)
            table['mm'].append(mm)

    return pd.DataFrame(table)
#
# def calc_gene_snp_density(gdb, ldb):
#     '''
#     Gene-level and mm-level clonality
#     '''
#     table = defaultdict(list)
#
#     for mm in sorted(ldb['mm'].unique()):
#         db = ldb[ldb['mm'] <= mm].drop_duplicates(subset=['scaffold', 'position'], keep='last')
#         cov = db.set_index('position')['ref_base'].sort_index()
#         if len(cov) == 0:
#             continue
#
#         for i, row in gdb.iterrows():
#             gcov = cov.loc[int(row['start']):int(row['end'])]
#             gLen = abs(row['end'] - row['start']) + 1
#
#             table['gene'].append(row['gene'])
#             table['SNPs_per_bp'].append(len(gcov) / gLen)
#             table['mm'].append(mm)
#
#     return pd.DataFrame(table)



def count_sites(seq, k=1, codon_table=None):
    '''
    From a nucleotide sequence and codon table, calculate S and N sites
    '''
    codon_lst = convert_to_codons(seq)

    if codon_table is None:
        codon_table = default_codon_table
    S_site = 0.0  # synonymous sites
    N_site = 0.0  # non-synonymous sites
    purine = ('A', 'G')
    pyrimidine = ('T', 'C')
    base_tuple = ('A', 'T', 'C', 'G')
    for codon in codon_lst:
        neighbor_codon = {'transition': [], 'transversion': []}
        # classify neighbor codons
        codon = codon.replace('U', 'T')
        if codon == '---':
            continue
        if 'N' in codon:
            continue
        for n, i in enumerate(codon):
            for j in base_tuple:
                if i == j:
                    pass
                elif i in purine and j in purine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                elif i in pyrimidine and j in pyrimidine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                else:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transversion'].append(this_codon)
        # count synonymous and non-synonymous sites
        #codon = codon.replace('T', 'U')
        if (codon in ['TAG', 'TAA', 'TGA']):
            #print("STOP DETECTED")
            continue
        aa = codon_table.forward_table[codon]
        this_codon_N_site = this_codon_S_site = 0
        for neighbor in neighbor_codon['transition']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += 1
            else:
                this_codon_N_site += 1
        for neighbor in neighbor_codon['transversion']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += k
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += k
            else:
                this_codon_N_site += k
        norm_const = (this_codon_N_site + this_codon_S_site)/3
        S_site += float(this_codon_S_site) / float(norm_const)
        N_site += float(this_codon_N_site) / float(norm_const)
    return (S_site, N_site)

def convert_to_codons(seq):
    codons = []
    for c in zip(*(iter(seq),) * 3):
        co = ''.join(c)
        assert len(co) == 3
        codons.append(co)
    return codons

def calc_gene_snp_counts(gdb, ldb, sdb, gene2sequence, scaffold=None):
    '''
    Count the number of SNPs in each gene, as well as N and S sites

    RELIES ON HAVING gene2sequence AS A GLOBAL (needed for multiprocessing speed)

    Argumnets:
        gdb = table of genes
        ldb = Raw cumulative snp table for a single scaffold (mm-level)
        sdb = SNP table with N and S and I annotated
    '''
    if len(ldb) == 0:
        return pd.DataFrame(), ''

    pid = os.getpid()

    # Merge ldb and sdb
    xdb = pd.merge(ldb, sdb[['position', 'mutation_type', 'gene']],
            on=['position'], how='left').reset_index(drop=True)

    # Calculate counts of N and S sites
    log_message = "\nSpecialPoint_genes {0} PID {1} SNP_counts_SiteCalc start {2}".format(scaffold, pid, time.time())
    table = defaultdict(list)
    for i, row in gdb.iterrows():
        try:
            S_site, N_site = count_sites(gene2sequence[row['gene']])
        except:
            S_site = np.nan
            N_site = np.nan

        table['gene'].append(row['gene'])
        table['S_sites'].append(S_site)
        table['N_sites'].append(N_site)
    SiteDb = pd.DataFrame(table)
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts_SiteCalc end {2}".format(scaffold, pid, time.time())

    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts_geneCalc start {2}".format(scaffold, pid, time.time())
    table = defaultdict(list)
    for mm in sorted(xdb['mm'].unique()):
        # Filter to this mm level and set up for quick indexing
        fdb = xdb[xdb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['scaffold', 'position'], keep='last').sort_values('position').set_index('position')

        for i, row in gdb.iterrows():
            # Calculate gene length
            gLen = abs(row['end'] - row['start']) + 1

            # Subset to this gene
            db = fdb.loc[int(row['start']):int(row['end'])]

            # Report summary stuff
            table['mm'].append(mm)
            table['gene'].append(row['gene'])
            table['gene_length'].append(gLen)
            table['divergent_site_count'].append(len(db))

            # Report type counts
            for allele_count, name in zip([1, 2], ['SNS', 'SNV']):
                table['{0}_count'.format(name)].append(len(db[db['allele_count'] == allele_count]))

                for snp_type in ['N', 'S']:
                    table["{0}_{1}_count".format(name, snp_type)].append(
                    len(db[(db['allele_count'] == allele_count) & (db['mutation_type'] == snp_type)]))

    GGdb = pd.DataFrame(table).merge(SiteDb, on='gene', how='left').reset_index(drop=True)
    log_message += "\nSpecialPoint_genes {0} PID {1} SNP_counts_geneCalc end {2}".format(scaffold, pid, time.time())

    # Calculate dn/ds
    GGdb['dNdS_substitutions'] = [((nC/nS) / (sC/sS)) if ((sC > 0) & (sS > 0)) else np.nan for nC, nS, sC, sS in zip(
                                GGdb['SNS_N_count'], GGdb['N_sites'],
                                GGdb['SNS_S_count'], GGdb['S_sites'])]
    GGdb['pNpS_variants'] = [((nC/nS) / (sC/sS)) if ((sC > 0) & (sS > 0)) else np.nan for nC, nS, sC, sS in zip(
                                    GGdb['SNV_N_count'], GGdb['N_sites'],
                                    GGdb['SNV_S_count'], GGdb['S_sites'])]
    # GGdb['SNPs_per_bp'] = [x/y if y > 0 else np.nan for x, y in \
    #                     zip(GGdb['divergent_site_count'], GGdb['gene_length'])]

    return GGdb, log_message

def Characterize_SNPs_wrapper(Ldb, gdb, gene2sequence):
    '''
    A wrapper for characterizing SNPs

    RELIES ON HAVING gene2sequence AS A GLOBAL (needed for multiprocessing speed)

    Arguments:
        Ldb = CumulativeSNVtable for a single scaffold
        gdb = table of genes

    Returns:
        Sdb = The Cumulative SNV table with extra information added
    '''
    if len(Ldb) == 0:
        return pd.DataFrame()

    # Get a non-nonredundant list of SNPs
    Sdb = Ldb.drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                .sort_index().drop(columns=['mm'])
    Sdb['position'] = Sdb['position'].astype(int)

    # Filter out SNPs that shouldn't be profiled like this
    Sdb = Sdb[Sdb['cryptic'] == False]
    Sdb = Sdb.drop(columns="cryptic")

    if 'morphia' in Sdb.columns:
        col = 'morphia'
    else:
        col = 'allele_count'
    Sdb[col] = Sdb[col].astype(int)
    Sdb = Sdb[(Sdb[col] > 0) & (Sdb[col] <= 2)]

    # Make sure some SNPs to profile remain
    if len(Sdb) == 0:
        return pd.DataFrame()

    # Characterize
    sdb = characterize_SNPs(gdb, Sdb, gene2sequence)
    assert len(Sdb) == len(sdb)
    sdb = pd.merge(Sdb, sdb, on=['position'], how='left').reset_index(drop=True)

    # Return
    return sdb

def characterize_SNPs(gdb, Sdb, gene2sequence):
    '''
    Determine the type of SNP (synonymous, non-synynomous, or intergenic)

    RELIES ON HAVING gene2sequence AS A GLOBAL (needed for multiprocessing speed)
    '''
    table = defaultdict(list)
    for i, row in Sdb.iterrows():
        db = gdb[(gdb['start'] <= row['position']) & (gdb['end'] >= row['position'])]
        if len(db) == 0:
            table['position'].append(row['position'])
            table['mutation_type'].append('I')
            table['mutation'].append('')
            table['gene'].append('')
        elif len(db) > 1:
            table['position'].append(row['position'])
            table['mutation_type'].append('M')
            table['mutation'].append('')
            table['gene'].append(','.join(db['gene'].tolist()))
        else:
            # Get the original sequence
            original_sequence = gene2sequence[db ['gene'].tolist()[0]]
            if db['direction'].tolist()[0] == '-1':
                original_sequence = original_sequence.reverse_complement()

            # Make the new sequence
            snp_start = row['position'] - db['start'].tolist()[0]
            new_sequence = original_sequence.tomutable()
            new_sequence[snp_start] = row['var_base']
            if new_sequence[snp_start] == original_sequence[snp_start]:
                new_sequence[snp_start] = row['con_base']
            new_sequence = new_sequence.toseq()

            # Translate
            if db['direction'].tolist()[0] == '-1':
                old_aa_sequence = original_sequence.reverse_complement().translate()
                new_aa_sequence = new_sequence.reverse_complement().translate()
            else:
                old_aa_sequence = original_sequence.translate()
                new_aa_sequence = new_sequence.translate()

            # old_aa_sequence = original_sequence.translate()
            # new_aa_sequence = new_sequence.translate()

            # Find mutation
            mut_type = 'S'
            mut = 'S:' + str(snp_start)
            for aa in range(0, len(old_aa_sequence)):
                if new_aa_sequence[aa] != old_aa_sequence[aa]:
                    mut_type = 'N'
                    mut = 'N:' + str(old_aa_sequence[aa]) + str(snp_start) + str(new_aa_sequence[aa])
                    break

            # Store
            table['position'].append(row['position'])
            table['mutation_type'].append(mut_type)
            table['mutation'].append(mut)
            table['gene'].append(db['gene'].tolist()[0])

    return pd.DataFrame(table)

def iterate_commands(scaffolds_to_profile, Gdb, kwargs):
    '''
    Break into individual scaffolds
    '''
    processes = kwargs.get('processes', 6)
    s2g = kwargs.get('s2g', None)
    SECONDS = min(60, sum(calc_estimated_runtime(s2g[scaffold]) for scaffold in scaffolds_to_profile)/(processes+1))

    cmds = []
    seconds = 0
    for scaffold, gdb in Gdb.groupby('scaffold'):
        if scaffold not in scaffolds_to_profile:
            continue

        # make this comammand
        cmd = Command()
        cmd.scaffold = scaffold
        cmd.arguments = kwargs

        # Add estimated seconds
        seconds += calc_estimated_runtime(s2g[scaffold])
        cmds.append(cmd)

        # See if you're done
        if seconds >= SECONDS:
            yield cmds
            seconds = 0
            cmds = []

    yield cmds

def calc_estimated_runtime(pairs):
    SLOPE_CONSTANT = 0.01
    return pairs * SLOPE_CONSTANT



class Command():
    def __init__(self):
        pass

def parse_genes(gene_file_loc, **kwargs):
    '''
    Parse a file of genes based on the file extention.

    Currently supported extentions are:
        .fna (prodigal)
        .gb / .gbk (genbank)

    Methods return a table of genes (Gdb) and a dictionary of gene -> sequence
    '''
    if ((gene_file_loc[-4:] == '.fna') | (gene_file_loc[-3:] == '.fa')):
        return parse_prodigal_genes(gene_file_loc)

    elif ((gene_file_loc[-3:] == '.gb') | (gene_file_loc[-4:] == '.gbk')):
        return parse_genbank_genes(gene_file_loc)

    else:
        print("I dont know how to process {0}".format(gene_file_loc))
        raise Exception

def parse_prodigal_genes(gene_fasta):
    '''
    Parse the prodigal .fna file

    Return a datatable with gene info and a dictionary of gene -> sequence
    '''
    table = defaultdict(list)
    gene2sequence = {}
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)

        table['gene'].append(gene)
        table['scaffold'].append("_".join(gene.split("_")[:-1]))
        table['direction'].append(record.description.split("#")[3].strip())
        table['partial'].append('partial=01' in record.description)

        # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE
        table['start'].append(int(record.description.split("#")[1].strip())-1)
        table['end'].append(int(record.description.split("#")[2].strip())-1)

        gene2sequence[gene] = record.seq

    Gdb = pd.DataFrame(table)
    logging.debug("{0:.1f}% of the input {1} genes were marked as incomplete".format((len(Gdb[Gdb['partial'] == True])/len(Gdb))*100, len(Gdb)))

    return Gdb, gene2sequence

def parse_genbank_genes(gene_file, gene_name='gene'):
    '''
    Parse a genbank file. Gets features marked as CDS
    '''
    table = defaultdict(list)
    gene2sequence = {}
    for record in SeqIO.parse(gene_file, 'gb'):
        scaffold = record.id
        for feature in record.features:
            if feature.type == 'CDS':
                gene = feature.qualifiers[gene_name][0]
                loc = feature.location
                if type(loc) is Bio.SeqFeature.CompoundLocation:
                    partial = 'compound'
                else:
                    partial = False

                table['gene'].append(gene)
                table['scaffold'].append(scaffold)
                table['direction'].append(feature.location.strand)
                table['partial'].append(partial)

                table['start'].append(loc.start)
                table['end'].append(loc.end - 1)

                gene2sequence[gene] = feature.location.extract(record).seq

    Gdb = pd.DataFrame(table)
    logging.debug("{0:.1f}% of the input {1} genes were marked as compound".format((len(Gdb[Gdb['partial'] != False])/len(Gdb))*100, len(Gdb)))

    return Gdb, gene2sequence

def get_gene_info(IS,  ANI_level=0):
    #IS = inStrain.SNVprofile.SNVprofile(IS_loc)

     # Get the mm level
    mm = _get_mm(IS, ANI_level)

    # Load all genes
    Gdb = IS.get('genes_table')

    # Load coverage, clonality, and SNPs
    for thing in ['genes_coverage', 'genes_clonality', 'genes_SNP_count']:
        db = IS.get(thing)
        if db is None:
            logging.debug('Skipping {0} gene calculation; you have none'.format(thing))
            continue
        if len(db) == 0:
            logging.debug('Skipping {0} gene calculation; you have none'.format(thing))
            continue
        db = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['gene'], keep='last')
        del db['mm']
        Gdb = pd.merge(Gdb, db, on='gene', how='left')

    Gdb['min_ANI'] = ANI_level

    return Gdb

def _get_mm(IS, ANI):
    '''
    Get the mm corresponding to an ANI level in an IS
    '''
    if ANI > 1:
        ANI = ANI / 100

    rLen = IS.get('mapping_info')['mean_pair_length'].tolist()[0]
    mm = int(round((rLen - (rLen * ANI))))
    return mm
