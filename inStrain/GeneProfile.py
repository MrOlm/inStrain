#!/usr/bin/env python

import os
import csv
import sys
import glob
import logging
import argparse
import traceback

import Bio
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

def calculate_gene_metrics(IS, Gdb, gene2sequence, **kwargs):
    '''
    Calculate the metrics of all genes on a parallelized scaffold-level basis
    '''
    # get arguments for the wrapper
    p = int(kwargs.get('processes', 6))

    # figure out your overlap level
    scaffolds_with_genes = set(Gdb['scaffold'].unique())
    scaffolds_in_IS = set(IS._get_covt_keys())
    scaffolds_to_profile = scaffolds_with_genes.intersection(scaffolds_in_IS)
    logging.info("{0} scaffolds with genes, {1} in the IS, {2} to compare".format(
            len(scaffolds_with_genes), len(scaffolds_in_IS), len(scaffolds_to_profile)))
    # print("{0} scaffolds with genes, {1} in the IS, {2} to compare".format(
    #         len(scaffolds_with_genes), len(scaffolds_in_IS), len(scaffolds_to_profile)))

    # iterate
    GeneProfiles = []

    # Make a global scaff table
    global CumulativeSNVtable
    CumulativeSNVtable = IS.get('cumulative_snv_table')
    if len(CumulativeSNVtable) > 0:
        CumulativeSNVtable = CumulativeSNVtable.sort_values('mm')
    else:
        CumulativeSNVtable = pd.DataFrame(columns=['scaffold'])

    if p > 1:
        ex = concurrent.futures.ProcessPoolExecutor(max_workers=p)
        total_cmds = len([x for x in iterate_commands(IS, scaffolds_to_profile, Gdb, gene2sequence, kwargs)])
        wait_for = [ex.submit(profile_genes_wrapper, cmd) for cmd in iterate_commands(IS, scaffolds_to_profile, Gdb, gene2sequence, kwargs)]
        for f in tqdm(futures.as_completed(wait_for), total=total_cmds, desc='Running gene-level calculations on scaffolds'):
            try:
                results = f.result()
                GeneProfiles.append(results)
            except:
                logging.error("We had a failure! Not sure where!")

    else:
        for cmd in tqdm(iterate_commands(IS, scaffolds_to_profile, Gdb, gene2sequence, kwargs),
                        desc='Running gene-level calculations on scaffolds',
                        total = len(scaffolds_to_profile)):
            GeneProfiles.append(profile_genes_wrapper(cmd))

    # Return
    name2result = {}
    for i, name in enumerate(['coverage', 'clonality', 'SNP_density', 'SNP_mutation_types']):
        name2result[name] = pd.concat([G[i] for G in GeneProfiles])

    return name2result

def profile_genes(IS, scaffold, gdb, gene2sequence, **kwargs):
    '''
    This is the money that gets multiprocessed

    Relies on having a global "CumulativeSNVtable"

    * Calculate the clonality, coverage, linkage, and SNV_density for each gene
    * Determine whether each SNP is synynomous or nonsynonymous
    '''
    # Calculate gene-level coverage
    covTs = IS.get('covT', scaffolds=[scaffold])
    if scaffold not in covTs:
        logging.info("{0} isnt in covT!".format(scaffold))
        cdb = pd.DataFrame()
    else:
        covT = covTs[scaffold]
        cdb = calc_gene_coverage(gdb, covT)
        del covT

    # Calculate gene-level clonality
    covTs = IS.get('clonT', scaffolds=[scaffold])
    if (covTs is None):
        logging.info("{0} isnt in clovT!".format(scaffold))
        cldb = pd.DataFrame()
    if scaffold not in covTs:
        logging.info("{0} isnt in clovT!".format(scaffold))
        cldb = pd.DataFrame()
    else:
        covT = covTs[scaffold]
        cldb = calc_gene_clonality(gdb, covT)
        del covT

    # I DECIDED THIS DOESN'T MAKE SENSE TO DO

    # # Calculate gene-level linkage
    # Ldb = IS.get('raw_linkage_table')
    # Ldb = Ldb[Ldb['scaffold'] == 'scaffold'].sort_values('mm')
    # if len(Ldb) == 0:
    #     ldb = pd.DataFrame()
    # else:
    #     ldb = calc_gene_linkage(gdb, ldb)

    # Calculate gene-level SNP densisty
    # Ldb = IS.get('cumulative_snv_table')
    # Ldb = Ldb[Ldb['scaffold'] == scaffold].sort_values('mm')
    # THIS IS A GLOBAL NOW

    Ldb = CumulativeSNVtable[CumulativeSNVtable['scaffold'] == scaffold]

    if len(Ldb) == 0:
        ldb = pd.DataFrame()
    else:
        ldb = calc_gene_snp_density(gdb, Ldb)

    # Determine whether SNPs are synonmous or non-synonmous
    # Sdb = IS.get_nonredundant_snv_table()
    scdb = Ldb
    if (scdb is None) or (len(scdb) == 0):
        Sdb = pd.DataFrame()
        sdb = pd.DataFrame()
    else:
        Sdb = scdb.drop_duplicates(subset=['scaffold', 'position'], keep='last')\
                    .sort_index().drop(columns=['mm'])

        Sdb = Sdb[Sdb['scaffold'] == scaffold]
        Sdb = Sdb[Sdb['cryptic'] == False]

        if 'morphia' in Sdb.columns:
            Sdb['morphia'] = Sdb['morphia'].astype(int)
            # Sdb = Sdb[Sdb['morphia'] == 2]
            Sdb = Sdb[(Sdb['morphia'] > 0) & (Sdb['morphia'] <= 2)]
            #Sdb = Sdb.drop(columns="morphia")

        elif 'allele_count' in Sdb.columns:
            Sdb['allele_count'] = Sdb['allele_count'].astype(int)
            # Sdb = Sdb[Sdb['allele_count'] == 2]
            Sdb = Sdb[(Sdb['allele_count'] > 0) & (Sdb['allele_count'] <= 2)]
            #Sdb = Sdb.drop(columns="allele_count")

        Sdb = Sdb.drop(columns="cryptic")
        Sdb['position'] = Sdb['position'].astype(int)

        if len(Sdb) == 0:
            sdb = pd.DataFrame()
        else:
            sdb = characterize_SNPs(gdb, Sdb, gene2sequence)

    assert len(Sdb) == len(sdb)
    if len(Sdb) > 0:
        sdb = pd.merge(Sdb, sdb, on=['position'], how='left').reset_index(drop=True)

    results = (cdb, cldb, ldb, sdb)

    return results

def calc_gene_coverage(gdb, covT):
    '''
    Gene-level and mm-level coverage
    '''
    table = defaultdict(list)

    # for mm in sorted(list(covT.keys())):
    #     cov = inStrain.profileUtilities._mm_counts_to_counts_shrunk(covT, mm)
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
    # p2c = {}
    # mms = sorted([int(mm) for mm in list(clonT.keys())])
    # for mm in mms:
    #     for pos, val in clonT[mm].items():
    #         p2c[pos] = val
    #
    #     inds = []
    #     vals = []
    #     for ind in sorted(p2c.keys()):
    #         inds.append(ind)
    #         vals.append(p2c[ind])
    #
    #     yield mm, pd.Series(data = vals, index = np.array(inds).astype('int'))

def calc_gene_clonality(gdb, clonT):
    '''
    Gene-level and mm-level clonality
    '''
    table = defaultdict(list)

    for mm, cov in iterate_clonT_mms(clonT):
    # for mm in sorted(list(clonT.keys())):
    #     cov = inStrain.profileUtilities._get_basewise_clons3(clonT, mm)
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

            table['clonality'].append(gcov.mean())
            table['microdiversity'].append(microdiversity)
            table['masked_breadth'].append(len(gcov) / gLen)
            table['mm'].append(mm)

    return pd.DataFrame(table)

# def calc_gene_linkage(gdb, ldb, on='r2'):
#     '''
#     Gene-level and mm-level clonality
#     '''
#     table = defaultdict(list)
#
#     for mm in sorted(ldb['mm'].unique()):
#         db = ldb[ldb['mm'] <= mm].drop_duplicates(subset=['scaffold', 'position_A', 'position_B'], keep='last')
#         cov = db.set_index('position_A')[on].sort_index()
#         if len(cov) == 0:
#             continue
#
#         for i, row in gdb.iterrows():
#             gcov = cov.loc[int(row['start']):int(row['end'])]
#             gLen = abs(row['end'] - row['start']) + 1
#
#             table['gene'].append(row['gene'])
#             table['linkage'].append(gcov.mean())
#             table['mm'].append(mm)
#
#     return pd.DataFrame(table)

def calc_gene_snp_density(gdb, ldb):
    '''
    Gene-level and mm-level clonality
    '''
    table = defaultdict(list)

    for mm in sorted(ldb['mm'].unique()):
        db = ldb[ldb['mm'] <= mm].drop_duplicates(subset=['scaffold', 'position'], keep='last')
        cov = db.set_index('position')['refBase'].sort_index()
        if len(cov) == 0:
            continue

        for i, row in gdb.iterrows():
            gcov = cov.loc[int(row['start']):int(row['end'])]
            gLen = abs(row['end'] - row['start']) + 1

            table['gene'].append(row['gene'])
            table['SNPs_per_bp'].append(len(gcov) / gLen)
            table['mm'].append(mm)

    return pd.DataFrame(table)

def characterize_SNPs(gdb, Sdb, gene2sequence):
    '''
    Determine the type of SNP (synonymous, non-synynomous, or intergenic)
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
            original_sequence = gene2sequence[db['gene'].tolist()[0]]
            if db['direction'].tolist()[0] == '-1':
                original_sequence = original_sequence.reverse_complement()

            # Make the new sequence
            snp_start = row['position'] - db['start'].tolist()[0]
            new_sequence = original_sequence.tomutable()
            new_sequence[snp_start] = row['varBase']
            if new_sequence[snp_start] == original_sequence[snp_start]:
                new_sequence[snp_start] = row['conBase']
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

def iterate_commands(IS, scaffolds_to_profile, Gdb, gene2sequence, kwargs):
    '''
    Break into individual scaffolds
    '''
    for scaffold, gdb in Gdb.groupby('scaffold'):
        if scaffold not in scaffolds_to_profile:
            continue

        cmd = Command()
        cmd.IS = IS
        cmd.scaffold = scaffold
        cmd.gdb = gdb
        cmd.arguments = kwargs
        cmd.gene2sequence = {g:gene2sequence[g] for g in gdb['gene'].tolist()}

        yield cmd

def profile_genes_wrapper(cmd):
    '''
    Take a command and profile the scaffold
    '''
    logging.debug('running {0}'.format(cmd.scaffold))
    try:
        return profile_genes(cmd.IS, cmd.scaffold, cmd.gdb, cmd.gene2sequence, **cmd.arguments)
    except Exception as e:
        print(e)
        traceback.print_exc()
        logging.error("whole scaffold exception- {0}".format(str(cmd.scaffold)))
        return pd.DataFrame({'Failed':[True]})

class Command():
    def __init__(self):
        pass

def parse_genes(gene_file, **kwargs):
    if gene_file[-4:] == '.fna':
        return parse_prodigal_genes(gene_file)
    elif ((gene_file[-3:] == '.gb') | (gene_file[-4:] == '.gbk')):
        return parse_genbank_genes(gene_file)
    else:
        print("I dont know how to process {0}".format(gene_file))
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
        table['partial'].append('partial=00' not in record.description)

        # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE
        table['start'].append(int(record.description.split("#")[1].strip())-1)
        table['end'].append(int(record.description.split("#")[2].strip())-1)

        gene2sequence[gene] = record.seq

    Gdb = pd.DataFrame(table)
    logging.info("{0}% of the input {1} genes were marked as incomplete".format((len(Gdb[Gdb['partial'] == True])/len(Gdb))*100, len(Gdb)))

    return Gdb, gene2sequence

def parse_genbank_genes(gene_file, gene_name='gene'):
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
                table['end'].append(loc.end)

                gene2sequence[gene] = record.seq

    Gdb = pd.DataFrame(table)
    logging.info("{0}% of the input {1} genes were marked as compound".format((len(Gdb[Gdb['partial'] != False])/len(Gdb))*100, len(Gdb)))

    return Gdb, gene2sequence

def get_gene_info(IS,  ANI_level=0):
    #IS = inStrain.SNVprofile.SNVprofile(IS_loc)

     # Get the mm level
    mm = _get_mm(IS, ANI_level)

    # Load all genes
    Gdb = IS.get('genes_table')

    # Get the coverage
    db = IS.get('genes_coverage')
    db = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['gene'], keep='last')
    del db['mm']
    Gdb = pd.merge(Gdb, db, on='gene', how='left')

    # Get the clonality
    db = IS.get('genes_clonality')
    db = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['gene'], keep='last')
    del db['mm']
    Gdb = pd.merge(Gdb, db, on='gene', how='left')

    # Get the SNP density
    db = IS.get('genes_SNP_density')
    if len(db) > 0:
        db = db[db['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['gene'], keep='last')
        del db['mm']
        Gdb = pd.merge(Gdb, db, on='gene', how='left')
    else:
        logging.info('You have no SNPs! Skipping SNP density calculation')

    Gdb['min_ANI'] = ANI_level

    return Gdb

def _get_mm(IS, ANI):
    '''
    Get the mm corresponding to an ANI level in an IS
    '''
    if ANI > 1:
        ANI = ANI / 100

    rLen = IS.get('read_report')['mean_pair_length'].tolist()[0]
    mm = int(round((rLen - (rLen * ANI))))
    return mm

class Controller():

    def main(self, args):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        #args = self.parse_arguments(sys_args)
        args = self.validate_input(args)

        vargs = vars(args)
        IS = vargs.pop('IS')
        GF = vargs.pop('gene_file')

        # Read the genes file
        logging.debug('Loading genes')
        Gdb, gene2sequence = parse_genes(GF, **vargs)

        # Calculate all your parallelized gene-level stuff
        name2result = calculate_gene_metrics(IS, Gdb, gene2sequence, **vargs)

        # Store information
        IS.store('genes_fileloc', GF, 'value', 'Location of genes .faa file that was used to call genes')
        IS.store('genes_table', Gdb, 'pandas', 'Location of genes in the associated genes_file')
        IS.store('genes_coverage', name2result['coverage'], 'pandas', 'Coverage of individual genes')
        IS.store('genes_clonality', name2result['clonality'], 'pandas', 'Clonality of individual genes')
        IS.store('genes_SNP_density', name2result['SNP_density'], 'pandas', 'SNP density of individual genes')
        IS.store('SNP_mutation_types', name2result['SNP_mutation_types'], 'pandas', 'The mutation types of SNPs')

        # Store the output
        out_base = IS.get_location('output') + os.path.basename(IS.get('location')) + '_'
        Gdb = get_gene_info(IS)
        Gdb.to_csv(out_base + 'gene_info.tsv', index=False, sep='\t')

        try:
            name2result['SNP_mutation_types'].to_csv(out_base + 'SNP_mutation_types.tsv', index=False, sep='\t')
        except:
            pass

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

    def parse_arguments(self, args):
        '''
        Argument parsing
        '''
        parser = argparse.ArgumentParser(description= """
            A script that runs gene-based analyses on inStrain SNV profiles\n

            Required input: a prodigal .fna gene calls file (created with prodigal -d). Going to be opinionated about this and not support alternatives.
            This file is easy to recreate from other / custom data - please look under examples to see how it should be formatted.
            """, formatter_class=argparse.RawTextHelpFormatter)

        # Required positional arguments
        parser.add_argument("-i", '--IS', help="an inStrain object", required=True)
        parser.add_argument("-g", "--gene_file", action="store", required=True, \
            help='Path to prodigal .fna genes file.')
        parser.add_argument("-p", "--processes", action="store", default=6, type=int, \
                            help='Threads to use for multiprocessing')

        # Parse
        if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
            parser.print_help()
            sys.exit(0)
        else:
            return parser.parse_args(args)

# if __name__ == '__main__':
#     Controller().main(sys.argv[1:])
