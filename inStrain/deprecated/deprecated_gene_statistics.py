#!/usr/bin/env python

import csv
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from inStrain import SNVprofile
from collections import defaultdict

def calculate_pi(counts, min_cov = 5):
    '''
    Calculates 1-the probability that two reads have the same allele at a position (per site nucleotide diversity)
    '''
    s = np.sum(counts)
    if s >= min_cov:
        prob = (float(counts[0]) / s) * (float(counts[0]) / s) + (float(counts[1]) / s) * (float(counts[1]) / s) + (float(counts[2]) / s) * (float(counts[2]) / s) + (float(counts[3]) / s) * (float(counts[3]) / s)
        return 1.0-prob

def characterize_snp(gene_index, gene_starts, seqs, gene_direction, snv_table):
    ''' tests if SNPs are synonymous or non-synynomous with respect to a reference
    gene file
    '''

    print("Characterizing SNPs")
    freq = snv_table
    snp_types = {}

    for index,snp in freq.iterrows():
        #get gene for this snp
        if snp['scaffold'] + ":" + str(snp['position']) in gene_index:
            gene = gene_index[snp['scaffold'] + ":" + str(snp['position'])]

            gene_absolute_start = gene_starts[gene]

            #calculate position of SNP within gene (0 based index)
            snp_start = snp['position'] - gene_absolute_start
            original_sequence = seqs[gene]
            if gene_direction[gene] == '-1': #need to do some flipping if we want to have the positions right....
                original_sequence = original_sequence.reverse_complement()

            new_sequence = original_sequence.tomutable()
            new_sequence[snp_start] = snp['varBase']
            if new_sequence[snp_start] == original_sequence[snp_start]:
                new_sequence[snp_start] = snp['conBase']
            new_sequence = new_sequence.toseq()

            if gene_direction[gene] == '-1':
                old_aa_sequence = original_sequence.reverse_complement().translate()
                new_aa_sequence = new_sequence.reverse_complement().translate()
            else:
                old_aa_sequence = original_sequence.translate()
                new_aa_sequence = new_sequence.translate()
            mut_type='S'
            mut = 'S:' + str(snp_start)
            for aa in range(0, len(old_aa_sequence)):
                if new_aa_sequence[aa] != old_aa_sequence[aa]:
                    mut_type = 'N'
                    mut = 'N:' + str(old_aa_sequence[aa]) + str(snp_start) + str(new_aa_sequence[aa])
                    break
        else:
            mut_type = 'I'
            mut = ''
            gene = 'None'

        snp_types[snp['scaffold'] + ":" + str(snp['position'])] = mut_type
        freq.at[index,'mutation_type'] = mut_type
        freq.at[index,'mutation'] = mut
        freq.at[index,'gene'] = gene

    return freq, snp_types

def create_gene_index(gene_fasta):
    gene_index = {}
    gene_starts = {}
    seqs = {}
    gene_direction = {}
    complete_genes = 0.0
    partial_genes = 0.0
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)
        # if 'partial=00' in record.description:
        #     complete_genes += 1
        #     gene_scaf = "_".join(gene.split("_")[:-1])
        #     gene_direction[gene] = record.description.split("#")[3].strip()
        #     # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE
        #     gene_start = int(record.description.split("#")[1].strip())-1
        #     gene_end = int(record.description.split("#")[2].strip())-1
        #     gene_starts[gene] = gene_start
        #     seqs[gene] = record.seq
        #     for i in range(gene_start, gene_end+1):
        #         gene_index[gene_scaf + ":" + str(i)] = gene
        # else:
        #     partial_genes += 1
        complete_genes += 1
        gene_scaf = "_".join(gene.split("_")[:-1])
        gene_direction[gene] = record.description.split("#")[3].strip()
        # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE
        gene_start = int(record.description.split("#")[1].strip())-1
        gene_end = int(record.description.split("#")[2].strip())-1
        gene_starts[gene] = gene_start
        seqs[gene] = record.seq
        for i in range(gene_start, gene_end+1):
            gene_index[gene_scaf + ":" + str(i)] = gene

    #print("Notice: " + str(round(partial_genes *100 / (complete_genes + partial_genes))) + "% of genes were incomplete and snps in these genes were marked I. Pi is not calculated for incomplete genes.")
    return gene_index, gene_starts, gene_direction, seqs


def get_gene_index(x, gene_index):
    try:
        return(gene_index[x])
    except:
        return None

def calculate_gene_pi(gene_index, pi_table, min_cov = 5, min_gene_size = 120):
    '''
    gene_index = a dictionary of scafs and positions to gene IDs calculated by create_gene_index()
    pi_table = pandas dataframe of nucleotide diversities
    min_cov = minimum coverage required for a position to calculate clonality
    min_gene_size = minimum number of positions > min_cov in a gene to calculate clonality for that gene.
    '''

    # Read in clonalities
    print("Calculating gene clonalities")
    # Create gene column
    pi_table['gene'] = pi_table['position'].apply(lambda x: get_gene_index(x, gene_index) )
    # Filter positions to at least min_cov, only positions within genes
    pi_filtered = pi_table[pi_table['coverage'] >= min_cov]
    pi_filtered = pi_filtered[pi_filtered['gene'] != None]
    # calculate means of clonality and coverage per gene
    gene_table = pi_filtered.groupby('gene').agg({'pi':'mean', 'coverage':'mean','position':'count'}).reset_index()
    # remove genes with fewer than min_gene_size positions
    gene_table = gene_table[gene_table['position'] >= min_gene_size]
    # add this sample name to the sample column
    #gene_table['sample'] = sample
    gene_table.rename(columns={'position': 'length_with_coverage'}, inplace=True)

    return gene_table

def calculate_diversity(s, min_cov = 5):
    pi_table = defaultdict(list)
    c = 0
    scaff_list = s.get('scaffold_list')
    counts_table = s.get('counts_table')
    #print(counts_table)
    for scaf in counts_table:
        i = 0
        for counts in scaf:
            # print(counts)
            # break
            pi_table['position'].append(scaff_list[c] + ":" + str(i))
            pi_table['pi'].append(calculate_pi(counts))
            pi_table['coverage'].append(np.sum(counts))
            i += 1
        c += 1

    pi_table = pd.DataFrame(pi_table)
    return(pi_table)

def main(args):
    gene_index, gene_starts, gene_direction, seqs = create_gene_index(args.gene_file)

    # load inStrain object
    print("loading " + args.input)
    print(args.input)
    s = SNVprofile.SNVprofile(args.input)
    print(s.get('location'))
    print("calculating nucleotide diversity in " + args.input)
    pi_table = calculate_diversity(s)
    g = calculate_gene_pi(gene_index, pi_table, min_cov = 5, min_gene_size = 100)

    # Store genes
    s.store('genes', g, 'pandas', 'Gene calls and clonality')
    out_base = s.get_location('output')
    g.to_csv(out_base + 'genes.tsv', index=False, sep='\t')

    print("Determining function of SNPs (N vs S vs O)")
    snv_table = s.get_nonredundant_snv_table()
    snv_table = snv_table[snv_table.cryptic == False]
    if 'morphia' in snv_table.columns:
        snv_table = snv_table[(snv_table.morphia == 2) | (snv_table.morphia == '2')]
        snv_table = snv_table.drop(columns="morphia")
    snv_table = snv_table.drop(columns="cryptic")
    snv_table, snp_types = characterize_snp(gene_index, gene_starts, seqs, gene_direction, snv_table)

    # Store
    s.store('SNP-calls', snv_table, 'pandas', 'Function of SNPs (N vs S vs O)')
    snv_table.to_csv(out_base  + "aa-SNVs.tsv",sep='\t', quoting=csv.QUOTE_NONE)

    print("Updating linkage table...")
    linkage_table = s.get_nonredundant_linkage_table()

    for index,snp_pair in linkage_table.iterrows():
        if snp_pair['scaffold'] + ":" + str(snp_pair['position_A']) in snp_types:
            linkage_table.at[index,'mutation_type_A'] = snp_types[snp_pair['scaffold'] + ":" + str(snp_pair['position_A'])]
        if snp_pair['scaffold'] + ":" + str(snp_pair['position_B']) in snp_types:
            linkage_table.at[index,'mutation_type_B'] = snp_types[snp_pair['scaffold'] + ":" + str(snp_pair['position_B'])]


    if len(linkage_table) > 0:
        linkage_table.to_csv(out_base + "aa-linkage.tsv",sep='\t', quoting=csv.QUOTE_NONE)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """
        A script that runs two different gene-based analyses on inStrain output:\n
        (1) Calculates nucleotide diversity for each gene provided.
        (2) Determines if SNPs are non-synonymous or synonymous, and the AA changed, adds this info to .linkage and .snp tables.
        \n

        Required input: a prodigal .fna gene calls file (created with prodigal -d). Going to be opinionated about this and not support alternatives.
        This file is easy to recreate from other / custom data - please look under examples to see how it should be formatted.
        """, formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument('input', help="an inStrain object prefix. (don't include `.`, only prefixes - must be in the local directory).")
    parser.add_argument("-g", "--gene_file", action="store", required=True, \
        help='Path to prodigal .fna genes file.')

    args = parser.parse_args()
    main(args)
