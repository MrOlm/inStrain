#!/usr/bin/env python

import os
import sys
import shutil
import argparse
import subprocess
import pandas as pd

from Bio import SeqIO
from collections import defaultdict

import inStrain
import inStrain.genomeUtilities

def main(args):
    '''
    The main controller of the program
    '''
    # Parse and validate arguments
    coverm_loc, genome2length = parse_validate(args)

    # Run coverm
    Cdb = run_coverm(coverm_loc, args)

    # Calculate on the genome level
    CGdb = parse_coverm(Cdb, genome2length, args)

    # Produce output
    produce_output(CGdb, args)

def parse_validate(args):
    '''
    Make sure all is well
    '''
    # Make sure coverM is working
    loc = shutil.which('coverm')
    works = False
    if loc != None:
        try:
            o = subprocess.check_output([loc, '-h'],stderr= subprocess.STDOUT)
            works = True
        except:
            pass
    if not works:
        print("Cannot find coverm; make sure its installed")
        sys.exit()

    # Set up a results folder
    if args.output is None:
        args.output = os.path.basename(args.bam)[:-4]
    try:
        os.mkdir(args.output)
    except:
        print("Couldnt make the output directory {0}".format(args.output))
    if args.output[-1] != '/':
        args.output += '/'

    # Get genome to length
    scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    s2l = {s: len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}
    del scaff2sequence

    # Set up the stb
    #args.stb = parse_stb(args.stb)
    args.stb = inStrain.genomeUtilities.load_scaff2bin(args.stb)

    if args.stb == {}:
        args.stb = {s:'all_scaffolds' for s, l in s2l.items()}

    genome2length = {}
    for scaffold, length in s2l.items():
        if scaffold not in args.stb:
            continue
        genome = args.stb[scaffold]
        if genome not in genome2length:
            genome2length[genome] = 0
        genome2length[genome] += (s2l[scaffold])# - 150) # 150 to account for the ignored ends; but it doesn't do that anymore

    return loc, genome2length

def run_coverm(exe_loc, args):
    '''
    run coverm
    '''
    coverm_cmd = "coverm contig -b {0} -t {1} -m mean covered_bases length count --output-format sparse".format(
                args.bam, args.processes)
    coverm_cmd += " | awk 'NR == 1 {{print $0}} NR != 1 {{if ($4/$5 > {0}) print $0}}'".format(args.stringent_breadth_cutoff)

    coverm_cmd += " > {0}".format(args.output + 'coverm_raw.tsv')

    print(coverm_cmd)
    subprocess.call(coverm_cmd, shell=True)

    Cdb = pd.read_csv(args.output + 'coverm_raw.tsv', sep='\t')
    return Cdb

def parse_coverm(Cdb, genome2length, args):
    '''
    Make this genome-wide
    '''
    # Load the stb file
    stb = args.stb
    assert type(stb) == type({})

    # Apply it
    Cdb.loc[:,'genome'] = Cdb['Contig'].map(stb)

    # Print results
    inS = len(Cdb[~Cdb['genome'].isna()]['Contig'].unique())
    tS = len(Cdb['Contig'].unique())
    print("{0:,} of {1:,} scaffolds are in the .stb ({2:.1f}%)".format(inS, tS, (inS/tS)*100))

    inR = Cdb[~Cdb['genome'].isna()]['Read Count'].sum()
    tR = Cdb['Read Count'].sum()
    print("{0:,} of {1:,} mapping reads are in the .stb ({2:.1f}%)".format(inR, tR, (inR/tR)*100))

    # Calculate breadth, coverage, and total reads of genomes
    table = defaultdict(list)
    for genome, db in Cdb.groupby('genome'):
        genome_length = genome2length[genome]
        covered_bases = db['Covered Bases'].sum() / genome_length
        averaged_coverage = sum([b * l for b, l in zip(db['Mean'], db['Length'])]) / genome_length

        table['genome'].append(genome)
        table['length'].append(genome_length)
        table['breadth'].append(covered_bases)
        table['coverage'].append(averaged_coverage)
        table['reads'].append(db['Read Count'].sum())
    CGdb = pd.DataFrame(table)

    return CGdb

def produce_output(CGdb, args):
    '''
    Make stuff
    '''
    # Save the parsed file
    CGdb.to_csv(args.output + 'genomeCoverage.csv', index=False)

    # Save the scaffolds
    if len(CGdb) > 0:
        genomes = set(CGdb[CGdb['breadth'] >= args.breadth_cutoff]['genome'].unique())
        with open(args.output + 'scaffolds.txt', 'w') as o:
            for scaffold, bin in args.stb.items():
                if bin in genomes:
                    o.write(scaffold + '\n')


# if __name__ == '__main__':
#
#     parser = argparse.ArgumentParser(description= """
#         A script to quickly determine breadth and coverage of a .bam file\n
#         """, formatter_class=argparse.RawTextHelpFormatter)
#
#     # Required positional arguments
#     parser.add_argument('-b', '--bam', help="A bam file to profile",
#                         required=True)
#     parser.add_argument('-f', '--fasta', help="The .fasta file to profile",
#                         required=True)
#     parser.add_argument('-s', '--stb', help="Scaffold to bin file",
#                         required=True)
#     parser.add_argument("-o", "--output", action="store", \
#                         help='Output prefix')
#     parser.add_argument("-p", "--processes", action="store", default=6, type=int, \
#                         help='Threads to use for multiprocessing')
#
#     args = parser.parse_args()
#     main(args)
