'''
This houses functions that are from an old version of inStrain; useful for testing
'''

import inStrain.controller

import inStrain.deprecated.deprecated_filter_reads

import inStrain.profile
import inStrain.profile.profile_utilities
import inStrain.profile.snv_utilities
# import inStrain.profileUtilities

P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}  # base -> position


# THIS IS CC'S VERSION
def main(args):
    '''
    Main entry point
    '''

    # Set up Logger
    if not args.output:
        log = args.fasta.split("/")[-1].split(".")[0] + '.log'
    else:
        log = args.output + '.log'
    inStrain.controller.setup_logger(log)

    if args.filter_cutoff:
        if float(args.filter_cutoff) < 1:
            strain_pipeline(args, float(args.filter_cutoff))
        else:
            logging.error("Error: Strain level should be a number between 0 and 1.")
            sys.exit(1)
    else:
        strain_pipeline(args, 0.96)


def strain_pipeline(args, filter_cutoff):
    strains = SNVdata()

    # if args.testing:
    #     strains.testing = True

    if not args.output:
        strains.output = args.fasta.split("/")[-1].split(".")[0] + "_" + str(filter_cutoff)
    else:
        strains.output = args.output + "_" + str(filter_cutoff)

    logging.info("Running at resolution: (>" + str(filter_cutoff) + "%)")

    args.genes = None
    strains.get_scaffold_positions(args.genes, args.fasta)
    strains.run_strain_profiler(args.bam, min_cov=int(args.min_cov), min_freq=float(args.min_freq),
                                filter_cutoff=filter_cutoff)

    logging.info("\nCalculating linkage network...\n")
    strains.calc_linkage_network()

    logging.info("Calculating LD...\n")
    strains.calc_ld_all_sites(int(args.min_snp))
    strains.save()


class SNVdata:
    # The class "constructor" - It's actually an initializer
    def __init__(self):

        # Parameters
        self.fasta = ""
        self.bam = ""
        self.results = False
        self.output = None
        self.testing = False
        self.positions = []
        self.fasta_length = 0

        # Data structures
        self.snv_table = None  #
        self.read_to_snvs = None  #
        self.windows_to_snvs = None  #
        self.all_counts = None  # All counts of AGCT for each position
        self.snv_counts = None  # Counts of AGCT for each SNV position
        self.clonality_table = None  # Dict of clonality histograms (lists) by window
        self.snv_graph = None  # Weighted networkx graph that tracks SNVs( 100:A, 100:T are diff nodes) that are linked
        self.position_graph = None  # Unweighted networkx graph that tracks positions that are linked
        self.r2linkage_table = None  # Dict of r2 histograms (list of lists) by window

        # General statistics
        self.coverages = None
        self.total_positions = None
        self.total_snv_sites = None
        self.alpha_snvs = None

    def save(self, size=0):
        if self.results:
            # Generate tables
            if size == 0:
                self.read_to_snvs = None

            logging.info("printing to tables.")
            sample = self.bam.split("/")[-1].split(".bam")[0]

            genome = self.output
            logging.info(genome)
            self.snv_table.to_csv(genome + ".freq", sep='\t')  # , quoting=csv.QUOTE_NONE)
            self.clonality_table.to_csv(genome + ".clonal", sep='\t')  # , quoting=csv.QUOTE_NONE)
            self.r2linkage_table.to_csv(genome + ".linkage", sep='\t')  # , quoting=csv.QUOTE_NONE)
            self.mapping_info.to_csv(genome + ".readFiltering", sep='\t')  # ', quoting=csv.QUOTE_NONE)

            f = open(genome + ".data", 'wb')
            pickle.dump(self.__dict__, f, 2)
            f.close()
        else:
            logging.info("No data to save.")

    def load(self, name):
        self.__dict__.clear()
        f = open(name + ".data", 'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)

    def get_scaffold_positions(self, gene_list=None, fasta_file=None):
        ''' Returns a list of windows to record SNVs in'''
        if not fasta_file and not gene_list:
            #            logging.info("ERROR: REQUIRED TO SUPPLY FASTA or GENE FILE")
            sys.exit(1)

        if gene_list:
            f = open(gene_list)
            for line in f.readlines():
                self.positions.append([line.split(",")[1], int(line.split(",")[2]), int(line.split(",")[3])])
            f.close()
        else:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                self.fasta_length += len(rec.seq)
                start = 0
                while True:
                    chunk = start + 15000
                    if chunk <= len(rec.seq):
                        self.positions.append([str(rec.id), start, start + 10000])
                    else:
                        self.positions.append([str(rec.id), start, len(rec.seq)])
                        break
                    start += 10000
                    start += 1

        self.fasta = fasta_file

    def get_snps_per_gene(self, gene_file=None):
        pass

    def calc_linkage_network(self):
        ''' Calculates the SNV linkage network - saves it as a dictionary of edges in self.snv_net.
        Writes it out to a file, genome.net for reading in through other programs like node2vec
        '''
        G = nx.Graph()  # graph between SNPs  (100:A and 100:C are different nodes)
        G_pos = nx.Graph()  # graph between positions; unweighted

        logging.info("Calculating SNV linkage network...")
        for read in self.read_to_snvs:
            snvs = self.read_to_snvs[read]
            for pair in itertools.combinations(snvs, 2):
                if G.has_edge(pair[0], pair[1]):
                    G[pair[0]][pair[1]]['weight'] += 1
                else:
                    G.add_edge(pair[0], pair[1], weight=1)
                    if not G_pos.has_edge(pair[0].split(":")[0], pair[1].split(":")[0]):
                        G_pos.add_edge(pair[0].split(":")[0], pair[1].split(":")[0])

        logging.info(str("Of " + str(self.total_snv_sites) + " SNP-sites, there were " + str(
            len(G_pos)) + " SNPs that could be linked to at least one other SNP."))
        logging.info("The average SNP was linked to " + str(
            int(np.mean([int(x[1]) for x in list(G_pos.degree())]))) + " other SNPs.")
        self.snv_graph = G
        self.position_graph = G_pos

    def calc_ld(self, snp_a, snp_b, min_snp, report=False):
        '''
        A function that calculates the LD between two SNPs
        '''

        # distance between snps
        distance = abs(int(snp_a.split("_")[-1]) - int(snp_b.split("_")[-1]))

        # calculate allele frequencies
        allele_A = major_allele(snp_a, self.snv_counts[snp_a])  # get major allele
        allele_a = minor_allele(snp_a, self.snv_counts[snp_a])  # get minor allele
        freq_A = float(allele_A[1]) / (allele_A[1] + allele_a[1])
        freq_a = float(allele_a[1]) / (allele_A[1] + allele_a[1])

        # get snp B frequencies
        allele_B = major_allele(snp_b, self.snv_counts[snp_b])
        allele_b = minor_allele(snp_b, self.snv_counts[snp_b])
        freq_B = float(allele_B[1]) / (allele_B[1] + allele_b[1])
        freq_b = float(allele_b[1]) / (allele_B[1] + allele_b[1])

        # Get frequencies of linkages
        countAB, countAb, countaB, countab = 0, 0, 0, 0
        if self.snv_graph.has_edge(allele_A[0], allele_B[0]):
            countAB = self.snv_graph[allele_A[0]][allele_B[0]]['weight']
        if self.snv_graph.has_edge(allele_A[0], allele_b[0]):
            countAb = self.snv_graph[allele_A[0]][allele_b[0]]['weight']
        if self.snv_graph.has_edge(allele_a[0], allele_B[0]):
            countaB = self.snv_graph[allele_a[0]][allele_B[0]]['weight']
        if self.snv_graph.has_edge(allele_a[0], allele_b[0]):
            countab = self.snv_graph[allele_a[0]][allele_b[0]]['weight']

        total = countAB + countAb + countaB + countab

        if report:
            print("total {0}".format(total))
            print([countAB, countAb, countaB, countab])

        # Requires at least min_snp
        if total > min_snp:

            linkage_points_x = []
            linkage_points_y = []
            for point in range(0, countAB):
                linkage_points_x.append(1)
                linkage_points_y.append(1)
            for point in range(0, countAb):
                linkage_points_x.append(1)
                linkage_points_y.append(0)
            for point in range(0, countaB):
                linkage_points_x.append(0)
                linkage_points_y.append(1)
            for point in range(0, countab):
                linkage_points_x.append(0)
                linkage_points_y.append(0)

            freq_AB = float(countAB) / total
            freq_Ab = float(countAb) / total
            freq_aB = float(countaB) / total
            freq_ab = float(countab) / total

            freq_A = freq_AB + freq_Ab
            freq_a = freq_ab + freq_aB
            freq_B = freq_AB + freq_aB
            freq_b = freq_ab + freq_Ab

            linkD = freq_AB - freq_A * freq_B

            if freq_a == 0 or freq_A == 0 or freq_B == 0 or freq_b == 0:
                r2_book = np.nan
            else:
                r2_book = linkD * linkD / (freq_A * freq_a * freq_B * freq_b)

            linkd = freq_ab - freq_a * freq_b
            # calc D-prime
            if (linkd < 0):
                denom = max([(-freq_A * freq_B), (-freq_a * freq_b)])
                d_prime = linkd / denom

            elif (linkD > 0):
                denom = min([(freq_A * freq_b), (freq_a * freq_B)])
                d_prime = linkd / denom
            else:
                d_prime = 0

            # r2 = stats.pearsonr(linkage_points_x, linkage_points_y)[0]
            # r2 = r2 * r2
            # logging.info(r2)
            # logging.info(r2_book)
            r2 = r2_book
            return ([distance, r2, linkD, linkd, r2_book, total, countAB, countAb, countaB, countab, allele_A, allele_a,
                     allele_B, allele_b])
            # else:
            #     logging.info("nan")
            # else:
        else:
            return False

    def calc_ld_all_sites(self, min_snp):
        '''
        Calculates Linkage Disequilibrium for all SNVs in a window.
        '''
        r2_total = {}
        sigma_total = {}
        logging.info("Calculating Linkage Disequilibrium")

        for window in self.positions:
            r2_spectrum = []
            sigma_spectrum = []

            window_name = window[0] + ":" + str(window[1]) + ":" + str(window[2])
            window_snvs = self.windows_to_snvs[window_name]

            for edge in self.position_graph.edges(window_snvs):
                snp_a = edge[0]
                snp_b = edge[1]

                ld_result = self.calc_ld(snp_a, snp_b, min_snp)
                if ld_result:
                    r2_spectrum.append(ld_result)

            r2_total[window_name] = r2_spectrum

        # Create r2 linkage table
        r2linkage_table = defaultdict(list)

        for window in r2_total:
            for datum in r2_total[window]:
                r2linkage_table['Window'].append(window)
                r2linkage_table['Distance'].append(datum[0])
                r2linkage_table['r2'].append(datum[1])
                r2linkage_table['linkD'].append(datum[2])
                # r2linkage_table['linkd'].append(datum[3])
                # r2linkage_table['r2_book'].append(datum[4])
                r2linkage_table['total'].append(datum[5])
                r2linkage_table['count_AB'].append(datum[6])
                r2linkage_table['count_Ab'].append(datum[7])
                r2linkage_table['count_aB'].append(datum[8])
                r2linkage_table['count_ab'].append(datum[9])

                r2linkage_table['total_A'].append(datum[10])
                r2linkage_table['total_a'].append(datum[11])
                r2linkage_table['total_B'].append(datum[12])
                r2linkage_table['total_b'].append(datum[13])

        self.r2linkage_table = pd.DataFrame(r2linkage_table)

    def run_strain_profiler(self, bam, min_cov=5, min_freq=0.05, filter_cutoff=0):
        '''
        Main class for finding SNVs and generating data profile for a genome.
        '''
        minimum_mapq = 2
        P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
        C2P = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}

        # Set up variables
        alpha_snvs = 0
        total_positions = 0
        total_snv_sites = 0

        raw_counts_data = defaultdict(dict)  # Set up SNP table
        read_to_snvs = defaultdict(list)  # read name -> position : variant base
        snvs_frequencies = defaultdict(int)
        clonality_by_window = defaultdict(list)  # window -> [position, clonality]
        windows_to_snvs = defaultdict(list)  # window -> position

        all_counts = {}  # position -> coverage
        snv_counts = {}  # position -> [A, C, T, G]
        coverages = {}  # position -> coverage

        # read null model
        null_loc = os.path.dirname(__file__) + '/../helper_files/NullModel.txt'
        null_snp_model = inStrain.profile.snv_utilities.generate_snp_model(null_loc)

        ## Start reading BAM
        samfile = pysam.AlignmentFile(bam)
        sample = bam.split("/")[-1].split(".bam")[0]

        if self.testing:
            self.positions = self.positions[0:10]

        # Call read filtering function
        subset_reads, Rdb = inStrain.deprecated.deprecated_filter_reads.filter_reads(bam, self.positions, self.fasta_length,
                                                                          filter_cutoff, 3, 50, 2)

        # Save the read filtering report
        self.mapping_info = Rdb

        ### START SNP FINDING

        ## Start looping through each region
        for gene in tqdm(self.positions, desc='Finding SNVs ...'):
            scaff = gene[0]
            window = gene[0] + ":" + str(gene[1]) + ":" + str(gene[2])
            for pileupcolumn in samfile.pileup(scaff, truncate=True, max_depth=100000, stepper='nofilter',
                                               compute_baq=True, ignore_orphans=True, ignore_overlaps=True,
                                               min_base_quality=30):
                # for pileupcolumn in samfile.pileup(scaff, gene[1], gene[2], truncate = True, max_depth=100000, stepper = 'nofilter', compute_baq= True, ignore_orphans = True, ignore_overlaps = True,  min_base_quality = 30):
                ## Step 1: Are there any reads at this position?

                position = scaff + "_" + str(pileupcolumn.pos)
                counts = _get_base_counts(pileupcolumn, filtered_reads=subset_reads)

                consensus = False
                # Yes there were reads at this position
                if counts:
                    all_counts[position] = counts

                    if sum(counts) >= min_cov:
                        total_positions += 1
                        pos_clonality = inStrain.profile.snv_utilities.calculate_clonality(counts)
                        clonality_by_window[window].append([position, pos_clonality])
                        consensus = call_snv_site_old(counts, null_snp_model, min_cov=min_cov, min_freq=min_freq)
                        coverages[position] = sum(counts)

                # SEARCH HERE FOR DEBUG
                # if position == "N5_271_010G1_scaffold_114_990":
                #     print(counts)
                #     #counts2 = _get_base_counts_(pileupcolumn, filtered_reads = subset_reads)
                #     reads = []
                #     for pileupread in pileupcolumn.pileups:
                #         if not pileupread.is_del and not pileupread.is_refskip:
                #             if pileupread.alignment.query_name in subset_reads:
                #                 read_name = pileupread.alignment.query_name
                #                 reads.append((read_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                #
                #     print(len(reads))
                #     print('\n'.join([str(r) for r in reads]))
                #     print(consensus)

                ## Strep 2: Is there an SNV at this position?
                if consensus:
                    # min_snp for this site
                    local_min_snp = null_snp_model[sum(counts)]
                    # there's an SNV at this site
                    total_snv_sites += 1
                    windows_to_snvs[window].append(position)
                    snv_counts[position] = counts

                    # Get reads to snvs
                    for pileupread in pileupcolumn.pileups:
                        read_name = pileupread.alignment.query_name
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if read_name in subset_reads:
                                try:
                                    val = pileupread.alignment.query_sequence[pileupread.query_position]
                                    # if value is not the consensus value
                                    if counts[P2C[val]] >= local_min_snp and float(counts[P2C[val]]) / sum(
                                            counts) >= min_freq:
                                        # this is a variant read!
                                        read_to_snvs[read_name].append(position + ":" + val)
                                        if val != consensus:
                                            alpha_snvs += 1
                                except KeyError:  # This would be like an N or something not A/C/T/G
                                    pass

                    # Add to frequencies
                    nucl_count = 0
                    for nucl in counts:
                        if nucl >= local_min_snp:
                            freq = float(nucl) / float(sum(counts))
                            if freq >= min_freq:
                                snp = position + ":" + C2P[nucl_count]
                                snvs_frequencies[snp] = [freq, window]

                        nucl_count += 1

        #### SNP FINDING COMPLETE
        # Create SNV frequency table
        snv_table = defaultdict(list)
        for snv in snvs_frequencies:
            snv_table['SNV'].append(snv)
            snv_table['freq'].append(snvs_frequencies[snv][0])
            snv_table['Window'].append(snvs_frequencies[snv][1])
        snv_table = pd.DataFrame(snv_table)

        # Create clonality table
        clonality_table = defaultdict(list)
        for window in clonality_by_window:
            for position_pair in clonality_by_window[window]:
                clonality_table['scaffold'].append(window.split(":")[0])
                clonality_table['window_name'].append(window)
                clonality_table['position'].append(position_pair[0])
                clonality_table['clonality'].append(position_pair[1])
                clonality_table['coverage'].append(coverages[position_pair[0]])

        clonality_table_final = pd.DataFrame(clonality_table)

        # Final statistics
        logging.debug("Total SNVs-sites: " + str(total_snv_sites))
        logging.debug("Total SNV-bases: " + str(alpha_snvs))
        logging.debug("Mean clonality: " + str(
            float(sum(clonality_table['clonality'])) / float(len(clonality_table['clonality']))))
        logging.debug("Total sites above min-coverage: " + str(total_positions))
        logging.debug("Mean coverage: " + str(float(sum(coverages.values())) / len(coverages)))
        logging.debug("Total number of bases: " + str(sum(coverages.values())))

        self.coverages = coverages
        self.alpha_snvs = alpha_snvs
        self.total_snv_sites = total_snv_sites
        self.clonality_table = clonality_table_final
        self.snv_table = snv_table
        self.read_to_snvs = read_to_snvs
        self.total_positions = total_positions
        self.windows_to_snvs = windows_to_snvs
        self.all_counts = all_counts
        self.snv_counts = snv_counts

        self.results = True


def _get_base_counts(pileupcolumn, filtered_reads):
    '''
    From a pileupcolumn object, return a list with the counts of [A, C, T, G]
    '''
    P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
    C2P = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}

    counts = [0, 0, 0, 0]
    empty = [0, 0, 0, 0]

    for pileupread in pileupcolumn.pileups:
        # logging.info(pileupread.)

        if not pileupread.is_del and not pileupread.is_refskip:
            read_name = pileupread.alignment.query_name
            if read_name in filtered_reads:
                # if pileupread.indel != 0:
                #     logging.info(pileupread.alignment.get_aligned_pairs())
                try:
                    counts[P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
                except KeyError:  # This would be like an N or something not A/C/T/G
                    pass
    if counts != empty:
        return counts
    else:
        return False


def call_snv_site_old(counts, null_model, min_cov=5, min_freq=0.05):
    '''
    Determines whether a site has a variant based on its nucleotide count frequencies.
    '''
    P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
    C2P = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
    total = sum(counts)

    if total >= min_cov:
        i = 0
        for c in counts:
            if c >= null_model[total] and float(c) / total >= min_freq:
                i += 1
        if i > 1:
            return C2P[counts.index(max(counts))]
    else:
        return False


def major_allele(snv, counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3]}
    nucl = sorted(d, key=d.get, reverse=True)[0]
    return [snv + ":" + nucl, d[nucl]]


def major_allele2(counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3]}
    return sorted(d, key=d.get, reverse=True)[0]


def minor_allele(snv, counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3]}
    nucl = sorted(d, key=d.get, reverse=True)[1]
    return [snv + ":" + nucl, d[nucl]]


def minor_allele2(counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3]}
    return sorted(d, key=d.get, reverse=True)[1]


'''
           __..--''``---....___   _..._    __
 /// //_.-'    .-/";  `        ``<._  ``.''_ `. / // /
///_.-' _..--.'_    \                    `( ) ) // //
/ (_..-' // (< _     ;_..__               ; `' / ///
 / // // //  `-._,_)' // / ``--...____..-' /// / //

 Below this is super deprecated (used nowhere)
'''

# !/usr/bin/env python

'''
Possible names for this program:
* strainRep - not a legit acronym
* blockStrain - not a legit acronym
* destrain - not a legit acronym
* instrain - not a legit acronym
* metaIPD - metagenome interference of population structure - a legt acronym
'''

# Get the version
from inStrain._version import __version__

# Import
import os
import sys
import math
import pysam
import pickle
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict

## local imports
import inStrain.filter_reads


def entropy2(counts):
    probs = []
    total = sum(counts)
    for c in counts:
        probs.append(float(c) / total)

    ent = 0
    for i in probs:
        if i != 0:
            ent -= i * math.log(i, math.e)
    return ent


def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts


def iterate_windows(Wdb, snv_table, mm):
    '''
    Iterate windows
    '''
    Sdb = snv_table[snv_table['mm'] >= mm].sort_values('mm') \
        .drop_duplicates(subset=['scaffold', 'position'], keep='last') \
        .sort_index().drop(columns=['mm'])
    Sdb = Sdb.set_index('position', drop=False)

    for i, row in Wdb.iterrows():
        window = (int(row['start']), int(row['end']))
        sdb = Sdb[(Sdb['position'] >= window[0]) & (Sdb['position'] <= window[1])]
        snvs = list(sdb['position'].tolist())
        if len(snvs) > 0:
            yield snvs, window, sdb


def _get_series(s1, mm, lms1):
    # already a series
    if type(s1) == pd.Series:
        if s1['mm'] >= mm:
            return s1
        else:
            return lms1
    else:
        try:
            return s1.loc[mm]
        except:
            return lms1


def _update_counts(counts, ts1, ts2, mm):
    for b in ['A', 'C', 'T', 'G']:
        if len(ts1) > 0:
            counts[0][P2C[b]] += ts1[b]
        if len(ts2) > 0:
            counts[1][P2C[b]] += ts2[b]


def calc_SNV_linkage_network(read_to_snvs):
    '''
    Calculates the SNV linkage network

    Arguments:
        read_to_snvs = mm -> read -> SNVs

    saves it as a dictionary of edges in self.snv_net.
    Writes it out to a file, genome.net for reading in through other programs like node2vec
    '''
    mm_to_snv_graph = {}
    mm_to_position_graph = {}

    # logging.info("Calculating SNV linkage networks...")
    read2snvs = {}
    for mm, tread2snvs in read_to_snvs.items():
        G = nx.Graph()  # graph between SNPs  (100:A and 100:C are different nodes)
        G_pos = nx.Graph()  # graph between positions; unweighted

        read2snvs = {**read2snvs, **tread2snvs}
        for read in read2snvs:
            snvs = read2snvs[read]
            for pair in itertools.combinations(snvs, 2):
                if G.has_edge(pair[0], pair[1]):
                    G[pair[0]][pair[1]]['weight'] += 1
                else:
                    G.add_edge(pair[0], pair[1], weight=1)
                    if not G_pos.has_edge(int(pair[0].split(":")[0]), int(pair[1].split(":")[0])):
                        G_pos.add_edge(int(pair[0].split(":")[0]), int(pair[1].split(":")[0]))

        mm_to_snv_graph[mm] = G
        mm_to_position_graph[mm] = G_pos

    return mm_to_snv_graph, mm_to_position_graph


def _merge_tables_special(Cdb, Adb):
    FIX_COLS = ['ANI', 'SNPs', 'breadth_minCov']

    # Handle edge-case where there are no SNPs
    if len(Adb) == 0:
        for c in FIX_COLS:
            Cdb[c] = 0
        return Cdb

    # Make sure anything with a coverage has a SNP
    assert len(set(Adb['mm']) - set(Cdb['mm'])) == 0

    # Do initial merge
    Cdb = pd.merge(Cdb, Adb, on=['scaffold', 'mm'], how='outer')
    Cdb = Cdb.sort_values(['scaffold', 'mm'])
    Cdb = Cdb.reset_index(drop=True)

    # Run up NaN values. For example, in the cases
    Fdb = run_up_NaN(Cdb, FIX_COLS, on='scaffold')

    # Might as well adjust some datatypes
    pass

    return Fdb


def run_up_NaN(odb, cols, on='scaffold'):
    '''
    Take NaN values and fill them with the column above.
    For example, if you have mm of [0,1,2,3], and breadth of [.9, .92, NaN, .94],
    change the breadth to [.9, .92, .92, .94]
    If you have [0,1,2,3] and [NaN, NaN, .92, .94], it'll change it to:
    [0, 0, .92, .94]
    NOTE: Must be sorted / indexed in the order that you want run up
    Args:
        odb: original dataframe
        cols: columns to do this filling on
        on: the "key" of the dataframe. If this is all NaN, fill with 0s
    Returns:
        DataFrame: new dataframe
    '''
    Fdb = odb.copy()
    for scaff, db in odb.groupby(on):
        # handle edge case where all are NaN
        if len(db[cols[0]].dropna()) == 0:
            for i, row in db.iterrows():
                for col in cols:
                    Fdb.at[i, col] = 0
            continue

        # do the actual run-up
        top = True
        for i, row in db.iterrows():
            # hangle edge case where top values are NaN
            if top & np.isnan(row[cols[0]]):
                for col in cols:
                    Fdb.at[i, col] = 0
                continue
            else:
                top = False

            # The normal run-up case
            if np.isnan(row['ANI']):
                for col in cols:
                    Fdb.at[i, col] = Fdb.at[i - 1, col]

    return Fdb


def _zeros(mLen=1):
    return np.zeros(mLen, dtype=int)


def _clonT_to_clons(clonT, maxMM):
    '''
    Return the max clonality at each position
    '''
    counts = None
    for mm, count in [(mm, count) for mm, count in clonT.items() if mm <= maxMM]:
        if counts is None:
            counts = count
        else:
            counts = np.maximum(counts, count)

    if counts is None:
        return np.full(1, fill_value=np.nan)

    else:
        return counts


def make_ANI_table(snpsCounted, basesCounted, lengt, scaff, covT,
                   min_cov, Wdb):
    '''
    Fill in the SNP table with the SNPs and breadth_minCov for each scaffold and mm
    '''
    table = defaultdict(list)
    # fill in all SNP information
    for mm in sorted(list(covT.keys())):
        covs = inStrain.profile.profile_utilities.mm_counts_to_counts(covT, mm)
        if covs == [0, 0, 0, 0]:
            counted_basesO = 0
        else:
            zeros = (covs < min_cov).sum()
            counted_basesO = lengt - zeros

        counted_snps = _calc_counted_bases(snpsCounted, mm)
        counted_bases = _calc_counted_bases(basesCounted, mm)

        table['scaffold'].append(scaff)
        table['mm'].append(mm)
        table['SNPs'].append(counted_snps)
        # table['breadth_minCov'].append(counted_bases / lengt)
        if counted_bases == 0:
            table['ANI'].append(0)
        else:
            table['ANI'].append((counted_bases - counted_snps) / counted_bases)

    return pd.DataFrame(table)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="inStrain version {0}".format(__version__),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required positional arguments
    parser.add_argument("bam", help="Sorted .bam file")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")

    # Optional arguments
    parser.add_argument("-o", "--output", action="store", default='inStrain', \
                        help='Output prefix')
    parser.add_argument("-p", "--processes", action="store", default=6, type=int, \
                        help='Threads to use for multiprocessing')
    parser.add_argument("-c", "--min_cov", action="store", default=5, \
                        help='Minimum SNV coverage')
    parser.add_argument("-s", "--min_snp", action="store", default=20, \
                        help='Absolute minimum number of reads connecting two SNPs to calculate LD between them.')
    parser.add_argument("-f", "--min_freq", action="store", default=0.05, \
                        help='Minimum SNP frequency to confirm a SNV (both this AND the 0.  001 percent FDR snp count cutoff must be true).')
    parser.add_argument("-fdr", "--fdr", action="store", default=1e-6, type=float, \
                        help='SNP false discovery rate- based on simulation data with a 0.1 percent error rate (Q30)')
    parser.add_argument("--min_scaffold_reads", action="store", default=0, type=int, \
                        help='Minimum number of reads mapping to a scaffold to proceed with profiling it')
    parser.add_argument("--scaffolds_to_profile", action="store", \
                        help='Path to a file containing a list of scaffolds to profile- if provided will ONLY profile those scaffolds')

    # Read filtering cutoffs
    parser.add_argument("-l", "--filter_cutoff", action="store", default=0.95, type=float, \
                        help='Minimum percent identity of read pairs to consensus to use the reads. Must be >, not >=')
    parser.add_argument("--min_mapq", action="store", default=-1, type=int, \
                        help='Minimum mapq score of EITHER read in a pair to use that pair. Must be >, not >=')
    parser.add_argument("--max_insert_relative", action="store", default=3, type=float, \
                        help='Multiplier to determine maximum insert size between two reads - default is to use 3x median insert size. Must be >, not >=')
    parser.add_argument("--min_insert", action="store", default=50, type=int, \
                        help='Minimum insert size between two reads - default is 50 bp. If two reads are 50bp each and overlap completely, their insert will be 50. Must be >, not >=')

    parser.add_argument('--store_everything', action='store_true', default=False, \
                        help="Store intermediate dictionaries in the pickle file; will result in significantly more RAM and disk usage")
    parser.add_argument('--skip_mm_profiling', action='store_true', default=False, \
                        help="Dont perform analysis on an mm level; saves RAM and time")

    # parser.add_argument("-g", "--genes", action="store", default=None, \
    #     help='Optional genes file')
    parser.add_argument('--debug', action='store_true', default=False, \
                        help="Produce some extra output helpful for debugging")

    # Parse
    if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
        parser.print_help()
        sys.exit(0)
    else:
        return parser.parse_args(args)


if __name__ == '__main__':
    """ This is executed when run from the command line """

    pass
    # Controller().main(sys.argv[1:])
    # args = parse_arguments(sys.argv[1:])
    # main(args)


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
