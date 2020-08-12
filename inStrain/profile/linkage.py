'''
Handles code related to linkage
'''

import itertools
import collections
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict

import inStrain.profile.profile_utilities

def calc_mm_SNV_linkage_network(read_to_snvs, scaff=False):
    '''
    Calculates the SNV linkage network

    Arguments:
        read_to_snvs = mm -> read -> SNVs

    saves it as a dictionary of edges in self.snv_net.
    Writes it out to a file, genome.net for reading in through other programs like node2vec
    '''

    G = nx.Graph() # graph between positions. has the attribute mm2combo2counts
    for mm, tread2snvs in read_to_snvs.items():

        for read in tread2snvs:
            snvs = tread2snvs[read]
            for pair in itertools.combinations(snvs, 2):
                p1 = int(pair[0].split(":")[0])
                p2 = int(pair[1].split(":")[0])

                c = "{0}:{1}".format(pair[0].split(":")[1], pair[1].split(":")[1])
                if not G.has_edge(p1, p2):
                    G.add_edge(p1, p2)
                    G[p1][p2]['mm2combo2counts'] = {}
                if mm not in G[p1][p2]['mm2combo2counts']:
                    G[p1][p2]['mm2combo2counts'][mm] = {}
                if c not in G[p1][p2]['mm2combo2counts'][mm]:
                    G[p1][p2]['mm2combo2counts'][mm][c] = 0
                G[p1][p2]['mm2combo2counts'][mm][c] += 1

    return G

def calculate_ld(mm_to_position_graph, min_snp,
                snv2mm2counts, scaffold):
    '''
    Calculates Linkage Disequilibrium for all SNVs in a window.
    '''
    r2_table = defaultdict(list)

    # This just gives potisions
    for edge in mm_to_position_graph.edges(list(snv2mm2counts.keys())):

        # For this position, now you have combo (allele combo) to mm (read mismatches) to counts (reads wih that combo and level of mismatches)
        mm2combo2counts = mm_to_position_graph.edges[edge[0],edge[1]]['mm2combo2counts']
        # Proceed with the full
        for mm, ld_args in _iterator_ld_sites(mm2combo2counts, min_snp, snv2mm2counts,
                    p1=int(edge[0]),p2=int(edge[1])):

            ld_result = _calc_ld_single(min_snp, **ld_args)

            if ld_result is False:
                continue

            ld_result['distance'] = abs(int(edge[0]) - int(edge[1]))
            ld_result['position_A'] = int(edge[0])
            ld_result['position_B'] = int(edge[1])

            _update_r2(ld_result, r2_table, mm, scaffold)

    # Create r2 linkage table
    r2linkage_table = pd.DataFrame(r2_table)
    return r2linkage_table

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
def _iterator_ld_sites(mm2combo2counts, min_snp, snv2mm2counts, p1, p2):
    '''
    pass
    '''
    # Get the tables ready
    s1 = snv2mm2counts[p1]
    s2 = snv2mm2counts[p2]

    # Figure out when you need to update
    updateMMs = set(s1.keys()).intersection(set(s2.keys()))

    counts = {0:np.zeros(4), 1:np.zeros(4)}
    combocounts = defaultdict(int)


    for mm, combo2counts in sorted(mm2combo2counts.items()):
        # update combo counts
        for combo, count in combo2counts.items():
            combocounts[combo] += count

        # update counts
        if mm not in updateMMs:
            continue
        counts[0] = inStrain.profile.profile_utilities.mm_counts_to_counts(s1, mm)
        counts[1] = inStrain.profile.profile_utilities.mm_counts_to_counts(s2, mm)

        if sum(counts[0]) + sum(counts[1]) >= min_snp:
            skip = False
            arg_dic = {}

            arg_dic['A'], arg_dic['a'] = major_minor_allele(counts[0])
            #assert arg_dic['A'] != arg_dic['a']

            arg_dic['B'], arg_dic['b'] = major_minor_allele(counts[1])
            #assert arg_dic['B'] != arg_dic['b']

            for val, key in zip(['A', 'a', 'B', 'b'], [0, 0, 1, 1]):
                count = int(counts[key][P2C[arg_dic[val]]])
                arg_dic[val + '_counts'] = count
                if count == 0:
                    skip = True

            if skip:
                continue

            total = 0
            for a in ['A', 'a']:
                for b in ['B', 'b']:
                    val = combocounts[arg_dic[a] + ':' + arg_dic[b]]
                    arg_dic[a + b] = val
                    total += val
            arg_dic['total'] = total

            yield mm, arg_dic

def major_minor_allele(counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3]}
    l = sorted(d, key=d.get, reverse=True)
    return l[0], l[1]

def _calc_ld_single(min_snp, **kwargs):
    '''
    A function that calculates the LD between two SNPs

    Argumements:
        snp_a = position of SNP A
        snp_b = position of SNP B

    Notes:
        allele_A = major allele at position A
        allele_a = minor allele at position A
    '''

    #distance between snps

    # figure out alleles
    allele_A = kwargs.get('A')
    allele_a = kwargs.get('a')
    allele_B = kwargs.get('B')
    allele_b = kwargs.get('b')
    countAB = kwargs.get('AB', 0)
    countAb = kwargs.get('Ab', 0)
    countaB = kwargs.get('aB', 0)
    countab = kwargs.get('ab', 0)

    total = kwargs.get('total')

    #Requires at least min_snp
    if total > min_snp:

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
            r2 = np.nan
        else:
            r2 = linkD*linkD / (freq_A * freq_a * freq_B * freq_b)

        linkd = freq_ab - freq_a * freq_b

        # calc D-prime
        d_prime = np.nan
        if (linkd < 0):
            denom = max([(-freq_A*freq_B),(-freq_a*freq_b)])
            d_prime = linkd / denom

        elif (linkD > 0):
            denom = min([(freq_A*freq_b), (freq_a*freq_B)])
            d_prime = linkd / denom

        ################
        # calc rarefied

        rareify = np.random.choice(['AB','Ab','aB','ab'], replace=True, p=[freq_AB,freq_Ab,freq_aB,freq_ab], size=min_snp)
        freq_AB = float(collections.Counter(rareify)['AB']) / min_snp
        freq_Ab = float(collections.Counter(rareify)['Ab']) / min_snp
        freq_aB = float(collections.Counter(rareify)['aB']) / min_snp
        freq_ab = float(collections.Counter(rareify)['ab']) / min_snp

        freq_A = freq_AB + freq_Ab
        freq_a = freq_ab + freq_aB
        freq_B = freq_AB + freq_aB
        freq_b = freq_ab + freq_Ab

        linkd_norm = freq_ab - freq_a * freq_b

        if freq_a == 0 or freq_A == 0 or freq_B == 0 or freq_b == 0:
            r2_normalized = np.nan
        else:
            r2_normalized = linkd_norm*linkd_norm / (freq_A * freq_a * freq_B * freq_b)



        # calc D-prime
        d_prime_normalized = np.nan
        if (linkd_norm < 0):
            denom = max([(-freq_A*freq_B),(-freq_a*freq_b)])
            d_prime_normalized = linkd_norm / denom

        elif (linkd_norm > 0):
            denom = min([(freq_A*freq_b), (freq_a*freq_B)])
            d_prime_normalized = linkd_norm / denom

        rt_dict = {}
        for att in ['r2', 'd_prime', 'r2_normalized', 'd_prime_normalized', 'total', 'countAB', \
                    'countAb', 'countaB', 'countab', 'allele_A', 'allele_a', \
                    'allele_B', 'allele_b']:
            rt_dict[att] = eval(att)


        return rt_dict

    else:
        return False

def _update_r2(ld_result, r2_table, mm, scaffold, window=False):
    if ld_result == False:
        return r2_table

    for key, value in ld_result.items():
        r2_table[key].append(value)
    r2_table['mm'].append(mm)
    r2_table['scaffold'].append(scaffold)
    #r2_table['window'].append("{0}:{1}_{2}".format(scaffold, window[0], window[1]))

    return r2_table

def update_linked_reads(read_to_snvs, pileupcolumn, MMcounts, position,
                        bases, R2M, min_freq=0.05, scaffold=False):
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

            try:
                val = pileupread.alignment.query_sequence[pileupread.query_position]

                # if value is not the consensus value
                if val in bases:
                    # this read is of a variant position
                    read_to_snvs[mm][read_name].append(position + ":" + val)
            except KeyError: # This would be like an N or something not A/C/T/G
                pass
