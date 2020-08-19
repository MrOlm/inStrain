'''
This has utilities related to SNVs and SNSs
'''

import numpy as np
import pandas as pd

import inStrain
import inStrain.profile.profile_utilities

P2C = {'A':0, 'C':1, 'T':2, 'G':3} # base -> position
C2P = {0:'A', 1:'C', 2:'T', 3:'G'} # position -> base

def generate_snp_model(model_file, fdr=1e-6):
    '''
    The model_file contains the columns "coverage", "probability of X coverage base by random chance"

    The model that it returns has the properties:
        model[position_coverage] = the minimum number of alternate bases to be legit
        model[-1] = use this value if coverage is not in dictionary
    '''
    f = open(model_file)
    model = {}
    for line in f.readlines():
        if 'coverage' in line:
            continue
        counts = line.split()[1:]
        coverage = line.split()[0]
        i = 0 # note 6.4.20 - this should be 1, not zero
        for count in counts:
            if float(count) < fdr:
                model[int(coverage)] = i
                break
            i += 1

    model[-1] = max(model.values())

    return model

def update_snp_table(Stable, clonT, clonTR, MMcounts, p2c,\
        pos, scaff, mLen, ref_base, null_model, min_cov=5, min_covR=50, min_freq=.05):
    '''
    Add information to SNP table. Update basesCounted and snpsCounted

    Also update clonality

    cryptic means that it's a SNP at lower mm levels, and then not again at higher levels

    6.10.18 - v0.4.1 - Notes

    What makes it into the Stable?
    *  There must be two bases that satisfy the following conditions:
        *  There are more than min_cov number of bases counted
        *  There are more than min_freq percent of reads at the variant base
        *  The number of variant bases is more than the null model for that coverage
    *  If those are satisfied twice, the base with the highest count is stored as the "var_base",
       and the the bases with the second highest is the "ref_base"

    What is anySNP?
    *  For at least one mm position, there was a SNP called

    What is bases?
    *  A set of bases that are the var_base or ref_base at least once

    What are the things that are being missed?
    *  The case where theres only one base that meets the chriteria above, but its
       different from the reference base
    *  The case where >2 variant bases are present

    What is SNP table used for?
    *  Looking for specific variant positions; performing dn/ds analysis
    *  Establishing linkage patterns

    What to do about SNPs that are fixed opposite of reference?
    *  Create a "allele_count" column (formerly morphia); 2 means there are 2 bases; 1 means theres 1
        (and its opposite of reference); 3/4 means there multiple, though only
        two will be logged
        *  Right now only 2+ are being shown; so it will be easy to filter others out
    *  Throw in a check for if there's only one bases that fits that chriteria,
       is it not the reference base?
    *  Add "ref_base" back into the table to make it easier to compare SNP tables

    '''
    anySNP = False
    bases = set()
    ret_counts = None
    for mm in sorted(list(MMcounts.keys())):
        counts = inStrain.profile.profile_utilities.mm_counts_to_counts(MMcounts, mm)
        ret_counts = counts
        snp, morphia = call_snv_site(counts, ref_base, null_model, min_cov=min_cov, min_freq=min_freq) # Call SNP

        # Update clonality
        if mm not in clonT:
            clonT[mm] = np.full(mLen, fill_value=np.nan, dtype="float32")
        if sum(counts) >= min_cov:
            clonT[mm][pos] = calculate_clonality(counts)

        # Update rarefied clonality
        if mm not in clonTR:
            clonTR[mm] = np.full(mLen, fill_value=np.nan, dtype="float32")
        if sum(counts) >= min_covR:
            clonTR[mm][pos] = calculate_rarefied_clonality(counts, rarefied_coverage=min_covR)

        if not snp: # means base was not counted
            continue

        elif snp != -1: # means this is a SNP

           # calculate var_base
           counts_temp = list(counts)
           counts_temp[P2C[snp]] = 0
           var_base = C2P[list(counts_temp).index(sorted(counts_temp)[-1])] # this fixes the var_base = ref_base error when there's a tie - alexcc 5/8/2019

           SNP_class = calc_snp_class(snp, ref_base, var_base, counts,
                                        morphia, null_model, min_cov=min_cov,
                                        min_freq=min_freq)

           Stable['scaffold'].append(scaff)
           Stable['position'].append(pos)
           Stable['ref_base'].append(ref_base)
           for b, c in zip(['A', 'C', 'T', 'G'], counts):
               Stable[b].append(c)
           Stable['con_base'].append(snp)
           Stable['var_base'].append(var_base)
           Stable['mm'].append(mm)
           Stable['allele_count'].append(morphia)
           Stable['class'].append(SNP_class)

           if morphia >= 2:
               anySNP = True
               bases.add(snp)
               bases.add(var_base)
               ret_counts = counts

           elif (morphia == 1) & (anySNP == True):
               p2c[pos] = True

        # if it's now not a SNP, but it was in the past, mark it cryptic
        elif (snp == -1) & (anySNP == True):
            p2c[pos] = True

    if ret_counts is None:
        ret_counts = np.zeros(4, dtype=int)

    return anySNP, bases, ret_counts

def call_snv_site(counts, ref_base, model_to_use, min_cov=5, min_freq=0.05):
    '''
    Determines whether a site has a variant based on its nucleotide count frequencies.

    Returns one of the following and the "morphia" (number of bases present in the reads):
        Base if SNP
        -1 if not SNP
        None if unCounted

    A base is considered present in the reads if:
    *  There are more than min_cov number of bases counted
    *  There are more than min_freq percent of reads at the variant base
    *  The number of variant bases is more than the null model for that coverage

    The morphia is the number of bases present in the reads. Things are returned as SNPs if:
    *  Morphia >= 2
    *  Morphia == 1, and the base present is not the genome reference base
    *  Morphia == 0 (this is a special case meaning no bases are really present)
    '''
    # Make sure you have the coverage
    total = sum(counts)
    if total < min_cov:
        return None, 0

    # Count how many bases are there
    i = 0
    for c in counts:
        if total in model_to_use:
            min_bases = model_to_use[total]
        else:
            min_bases = model_to_use[-1]

        if c >= min_bases and float(c) / total >= min_freq:
            i += 1

    # If there are 2, say that
    if i > 1:
        return C2P[np.argmax(counts)], i

    # If theres only 1, see if its the reference base
    elif (i == 1) & (C2P[np.argmax(counts)] != ref_base):
        return C2P[np.argmax(counts)], i

    # That means this is a super polymorphic position with no dominant bases
    elif i == 0:
        return C2P[np.argmax(counts)], i

    # This means its just not a SNP; one dominant reference base
    else:
        return -1, i

def calc_snp_class(con_base, ref_base, var_base, counts, allele_count, null_model, min_cov=5, min_freq=0.05):
    '''
    Return the class of SNP this is
    '''
    if ref_base not in ['A', 'C', 'T', 'G']:
        return 'AmbiguousReference'

    if allele_count == 0:
        return 'DivergentSite'

    elif allele_count == 1:
        assert con_base != ref_base
        return 'SNS'

    else:
        if ref_base == con_base:
            return 'SNV'

        else:
            if ref_base == var_base:
                return 'con_SNV'

            if inStrain.readComparer.is_present(counts[P2C[ref_base]], sum(counts), null_model, float(min_freq)):
                return 'con_SNV'

            return 'pop_SNV'

def calculate_clonality(counts):
    '''
    Calculates the probability that two reads have the same allele at a position (per site nucleotide diversity)
    '''
    s = sum(counts)
    prob = (float(counts[0]) / s) * (float(counts[0]) / s) + (float(counts[1]) / s) * (float(counts[1]) / s) + (float(counts[2]) / s) * (float(counts[2]) / s) + (float(counts[3]) / s) * (float(counts[3]) / s)
    return prob

def calculate_rarefied_clonality(counts, rarefied_coverage=50):
    '''
    Calculates the probability that two reads have the same allele at a position (per site nucleotide diversity)
    '''
    # Figure out the probability distribution
    s = sum(counts)
    p = [counts[i]/s for i in [0,1,2,3]]

    # Make rarefied counts
    rcounts_list = np.random.choice([0,1,2,3], rarefied_coverage, p=p)
    unique, ncounts = np.unique(rcounts_list, return_counts=True)
    item2counts = dict(zip(unique, ncounts))
    rcounts = [item2counts[i] if i in item2counts else 0 for i in [0,1,2,3]]

    return calculate_clonality(rcounts)

def calc_snps(Odb, mm, null_model, min_freq=0.05):
    '''
    Calculate the number of reference SNPs, bi-allelic SNPs, multi-allelic SNPs (>2), total SNPs, consensus_SNPs, and poplation_SNPs

    Vocabulary:

    SNS_count = (number of substitutions) (allele_count = 1 and ref != con)
    SNV_count = number of SNVs (allele_count > 1)

    con_snps = con != ref and allele count > 0
    pop_snps = confusing
    '''
    if len(Odb) == 0:
        return [0, 0, 0, 0, 0]

    db = Odb[Odb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['position'], keep='last')

    SNS_count = len(db[(db['allele_count'] == 1)])
    SNV_count = len(db[db['allele_count'] > 1])

    con_snps = len(db[db['class'].isin(['SNS', 'con_SNV', 'pop_SNV'])])
    pop_snps = len(db[db['class'].isin(['SNS', 'pop_SNV'])])

    return [SNS_count, SNV_count, len(db), con_snps, pop_snps]

def generate_snp_table(Stable, scaffold, p2c):
    '''
    Make a dataframe and parse it a little bit
    '''
    SNPTable = pd.DataFrame(Stable)
    if len(SNPTable) > 0:
        SNPTable['scaffold'] = scaffold

        # Add cryptic SNPs
        SNPTable['cryptic'] = SNPTable['position'].map(p2c)
        SNPTable['cryptic'] = SNPTable['cryptic'].fillna(False)

        # Calc base coverage
        SNPTable['position_coverage'] = [sum([a,c,t,g]) for a,c,t,g in zip(SNPTable['A'],SNPTable['C'],SNPTable['T'],SNPTable['G'])]

    del Stable
    return SNPTable
