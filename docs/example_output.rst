Example output
===================

InStrain produces a variety of output in the IS folder depending on which operations are run. Generally, output that is meant for human eyes to be easily interpretable is located in the ``output`` folder.

inStrain profile
---------------------

A typical run of inStrain will yield the following files in the output folder:

scaffold_info.tsv
+++++++++++++++++

This gives basic information about the scaffolds in your sample at the highest allowed level of read identity.

.. csv-table:: scaffold_info.tsv

    scaffold,length,coverage,breadth,nucl_diversity,coverage_median,coverage_std,coverage_SEM,breadth_minCov,breadth_expected,nucl_diversity_median,nucl_diversity_rarefied,nucl_diversity_rarefied_median,breadth_rarefied,conANI_reference,popANI_reference,SNS_count,SNV_count,divergent_site_count,consensus_divergent_sites,population_divergent_sites
    N5_271_010G1_scaffold_100,1148,1.89808362369338,0.9764808362369338,0.0,2,1.0372318863390368,0.030626273060932862,0.018292682926829267,0.8128805020451009,0.0,,,0.0,1.0,1.0,0,0,0,0,0
    N5_271_010G1_scaffold_102,1144,2.388986013986014,0.9956293706293706,0.003678160837326971,2,1.3042095721915248,0.038576628450898466,0.07604895104895107,0.8786983245100435,0.0,,,0.0,1.0,1.0,0,0,0,0,0
    N5_271_010G1_scaffold_101,1148,1.7439024390243902,0.9599303135888502,,2,0.8728918441975071,0.025773816178570358,0.0,0.7855901382035807,,,,0.0,0.0,0.0,0,00,0,0
    N5_271_010G1_scaffold_103,1142,2.039404553415061,0.9938704028021016,0.0,2,1.1288397384374758,0.03341869350286944,0.04028021015

scaffold
  The name of the :term:`scaffold` in the input .fasta file

length
  Full length of the :term:`scaffold` in the input .fasta file

coverage
  The average depth of coverage on the scaffold. If half the bases in a scaffold have 5 reads on them, and the other half have 10 reads, the coverage of the scaffold will be 7.5

breadth
  The percentage of bases in the scaffold that are covered by at least a single read. A breadth of 1 means that all bases in the scaffold have at least one read covering them

nucl_diversity
  The mean :term:`nucleotide diversity` of all bases in the scaffold that have a nucleotide diversity value calculated. So if only 1 base on the scaffold meats the minimum coverage to calculate nucleotide diversity, the nucl_diversity of the scaffold will be the nucleotide diversity of that base. Will be blank if no positions have a base over the minimum coverage.

coverage_median
  The median depth of coverage value of all bases in the scaffold, included bases with 0 coverage

coverage_std
  The standard deviation of all coverage values

coverage_SEM
  The standard error of the mean of all coverage values (calculated using `scipy.stats.sem <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sem.html>`_)

breadth_minCov
  The percentage of bases in the scaffold that have at least min_cov coverage (e.g. the percentage of bases that have a nucl_diversity value and meet the minimum sequencing depth to call SNVs)

breadth_expected
  This tells you the breadth that you should expect if reads are evenly distributed along the genome, given the reported coverage value. Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000. This is useful to establish whether or not the scaffold is actually in the reads, or just a fraction of the scaffold. If your coverage is 10x, the expected breadth will be ~1. If your actual breadth is significantly lower then the expected breadth, this means that reads are mapping only to a specific region of your scaffold (transposon, prophage, etc.)

nucl_diversity_median
  The median :term:`nucleotide diversity` value of all bases in the scaffold that have a :term:`nucleotide diversity` value calculated

nucl_diversity_rarefied
  The average :term:`nucleotide diversity` among positions that have at least ``--rarefied_coverage`` (50x by default). These values are also calculated by randomly subsetting the reads at that position to ``--rarefied_coverage`` reads

nucl_diversity_rarefied_median
  The median rarefied :term:`nucleotide diversity` (similar to that described above)

breadth_rarefied
  The percentage of bases in a scaffold that have at least ``--rarefied_coverage``

conANI_reference
  The :term:`conANI` between the reads and the reference genome

popANI_reference
    The :term:`popANI` between the reads and the reference genome

SNS_count
  The total number of :term:`SNSs<SNS>` called on this scaffold

SNV_count
  The total number of :term:`SNVs<SNV>` called on this scaffold

divergent_site_count
  The total number of :term:`divergent sites<divergent site>` called on this scaffold

consensus_divergent_sites
  The total number of :term:`divergent sites<divergent site>` in which the reads have a different consensus allele than the reference genome. These count as "differences" in the conANI_reference calculation, and ``breadth_minCov`` * ``length`` counts as the denominator.

population_divergent_sites
  The total number of :term:`divergent sites<divergent site>` in which the reads do not have the reference genome base as any allele at all (major or minor). These count as "differences" in the popANI_reference calculation, and ``breadth_minCov`` * ``length`` counts as the denominator.

mapping_info.tsv
+++++++++++++++++

This provides an overview of the number of reads that map to each scaffold, and some basic metrics about their quality. The header line (starting with #; not shown in the table below) describes the parameters that were used to filter the reads

.. csv-table:: mapping_info.tsv

    scaffold,pass_pairing_filter,filtered_pairs,unfiltered_priority_reads,filtered_priority_reads,pass_min_mapq,mean_insert_distance,median_insert,unfiltered_pairs,pass_min_read_ani,unfiltered_reads,mean_pair_length,mean_mapq_score,pass_max_insert,unfiltered_singletons,pass_min_insert,mean_PID,mean_mistmaches,filtered_singletons
    all_scaffolds,19293,7179,0,0,19293.0,307.724044990411,303.28290053387235,19293,7257.060551,253.4114963976572,15.254807443114085,19230.0,21965,19201.0,0.9388488368388499,15.244959311667445,0
    N5_271_010G1_scaffold_0,162,138,0,0,162.0,353.45061728395063,363.5,162,138.0,364,278.34567901234567,35.481481481481474,162.0,40,162.0,0.9829159042607164,4.697530864197532,0
    N5_271_010G1_scaffold_5,140,121,0,0,140.0,339.3142857142857,357.0,140,121.0,346,257.9214285714286,37.785714285714285,140.0,66,140.0,0.980420305410384,4.85,0

scaffold
  The name of the :term:`scaffold` in the input .fasta file. For the top row this will read ``all_scaffolds``, and it has the sum of all rows.

pass_pairing_filter
  The number of individual reads that pass the selecting pairing filter (only paired reads will pass this filter by default)

filtered_pairs
  The number of pairs of reads that pass all cutoffs

unfiltered_priority_reads
  The number of reads that pass the pairing filter because they were part of the ``priority_reads`` input file (will only be non-0 if a priority reads input file is provided).

filtered_priority_reads
  The number of priority reads that pass the rest of the filters (will only be non-0 if a priority reads input file is provided).

pass_min_mapq
  The number of pairs of reads mapping to this scaffold that pass the minimum mapQ score cutoff

mean_insert_distance
  Among all pairs of reads mapping to this scaffold, the mean insert distance. Note that the insert size is measured from the start of the first read to the end of the second read (2 perfectly overlapping 50bp reads will have an insert size of 50bp)

median_insert
  Among all pairs of reads mapping to this scaffold, the median insert distance.

unfiltered_pairs
  The raw number of pairs of reads that map to this scaffold. Only paired reads are used by inStrain

pass_min_read_ani
  The number of pairs of reads mapping to this scaffold that pass the min_read_ani cutoff

unfiltered_reads
  The raw number of reads that map to this scaffold

mean_pair_length
  Among all pairs of reads mapping to this scaffold, the average length of both reads in the pair summed together

mean_mapq_score
  Among all pairs of reads mapping to this scaffold, the average mapQ score

pass_max_insert
  The number of pairs of reads mapping to this scaffold that pass the maximum insert size cutoff- that is, their insert size is less than 3x the median insert size of all_scaffolds. Note that the insert size is measured from the start of the first read to the end of the second read (2 perfectly overlapping 50bp reads will have an insert size of 50bp)

unfiltered_singletons
  The number of reads detected in which only one read of the pair is mapped.

pass_min_insert
  The number of pairs of reads mapping to this scaffold that pass the minimum insert size cutoff

mean_PID
  Among all pairs of reads mapping to this scaffold, the average percentage ID of both reads in the pair to the reference .fasta file

mean_mistmaches
  Among all pairs of reads mapping to this scaffold, the mean number of mismatches

filtered_singletons
  The number of reads detected in which only one read of the pair is mapped AND which make it through to be considered. This will only be non-0 if the filtering settings allows non-paired reads.

SNVs.tsv
+++++++++++++++++

This describes the :term:`SNVs<SNV>` and :term:`SNSs<SNS>` that are detected in this mapping. While we should refer to these mutations as :term:`divergent sites<divergent site>`, sometimes SNV is used to refer to both SNVs and SNSs

.. csv-table:: SNVs.tsv

    scaffold,position,position_coverage,allele_count,ref_base,con_base,var_base,ref_freq,con_freq,var_freq,A,C,T,G,gene,mutation,mutation_type,cryptic,class
    N5_271_010G1_scaffold_120,174,5,2,C,C,A,0.6,0.6,0.4,2,3,0,0,,,I,False,SNV
    N5_271_010G1_scaffold_120,195,6,1,T,C,A,0.0,1.0,0.0,0,6,0,0,,,I,False,SNS
    N5_271_010G1_scaffold_120,411,8,2,A,A,C,0.75,0.75,0.25,6,2,0,0,N5_271_010G1_scaffold_120_1,N:V163G,N,False,SNV
    N5_271_010G1_scaffold_120,426,9,2,G,G,T,0.7777777777777778,0.7777777777777778,0.2222222222222222,0,0,2,7,N5_271_010G1_scaffold_120_1,N:S178Y,N,False,SNV
    N5_271_010G1_scaffold_120,481,6,2,C,T,C,0.3333333333333333,0.6666666666666666,0.3333333333333333,0,2,4,0,N5_271_010G1_scaffold_120_1,N:D233N,N,False,con_SNV
    N5_271_010G1_scaffold_120,484,6,2,G,A,G,0.3333333333333333,0.6666666666666666,0.3333333333333333,4,0,0,2,N5_271_010G1_scaffold_120_1,N:P236S,N,False,con_SNV
    N5_271_010G1_scaffold_120,488,5,1,T,C,T,0.2,0.8,0.2,0,4,1,0,N5_271_010G1_scaffold_120_1,S:240,S,False,SNS
    N5_271_010G1_scaffold_120,811,5,1,T,A,T,0.2,0.8,0.2,4,0,1,0,N5_271_010G1_scaffold_120_1,N:N563Y,N,False,SNS
    N5_271_010G1_scaffold_120,897,7,2,G,G,T,0.7142857142857143,0.7142857142857143,0.2857142857142857,0,0,2,5,,,I,False,SNV

See the :doc:`module_descriptions` for what constitutes a SNP (what makes it into this table)

scaffold
  The scaffold that the SNV is on

position
  The genomic position of the SNV

position_coverage
  The number of reads detected at this position

allele_count
  The number of bases that are detected above background levels (according to the :term:`null model`. An allele_count of 0 means no bases are supported by the reads, an allele_count of 1 means that only 1 base is supported by the reads, an allele_count of 2 means two bases are supported by the reads, etc.

ref_base
  The reference base in the .fasta file at that position

con_base
  The consensus base (the base that is supported by the most reads)

var_base
  Variant base; the base with the second most reads

ref_freq
  The fraction of reads supporting the ref_base

con_freq
  The fraction of reds supporting the con_base

var_freq
  The fraction of reads supporting the var_base

A, C, T, and G
  The number of mapped reads encoding each of the bases

gene
  If a gene file was included, this column will be present listing if the SNV is in the coding sequence of a gene

mutation
  Short-hand code for the amino acid switch. If synonymous, this will be S: + the position. If nonsynonymous, this will be N: + the old amino acid + the position + the new amino acid.

mutation_type
  What type of mutation this is. N = nonsynonymous, S = synonymous, I = intergenic, M = there are multiple genes with this base so you cant tell

cryptic
  If an SNV is cryptic, it means that it is detected when using a lower read mismatch threshold, but becomes undetected when you move to a higher read mismatch level. See "dealing with mm" in the advanced_use section for more details on what this means.

class
  The classification of this divergent site. The options are :term:`SNS` (meaning allele_count is 1 and con_base does not equal ref_base), :term:`SNV` (meaning allele_count is > 1 and con_base *does* equal ref_base), con_SNV (meaning allele_count is > 1, con_base does not equal ref_base, and ref_base *is* present in the reads; these count as differences in conANI calculations), pop_SNV (meaning allele_count is > 1, con_base does not equal ref_base, and ref_base *is not* present in the reads; these count as differences in popANI and conANI calculations), DivergentSite (meaning allele count is 0), and AmbiguousReference (meaning the ref_base is not A, C, T, or G)

linkage.tsv
+++++++++++++++++

This describes the :term:`linkage` between pairs of SNPs in the mapping that are found on the same read pair at least min_snp times.

.. csv-table:: linkage.tsv

    scaffold,position_A,position_B,distance,r2,d_prime,r2_normalized,d_prime_normalized,allele_A,allele_a,allele_B,allele_b,countab,countAb,countaB,countAB,total
    N5_271_010G1_scaffold_93,58,59,1,0.021739130434782702,1.0,0.031141868512110725,1.0,C,T,G,A,0,3,4,20,27
    N5_271_010G1_scaffold_93,58,70,12,0.012820512820512851,1.0,,,C,T,T,A,0,2,4,22,28
    N5_271_010G1_scaffold_93,58,80,22,0.016722408026755814,1.0,0.005847953216374271,1.0,C,T,G,A,0,2,5,21,28
    N5_271_010G1_scaffold_93,58,84,26,0.7652173913043475,1.0000000000000002,0.6296296296296297,1.0,C,T,G,C,4,0,1,22,27
    N5_271_010G1_scaffold_93,58,101,43,0.00907029478458067,1.0,,,C,T,C,A,0,2,2,19,23
    N5_271_010G1_scaffold_93,58,126,68,0.01754385964912257,1.0,0.002770083102493075,1.0,C,T,A,T,0,2,3,16,21
    N5_271_010G1_scaffold_93,58,133,75,0.008333333333333352,1.0,,,C,T,G,T,0,1,3,17,21
    N5_271_010G1_scaffold_93,59,70,11,0.010869565217391413,1.0,0.02777777777777779,1.0,G,A,T,A,0,2,3,21,26
    N5_271_010G1_scaffold_93,59,80,21,0.6410256410256397,1.0,1.0,1.0,G,A,G,A,2,0,1,25,28

Linkage is used primarily to determine if organisms are undergoing horizontal gene transfer or not. It's calculated for pairs of SNPs that can be connected by at least ``min_snp`` reads. It's based on the assumption that each SNP as two alleles (for example, a A and b B). This all gets a bit confusing and has a large amount of literature around each of these terms, but I'll do my best to briefly explain what's going on

scaffold
  The scaffold that both SNPs are on

position_A
  The position of the first SNP on this scaffold

position_B
  The position of the second SNP on this scaffold

distance
  The distance between the two SNPs

r2
  This is the r-squared linkage metric. See below for how it's calculated

d_prime
  This is the d-prime linkage metric. See below for how it's calculated

r2_normalized, d_prime_normalized
  These are calculated by rarefying to ``min_snp`` number of read pairs. See below for how it's calculated

allele_A
  One of the two bases at position_A

allele_a
  The other of the two bases at position_A

allele_B
  One of the bases at position_B

allele_b
  The other of the two bases at position_B

countab
  The number of read-pairs that have allele_a and allele_b

countAb
  The number of read-pairs that have allele_A and allele_b

countaB
  The number of read-pairs that have allele_a and allele_B

countAB
  The number of read-pairs that have allele_A and allele_B

total
  The total number of read-pairs that have have information for both position_A and position_B

Python code for the calculation of these metrics::

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

gene_info.tsv
+++++++++++++++++

This describes some basic information about the genes being profiled

.. csv-table:: gene_info.tsv

    scaffold,gene,gene_length,coverage,breadth,breadth_minCov,nucl_diversity,start,end,direction,partial,dNdS_substitutions,pNpS_variants,SNV_count,SNV_S_count,SNV_N_count,SNS_count,SNS_S_count,SNS_N_count,divergent_site_count
    N5_271_010G1_scaffold_0,N5_271_010G1_scaffold_0_1,141.0,0.7092198581560284,0.7092198581560284,0.0,,143,283,-1,False,,,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    N5_271_010G1_scaffold_0,N5_271_010G1_scaffold_0_2,219.0,4.849315068493151,1.0,0.45662100456620996,0.012312216758728069,2410,2628,-1,False,,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    N5_271_010G1_scaffold_0,N5_271_010G1_scaffold_0_3,282.0,7.528368794326241,1.0,0.9609929078014184,0.00805835530326815,3688,3969,-1,False,,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_1,336.0,2.7261904761904763,1.0,0.0625,0.0,0,335,-1,False,,,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_2,717.0,7.714086471408647,1.0,0.8926080892608089,0.011336830817162968,378,1094,-1,False,,0.554203539823008,9.0,2.0,6.0,0.0,0.0,0.0,9.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_3,114.0,13.105263157894735,1.0,1.0,0.016291986431991808,1051,1164,-1,False,,0.3956834532374099,4.0,1.0,2.0,0.0,0.0,0.0,4.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_4,111.0,11.342342342342342,1.0,1.0,0.02102806761458109,1164,1274,-1,False,,,5.0,0.0,5.0,0.0,0.0,0.0,5.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_5,174.0,9.057471264367816,1.0,1.0,0.006896087493019509,1476,1649,-1,False,,0.0,2.0,2.0,0.0,0.0,0.0,0.0,2.0
    N5_271_010G1_scaffold_1,N5_271_010G1_scaffold_1_6,174.0,6.195402298850576,1.0,0.7413793103448276,0.028698649055273976,1656,1829,-1,False,,0.5790697674418601,4.0,1.0,3.0,0.0,0.0,0.0,4.0

scaffold
  Scaffold that the gene is on

gene
  Name of the gene being profiled

gene_length
  Length of the gene in nucleotides

:term:`breadth`
  The number of bases in the gene that have at least 1x coverage

breadth_minCov
  The number of bases in the gene that have at least min_cov coverage

nucl_diversity
  The mean :term:`nucleotide diversity` of all bases in the gene that have a nucleotide diversity value calculated. So if only 1 base on the scaffold meats the minimum coverage to calculate nucleotide diversity, the nucl_diversity of the scaffold will be the nucleotide diversity of that base. Will be blank if no positions have a base over the minimum coverage.

start
  Start of the gene (position on scaffold; 0-indexed)

end
  End of the gene (position on scaffold; 0-indexed)

direction
  Direction of the gene (based on prodigal call). If -1, means the gene is not coded in the direction expressed by the .fasta file

partial
  If True this is a partial gene; based on not having `partial=00` in the record description provided by Prodigal

:term:`dNdS_substitutions<dN/dS>`
  The :term:`dN/dS` of :term:`SNSs<SNS>` detected in this gene. Will be blank if 0 N and/or 0 S substitutions are detected

:term:`pNpS_variants<pN/pS>`
  The :term:`pN/pS` of :term:`SNVs<SNV>` detected in this gene. Will be blank if 0 N and/or 0 S SNVs are detected

SNV_count
  Total number of :term:`SNVs<SNV>` detected in this gene

SNV_S_count
  Number of synonymous :term:`SNVs<SNV>` detected in this gene

SNV_N_count
  Number of non-synonymous :term:`SNVs<SNV>` detected in this gene

SNS_count
  Total number of :term:`SNSs<SNS>` detected in this gens

SNS_S_count
  Number of synonymous :term:`SNSs<SNS>` detected in this gens

SNS_N_count
  Number of non-synonymous :term:`SNSs<SNS>` detected in this gens

divergent_site_count
  Number of :term:`divergent sites<divergent site>` detected in this gens

genome_info.tsv
++++++++++++++++

Describes many of the above metrics on a genome-by-genome level, rather than a scaffold-by-scaffold level.

.. csv-table:: genome_info.tsv

    genome,coverage,breadth,nucl_diversity,length,true_scaffolds,detected_scaffolds,coverage_median,coverage_std,coverage_SEM,breadth_minCov,breadth_expected,nucl_diversity_rarefied,conANI_reference,popANI_reference,iRep,iRep_GC_corrected,linked_SNV_count,SNV_distance_mean,r2_mean,d_prime_mean,consensus_divergent_sites,population_divergent_sites,SNS_count,SNV_count,filtered_read_pair_count,reads_unfiltered_pairs,reads_mean_PID,reads_unfiltered_reads,divergent_site_count
    all_scaffolds,5.443497717712342,0.9638091828515172,0.012999760411488584,279325,178,178,4,13.042121218437627,0.02641794817118796,0.3605549091560011,0.9918244597989267,0.0017126380217433116,0.9968126936214156,0.999344665978235,,False,3000.0,92.27933333333333,0.1035315133374686,0.9402437184830787,321,66,63,1745,7179,19293,0.9802826482449184,60551,1808

genome
  The name of the genome being profiled. If all scaffolds were a single genome, this will read "all_scaffolds"

coverage
  Average :term:`coverage depth<coverage>` of all scaffolds of this genome

breadth
  The :term:`breadth` of all scaffolds of this genome

nucl_diversity
  The average :term:`nucleotide diversity` of all scaffolds of this genome

length
  The full length of this genome across all scaffolds

true_scaffolds
  The number of scaffolds present in this genome based off of the scaffold-to-bin file

detected_scaffolds
  The number of scaffolds with at least a single read-pair mapping to them

coverage_median
  The median :term:`coverage` amoung all bases in the genome

coverage_std
  The standard deviation of all coverage values

coverage_SEM
  The standard error of the mean of all coverage values (calculated using `scipy.stats.sem <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sem.html>`_)

breadth_minCov
  The percentage of bases in the scaffold that have at least min_cov coverage (e.g. the percentage of bases that have a nucl_diversity value and meet the minimum sequencing depth to call SNVs)

breadth_expected
  This tells you the breadth that you should expect if reads are evenly distributed along the genome, given the reported coverage value. Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000. This is useful to establish whether or not the scaffold is actually in the reads, or just a fraction of the scaffold. If your coverage is 10x, the expected breadth will be ~1. If your actual breadth is significantly lower then the expected breadth, this means that reads are mapping only to a specific region of your scaffold (transposon, prophage, etc.)

nucl_diversity_rarefied
  The average :term:`nucleotide diversity` among positions that have at least ``--rarefied_coverage`` (50x by default). These values are also calculated by randomly subsetting the reads at that position to ``--rarefied_coverage`` reads

conANI_reference
  The :term:`conANI` between the reads and the reference genome

popANI_reference
  The :term:`popANI` between the reads and the reference genome

iRep
  The :term:`iRep` value for this genome (if it could be successfully calculated)

iRep_GC_corrected
  A True / False value of whether the iRep value was corrected for GC bias

linked_SNV_count
  The number of :term:`divergent sites<divergent site>` that could be :term:`linked<linkage>` in this genome

SNV_distance_mean
  Average distance between linked :term:`divergent sites<divergent site>`

r2_mean
  Average r2 between linked SNVs (see explanation of linkage.tsv above for more info)

d_prime_mean
  Average d prime between linked SNVs (see explanation of linkage.tsv above for more info)

consensus_divergent_sites
  The total number of :term:`divergent sites<divergent site>` in which the reads have a different consensus allele than the reference genome. These count as "differences" in the conANI_reference calculation, and ``breadth_minCov`` * ``length`` counts as the denominator.

population_divergent_sites
  The total number of :term:`divergent sites<divergent site>` in which the reads do not have the reference genome base as any allele at all (major or minor). These count as "differences" in the popANI_reference calculation, and ``breadth_minCov`` * ``length`` counts as the denominator.

SNS_count
  The total number of :term:`SNSs<SNS>` called on this genome

SNV_count
  The total number of :term:`SNVs<SNV>` called on this genome

filtered_read_pair_count
  The total number of read pairs that pass filtering and map to this genome

reads_unfiltered_pairs
  The total number of pairs, filtered or unfiltered, that map to this genome

reads_mean_PID
  The average ANI of mapped read pairs to the reference genome for this genome

reads_unfiltered_reads
  The total number of reads, filtered or unfiltered, that map to this genome

divergent_site_count
  The total number of :term:`divergent sites<divergent site>` called on this genome

inStrain compare
---------------------

A typical run of inStrain will yield the following files in the output folder:

comparisonsTable.tsv
+++++++++++++++++++++

Summarizes the differences between two inStrain profiles on a scaffold by scaffold level

.. csv-table:: comparisonsTable.tsv

    scaffold,name1,name2,coverage_overlap,compared_bases_count,percent_genome_compared,length,consensus_SNPs,population_SNPs,popANI,conANI
    N5_271_010G1_scaffold_98,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,61,0.05290546400693842,1153,0,0,1.0,1.0
    N5_271_010G1_scaffold_133,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,78,0.0741444866920152,1052,0,0,1.0,1.0
    N5_271_010G1_scaffold_144,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,172,0.16715257531584066,1029,0,0,1.0,1.0
    N5_271_010G1_scaffold_158,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,36,0.035749751737835164,1007,0,0,1.0,1.0
    N5_271_010G1_scaffold_57,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,24,0.0183206106870229,1310,0,0,1.0,1.0
    N5_271_010G1_scaffold_139,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,24,0.023121387283236997,1038,0,0,1.0,1.0
    N5_271_010G1_scaffold_92,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,336,0.286934244235696,1171,0,0,1.0,1.0
    N5_271_010G1_scaffold_97,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,22,0.01901469317199654,1157,0,0,1.0,1.0
    N5_271_010G1_scaffold_100,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,21,0.018292682926829267,1148,0,0,1.0,1.0

scaffold
  The scaffold being compared

name1
  The name of the first `inStrain profile` being compared

name2
  The name of the second `inStrain profile` being compared

coverage_overlap
  The percentage of bases that are either covered or not covered in both of the profiles (covered = the base is present at at least min_snp coverage). The formula is length(coveredInBoth) / length(coveredInEither). If both scaffolds have 0 coverage, this will be 0.

compared_bases_count
  The number of considered bases; that is, the number of bases with at least min_snp coverage in both profiles. Formula is length([x for x in overlap if x == True]).

percent_genome_compared
  The percentage of bases in the scaffolds that are covered by both. The formula is length([x for x in overlap if x == True])/length(overlap). When ANI is np.nan, this must be 0. If both scaffolds have 0 coverage, this will be 0.

length
  The total length of the scaffold

consensus_SNPs
  The number of locations along the genome where both samples have the base at >= 5x coverage, and the consensus allele in each sample is different. Used to calculate :term:`conANI`

population_SNPs
  The number of locations along the genome where both samples have the base at >= 5x coverage, and no alleles are shared between either sample. Used to calculate :term:`popANI`

:term:`popANI`
  The average nucleotide identity among compared bases between the two scaffolds, based on population_SNPs. Calculated using the formula popANI = (compared_bases_count - population_SNPs) / compared_bases_count

:term:`conNI`
  The average nucleotide identity among compared bases between the two scaffolds, based on consensus_SNPs. Calculated using the formula conANI = (compared_bases_count - consensus_SNPs) / compared_bases_count

pairwise_SNP_locations.tsv
+++++++++++++++++++++++++++

Lists the locations of all differences between profiles

.. csv-table:: pairwise_SNP_locations.tsv

    scaffold,position,name1,name2,consensus_SNP,population_SNP,con_base_1,ref_base_1,var_base_1,position_coverage_1,A_1,C_1,T_1,G_1,con_base_2,ref_base_2,var_base_2,position_coverage_2,A_2,C_2,T_2,G_2
    N5_271_010G1_scaffold_9,823,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,G,G,A,10.0,3.0,0.0,0.0,7.0,A,G,G,6.0,3.0,0.0,0.0,3.0
    N5_271_010G1_scaffold_11,906,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,T,T,C,6.0,0.0,2.0,4.0,0.0,C,T,T,7.0,0.0,4.0,3.0,0.0
    N5_271_010G1_scaffold_29,436,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,C,T,T,6.0,0.0,3.0,3.0,0.0,T,T,C,7.0,0.0,3.0,4.0,0.0
    N5_271_010G1_scaffold_140,194,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,A,A,T,6.0,4.0,0.0,2.0,0.0,T,A,A,9.0,4.0,0.0,5.0,0.0
    N5_271_010G1_scaffold_24,1608,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,G,G,A,8.0,2.0,0.0,0.0,6.0,A,G,G,6.0,5.0,0.0,0.0,1.0
    N5_271_010G1_scaffold_112,600,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,A,G,G,6.0,4.0,0.0,0.0,2.0
    N5_271_010G1_scaffold_88,497,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,A,G,G,5.0,3.0,0.0,0.0,2.0
    N5_271_010G1_scaffold_53,1108,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,A,A,G,5.0,3.0,0.0,0.0,2.0,G,A,A,15.0,6.0,0.0,0.0,9.0
    N5_271_010G1_scaffold_46,710,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,True,False,A,C,C,6.0,4.0,2.0,0.0,0.0,C,C,A,6.0,2.0,4.0,0.0,0.0

scaffold
  The scaffold on which the difference is located

position
  The position where the difference is located (0-based)

name1
  The name of the first `inStrain profile` being compared

name2
  The name of the second `inStrain profile` being compared

consensus_SNP
  A True / False column listing whether or not this difference counts towards :term:`conANI` calculations

population_SNP
  A True / False column listing whether or not this difference counts towards :term:`popANI` calculations

con_base_1
  The consensus base of the profile listed in ``name1`` at this position

ref_base_1
  The reference base of the profile listed in ``name1`` at this position (will be the same as ``ref_base_2``)

var_base_1
  The variant base of the profile listed in ``name1`` at this position

position_coverage_1
  The number of reads mapping to this position in ``name1``

A_1, C_1, T_1, G_1
  The number of mapped reads with each nucleotide in ``name1``

con_base_2, ref_base_2, ...
  The above columns are also listed for the ``name2`` sample

genomeWide_compare.tsv
+++++++++++++++++++++++++++

A genome-level summary of the differences detected by inStrain compare. Generated by running ``inStrain genome_wide`` on the results of ``inStrain compare``

.. csv-table:: genomeWide_compare.tsv

    genome,name1,name2,coverage_overlap,compared_bases_count,consensus_SNPs,population_SNPs,popANI,conANI,percent_compared
    all_scaffolds,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,1.0,100712,0,0,1.0,1.0,0.3605549091560011
    all_scaffolds,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,0.6852932198159855,71900,196,,50.9999304589707928,0.9972739916550765,0.25740624720307886
    all_scaffolds,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam,1.0,145663,0,0,1.0,1.0,0.5214821444553835

genome
  The genome being compared

name1
  The name of the first `inStrain profile` being compared

name2
  The name of the second `inStrain profile` being compared

coverage_overlap
  The percentage of bases that are either covered or not covered in both of the profiles (covered = the base is present at at least min_snp coverage). The formula is length(coveredInBoth) / length(coveredInEither). If both scaffolds have 0 coverage, this will be 0.

compared_bases_count
  The number of considered bases; that is, the number of bases with at least min_snp coverage in both profiles. Formula is length([x for x in overlap if x == True]).

percent_genome_compared
  The percentage of bases in the scaffolds that are covered by both. The formula is length([x for x in overlap if x == True])/length(overlap). When ANI is np.nan, this must be 0. If both scaffolds have 0 coverage, this will be 0.

length
  The total length of the genome

consensus_SNPs
  The number of locations along the genome where both samples have the base at >= 5x coverage, and the consensus allele in each sample is different. Used to calculate :term:`conANI`

population_SNPs
  The number of locations along the genome where both samples have the base at >= 5x coverage, and no alleles are shared between either sample. Used to calculate :term:`popANI`

:term:`popANI`
  The average nucleotide identity among compared bases between the two scaffolds, based on population_SNPs. Calculated using the formula popANI = (compared_bases_count - population_SNPs) / compared_bases_count

:term:`conNI`
  The average nucleotide identity among compared bases between the two scaffolds, based on consensus_SNPs. Calculated using the formula conANI = (compared_bases_count - consensus_SNPs) / compared_bases_count

inStrain plot
----------------

This is what the results of inStrain plot look like.

1) Coverage and breadth vs. read mismatches
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/Example1.png
  :width: 800px
  :align: center

Breadth of coverage (blue line), coverage depth (red line), and expected breadth of coverage given the depth of coverage (dotted blue line) versus the minimum ANI of mapped reads. Coverage depth continues to increase while breadth of plateaus, suggesting that all regions of the reference genome are not present in the reads being mapped.

2) Genome-wide microdiversity metrics
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/genomeWide_microdiveristy_metrics_1.png
  :width: 800px
  :align: center

.. figure:: images/ExampleIS_plots/genomeWide_microdiveristy_metrics_2.png
  :width: 800px
  :align: center

SNV density, coverage, and nucleotide diversity. Spikes in nucleotide diversity and SNV density do not correspond with increased coverage, indicating that the signals are not due to read mis-mapping. Positions with nucleotide diversity and no SNV-density are those where diversity exists but is not high enough to call a SNV

3) Read-level ANI distribution
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/readANI_distribution.png
  :width: 800px
  :align: center

Distribution of read pair ANI levels when mapped to a reference genome; this plot suggests that the reference genome is >1% different than the mapped reads

4) Major allele frequencies
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/MajorAllele_frequency_plot.png
  :width: 800px
  :align: center

Distribution of the major allele frequencies of bi-allelic SNVs (the Site Frequency Spectrum). Alleles with major frequencies below 50% are the result of multiallelic sites. The lack of distinct puncta suggest that more than a few distinct strains are present.

5) Linkage decay
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/LinkageDecay_plot.png
  :width: 800px
  :align: center

.. figure:: images/ExampleIS_plots/Example5.png
  :width: 800px
  :align: center

Metrics of SNV linkage vs. distance between SNVs; linkage decay (shown in one plot and not the other) is a common signal of recombination.

6) Read filtering plots
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/ReadFiltering_plot.png
  :width: 800px
  :align: center

Bar plots showing how many reads got filtered out during filtering. All percentages are based on the number of paired reads; for an idea of how many reads were filtered out for being non-paired, compare the top bar and the second to top bar.

7) Scaffold inspection plot (large)
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/ScaffoldInspection_plot.png
  :width: 800px
  :align: center

This is an elongated version of the genome-wide microdiversity metrics that is long enough for you to read scaffold names on the y-axis

8) Linkage with SNP type (GENES REQUIRED)
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/LinkageDecay_types_plot.png
  :width: 800px
  :align: center

Linkage plot for pairs of non-synonymous SNPs and all pairs of SNPs

9) Gene histograms (GENES REQUIRED)
++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/GeneHistogram_plot.png
  :width: 800px
  :align: center

Histogram of values for all genes profiled

10) Compare dendrograms (RUN ON COMPARE; NOT PROFILE)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: images/ExampleIS_plots/Example10.png
  :width: 800px
  :align: center

A dendrogram comparing all samples based on popANI and based on shared_bases.
