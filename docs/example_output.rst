Example Output
===================

InStrain produces a variety of output in the IS folder depending on which operations are run. Generally, output that is meant for human eyes to be easily interpretable is located in the ``output`` folder.

inStrain profile output
-------

A typical run of inStrain will yield the following files in the output folder:

scaffold_info.tsv
+++++++++++++++++

This gives basic information about the scaffolds in your sample at the highest allowed level of read identity. Example below ::

  scaffold        length  breadth coverage        median_cov      std_cov bases_w_0_coverage      mean_clonality  median_clonality        unmaskedBreadth SNPs    expected_breadth    ANI
  1608A0826_scaffold_3621 5866    0.18615751789976132     0.27395158540743264     0       0.6119560844337848      4774    0.0     0.0     0.0     0       0.21486472398684486     0.0
  1608A0826_scaffold_3624 2829    0.30505478967833155     0.3336868151290209      0       0.5287759217659689      1966    0.0     0.0     0.0     0       0.25520439783710397     0.0
  1608A0826_scaffold_3631 3282    0.22516758074344911     0.32693479585618523     0       0.6823711580458389      2543    0.0     0.0     0.0     0       0.2507506383274967      0.0
  1608A0826_scaffold_3633 1893    0.6724775488642366      1.5884838880084522      1       1.664144253208806       620     1.0     1.0     0.104067617538299       0       0.7540510470980886   1.0
  1608A0826_scaffold_3642 2661    0.5606914693724164      1.4753851935362643      1       1.936771851529394       1169    1.0     1.0     0.10973318301390454     0       0.7282207507647787   1.0
  1608A0826_scaffold_3645 1887    0.4843667196608373      0.7456279809220986      0       0.9161743431517859      973     0.0     0.0     0.0     0       0.48231560949644997     0.0
  1608A0826_scaffold_3646 2168    0.6342250922509225      1.084870848708487       1       1.0875816472378086      793     0.0     0.0     0.0     0       0.6163179098079582      0.0
  1608A0826_scaffold_3647 1887    0.255431902490726       0.255431902490726       0       0.4361037097763493      1405    0.0     0.0     0.0     0       0.20191994730181093     0.0
  1608A0826_scaffold_3652 1884    0.7547770700636943      2.004777070063694       2       1.5722104009266933      462     0.9955963301002432      1.0     0.057855626326963915    1   0.8297041657514823       0.9908256880733946

scaffold
  The name of the scaffold in the input .fasta file

length
  Full length of the scaffold in the input .fasta file

breadth
  The percentage of bases in the scaffold that are covered by at least a single read. A breadth of 1 means that all bases in the scaffold have at least one read covering them

coverage
  The average depth of coverage on the scaffold. If half the bases in a scaffold have 5 reads on them, and the other half have 10 reads, the coverage of the scaffold will be 7.5

median_cov
  The median coverage value of all bases in the scaffold, included bases with 0 coverage

bases_w_0_coverage
  The number of bases with 0 coverage

mean_clonality
  The mean clonality value of all bases in the scaffold that have a clonality value calculated. So if only 1 base on the scaffold meats the minimum coverage to calculate clonality, the mean_clonality of the scaffold will be the clonality of that base

median_clonality
  The median clonality value of all bases in the scaffold that have a clonality value calculated

unmaskedBreadth
  The percentage of bases in the scaffold that have at least the min_cov number of bases. This value multiplied by the length of the scaffold gives the percentage of bases for which clonality is calculated and on which SNPs can be called

SNPs
  The number of SNPs called on this scaffold

expected_breadth
  This tells you the breadth that you should expect if reads are evenly distributed along the genome, given the reported coverage value. Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000. This is useful to establish whether or not the scaffold is actually in the reads, or just a fraction of the scaffold. If your coverage is 10x, the expected breadth will be ~1. If your actual breadth is significantly lower, this means that reads are mapping only to a specific region of your scaffold (transposon, etc.)

ANI
  The average nucleotide identity between the reads in the sample and the .fasta file. Calculated using the formula ANI = (unmaskedBreadth * length) - SNPs)/ (unmaskedBreadth * length))

read_report.tsv
+++++++++++++++++

This provides an overview of the number of reads that map to each scaffold, and some basic metrics about their quality. Example below ::

  # filter_cutoff:0.95 max_insert_relative:3 min_insert:50 min_mapq:-1
  scaffold        unfiltered_reads        unfiltered_pairs        pass_filter_cutoff      pass_max_insert pass_min_insert pass_min_mapq   filtered_pairs  mean_mistmaches mean_insert_distance mean_mapq_score mean_pair_length        median_insert   mean_PID
  all_scaffolds   205082276       96617528        94929495        96590231        96573054        96617524        94862959        1.9608805247014809      335.7611296886006       15.170048420199695   286.27552505793767      329.0   0.9930197361098816
  b003-d052_scaffold_330  326     101     96      101     99      101     94      2.118811881188119       194.5049504950495       1.0693069306930694      251.3465346534653       210.00.9912149801106098
  b003-d052_scaffold_468  0       0       0       0       0       0       0
  b003-d052_scaffold_140  0       0       0       0       0       0       0
  b003-d086_scaffold_1    39      18      18      18      18      18      18      1.8333333333333333      326.72222222222223      1.8888888888888888      271.94444444444446      341.00.9925481073725309
  b003-d052_scaffold_687  0       0       0       0       0       0       0
  b003-d052_scaffold_64   672     225     210     225     223     225     208     2.7644444444444445      219.26222222222225      2.0311111111111106      274.3688888888889       212.00.9901429258972524
  b003-d052_scaffold_312  400     177     176     177     177     177     176     4.59322033898305        281.090395480226        5.175141242937853       282.0169491525424       287.00.9839838132309756

The following metrics are provided for all individual scaffolds, and for all scaffolds together (scaffold "all_scaffolds"). For the max insert cutoff, the median_insert for all_scaffolds is used

header line
  The header line (starting with #) describes the parameters that were used to filter the reads

scaffold
  The name of the scaffold in the input .fasta file

unfiltered_reads
  The raw number of reads that map to this scaffold

unfiltered_pairs
  The raw number of pairs of reads that map to this scaffold. Only paired reads are used by inStrain

pass_filter_cutoff
  The number of pairs of reads mapping to this scaffold that pass the ANI filter cutoff (specified in the header as "filter_cutoff")

pass_max_insert
  The number of pairs of reads mapping to this scaffold that pass the maximum insert size cutoff- that is, their insert size is less than 3x the median insert size of all_scaffolds. Note that the insert size is measured from the start of the first read to the end of the second read (2 perfectly overlapping 50bp reads will have an insert size of 50bp)

pass_min_insert
  The number of pairs of reads mapping to this scaffold that pass the minimum insert size cutoff

pass_min_mapq
  The number of pairs of reads mapping to this scaffold that pass the minimum mapQ score cutoff

filtered_pairs
  The number of pairs of reads that pass all cutoffs

mean_mistmaches
  Among all pairs of reads mapping to this scaffold, the mean number of mismatches

mean_insert_distance
  Among all pairs of reads mapping to this scaffold, the mean insert distance. Note that the insert size is measured from the start of the first read to the end of the second read (2 perfectly overlapping 50bp reads will have an insert size of 50bp)

mean_mapq_score
  Among all pairs of reads mapping to this scaffold, the average mapQ score

mean_pair_length
  Among all pairs of reads mapping to this scaffold, the average length of both reads in the pair summed together

median_insert
  Among all pairs of reads mapping to this scaffold, the median insert distance.

mean_PID
  Among all pairs of reads mapping to this scaffold, the average percentage ID of both reads in the pair to the reference .fasta file

SNVs.tsv
+++++++++++++++++

This describes the SNPs that are detected in this mapping. Example below ::

  scaffold        position        refBase A       C       T       G       conBase varBase morphia cryptic baseCoverage    varFreq refFreq
  1608A0826_scaffold_3652 941     G       2       0       0       3       G       A       0       False   5       0.4     0.6
  1723A1010_scaffold_464  3543    C       0       3       2       0       C       T       0       False   5       0.4     0.6
  1723A1010_scaffold_464  3576    G       4       0       0       1       A       G       1       False   5       0.2     0.8
  1723A1010_scaffold_472  418     T       0       5       0       0       C       A       1       False   5       0.0     1.0
  1723A1010_scaffold_472  765     G       5       0       0       0       A       A       1       False   5       1.0     1.0
  1723A1010_scaffold_472  1269    C       0       0       6       0       T       A       1       False   6       0.0     1.0
  1723A1010_scaffold_472  1287    C       0       0       8       0       T       A       1       False   8       0.0     1.0
  1723A1010_scaffold_472  2034    A       1       0       0       4       G       A       1       False   5       0.2     0.8
  1723A1010_scaffold_476  1139    A       0       0       0       5       G       A       1       False   5       0.0     1.0

See the section "module_descriptions" for what constitutes a SNP (what makes it into this table)

scaffold
  The scaffold that the SNP is on

position
  The genomic position of the SNP

refBase
  The reference base in the .fasta file at that position

A, C, T, and G
  The number of mapped reads encoding each of the bases

conBase
  The consensus base; the base that is supported by the most reads

varBase
  Variant base; the base with the second most reads

morphia
  The number of bases that are detected above background levels. In order to be detected above background levels, you must pass an fdr filter. See module descriptions for a description of how that works. A morphia of 0 means no bases are supported by the reads, a morphia of 1 means that only 1 base is supported by the reads, a morphia of 2 means two bases are supported by the reads, etc.

cryptic
  If a SNP is cryptic, it means that it is detected when using a lower read mismatch threshold, but becomes undetected when you move to a higher read mismatch level. See "dealing with mm" in the advanced_use section for more details on what this means.

baseCoverage
  The total number of reads at this position

varFreq
  The fraction of reads supporting the varBase

refFreq
  The fraction of reds supporting the refBase

linkage.tsv
+++++++++++++++++

This describes the linkage between pairs of SNPs in the mapping that are found on the same read pair at least min_snp times. Example below ::

  r2      d_prime r2_normalized   d_prime_normalized      total   countAB countAb countaB countab allele_A        allele_a        allele_B        allele_b        distance        position_A   position_B      scaffold
  1.0000000000000009      1.0     0.9999999999999998      1.0     80      68      0       0       12      G       A       C       A       8       187     195     1727A1014_scaffold_559
  0.9999999999999992      1.0     1.0     1.0     78      65      0       0       13      G       A       T       C       9       187     196     1727A1014_scaffold_559
  0.9999999999999996      1.0     0.9999999999999998      0.9999999999999998      67      54      0       0       13      G       A       A       G       21      187     208     1727A1014_scaffold_559
  0.9999999999999988      1.0     0.9999999999999998      0.9999999999999998      59      48      0       0       11      G       A       G       T       22      187     209     1727A1014_scaffold_559
  1.0     1.0     0.9999999999999998      0.9999999999999998      43      30      0       0       13      G       A       T       A       36      187     223     1727A1014_scaffold_559
  1.0000000000000002      1.0     1.0     1.0     39      27      0       0       12      G       A       C       T       46      187     233     1727A1014_scaffold_559
  1.0000000000000002      1.0000000000000002      0.9999999999999998      1.0     27      14      0       0       13      G       A       C       A       61      187     248     1727A1014_scaffold_559
  0.9999999999999998      1.0     0.9999999999999998      1.0     78      66      0       0       12      C       A       T       C       1       195     196     1727A1014_scaffold_559
  1.0     1.0                     68      56      0       0       12      C       A       A       G       13      195     208     1727A1014_scaffold_559

Linkage is used primarily to determine if organisms are undergoing horizontal gene transfer or not. It's calculated for pairs of SNPs that can be connected by at least ``min_snp`` reads. It's based on the assumption that each SNP as two alleles (for example, a A and b B). This all gets a bit confusing and has a large amount of literature around each of these terms, but I'll do my best to briefly explain what's going on

scaffold
  The scaffold that both SNPs are on

position_A
  The position of the first SNP on this scaffold

position_B
  The position of the second SNP on this scaffold

distance
  The distance between the two SNPs

allele_A
  One of the two bases at position_A

allele_a
  The other of the two bases at position_A

allele_B
  One of the bases at position_B

allele_b
  The other of the two bases at position_B

countAB
  The number of read-pairs that have allele_A and allele_B

countAb
  The number of read-pairs that have allele_A and allele_b

countaB
  The number of read-pairs that have allele_a and allele_B

countab
  The number of read-pairs that have allele_a and allele_b

total
  The total number of read-pairs that have have information for both position_A and position_B

r2
  This is the r-squared linkage metric. See below for how it's calculated

d_prime
  This is the d-prime linkage metric. See below for how it's calculated

r2_normalized, d_prime_normalized
  These are calculated by rarefying to ``min_snp`` number of read pairs. See below for how it's calculated

Code for the calculation of these metrics::

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

inStrain compare output
-------

A typical run of inStrain will yield the following files in the output folder:
