Example Output
===================

InStrain produces a variety of output in the IS folder depending on which operations are run. Generally, output that is meant for human eyes to be easily interpretable is located in the ``output`` folder.

inStrain profile
-------

A typical run of inStrain will yield the following files in the output folder:

scaffold_info.tsv
+++++++++++++++++

This gives basic information about the scaffolds in your sample at the highest allowed level of read identity.

.. csv-table:: scaffold_info.tsv

  scaffold,length,breadth,coverage,median_cov,std_cov,bases_w_0_coverage,mean_clonality,median_clonality,mean_microdiversity,median_microdiversity,unmaskedBreadth,SNPs,expected_breadth,ANI
  1820A1025_scaffold_114,27197,1.0,117.53612530793836,119,22.517049757019283,0,0.9995305505412146,1.0,0.0004694494587853537,0.0,0.9996323123873956,5,1.0,0.999816088571744
  1820A1025_scaffold_106,28546,1.0,115.79373642541864,117,22.0981482314861,0,0.9995114336856696,1.0,0.0004885663143303631,0.0,0.9997898129335108,4,1.0,0.9998598458304134
  1820A1025_scaffold_112,27554,1.0,159.2970530594469,159,26.559269091036462,0,0.9994393799204232,1.0,0.0005606200795766902,0.0,0.9998911228859694,4,1.0,0.99985481470727
  1820A1025_scaffold_115,26897,1.0,171.2598059263115,167,36.6364045667205,0,0.999117155309884,1.0,0.0008828446901159025,0.0,1.0,46,1.0,0.998289772093542
  1820A1025_scaffold_118,25342,1.0,137.14375345276616,139,21.175189619940546,0,0.9994166977631538,1.0,0.0005833022368462171,0.0,0.9999605398153264,15,1.0,0.9994080738723808
  1820A1025_scaffold_130,22145,1.0,168.69356513885754,172,31.331335540083884,0,0.9993979713803196,1.0,0.0006020286196805058,0.0,1.0,8,1.0,0.9996387446376156
  1820A1025_scaffold_142,19159,1.0,163.43525236181432,164,26.75316027430802,0,0.9992356436742522,1.0,0.0007643563257477838,0.0,1.0,15,1.0,0.999217078135602
  1820A1025_scaffold_152,17330,1.0,143.43064050778997,144,26.235443720406167,0,0.9996060369840994,1.0,0.0003939630159005558,0.0,0.9978649740334681,0,1.0,1.0
  1820A1025_scaffold_158,16158,1.0,176.14370590419603,178,30.332807934861112,0,0.9993612181655604,1.0,0.0006387818344396612,0.0,0.9999381111523704,6,1.0,0.9996286439314228

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

std_cov
  The standard deviation of all coverage values

bases_w_0_coverage
  The number of bases with 0 coverage

mean_clonality
  The mean clonality value of all bases in the scaffold that have a clonality value calculated. So if only 1 base on the scaffold meats the minimum coverage to calculate clonality, the mean_clonality of the scaffold will be the clonality of that base

median_clonality
  The median clonality value of all bases in the scaffold that have a clonality value calculated

mean_microdiversity
  The mean mean_microdiversity value of all bases in the scaffold that have a mean_microdiversity value calculated (microdiveristy = 1 - clonality)

median_microdiversity
  The median microdiversity value of all bases in the scaffold that have a microdiversity value calculated

unmaskedBreadth
  The percentage of bases in the scaffold that have at least the min_cov number of bases. This value multiplied by the length of the scaffold gives the percentage of bases for which clonality is calculated and on which SNPs can be called

SNPs
  The total number of SNPs called on this scaffold

expected_breadth
  This tells you the breadth that you should expect if reads are evenly distributed along the genome, given the reported coverage value. Based on the function breadth = -1.000 * e^(0.883 * coverage) + 1.000. This is useful to establish whether or not the scaffold is actually in the reads, or just a fraction of the scaffold. If your coverage is 10x, the expected breadth will be ~1. If your actual breadth is significantly lower then the expected breadth, this means that reads are mapping only to a specific region of your scaffold (transposon, etc.)

ANI
  The average nucleotide identity between the reads in the sample and the .fasta file. Calculated using the formula ANI = (unmaskedBreadth * length) - SNPs)/ (unmaskedBreadth * length))

.. warning::

  As of inStrain v1.0.0, ANI includes consideration of all SNPs, including those that signify polymorphic positions and those that signify deviations from the reference base.

read_report.tsv
+++++++++++++++++

This provides an overview of the number of reads that map to each scaffold, and some basic metrics about their quality.

.. csv-table:: read_report.tsv

  scaffold,unfiltered_reads,unfiltered_pairs,pass_filter_cutoff,pass_max_insert,pass_min_insert,pass_min_mapq,filtered_pairs,mean_mistmaches,mean_insert_distance,mean_mapq_score,mean_pair_length,median_insert,mean_PID
  all_scaffolds,4464152,2168391,2154253,2168102,2167792,2168391,2153442,0.5102327024969205,325.4543267335089,41.506261555226885,293.22514574170435,313.0,0.9981314852335254
  JBBB007E_scaffold_233,10605,5062,5048,5062,5062,5062,5048,0.3832477281706835,312.3638877913868,1.3024496246542872,293.6845120505729,308.0,0.998581261373412
  1820A1025_scaffold_11,200651,98185,97870,98174,98150,98185,97825,0.4297295920965524,320.3908336303916,41.902917960992006,293.49184702347617,314.0,0.9984105714013544
  1820A1025_scaffold_12,192161,94169,93860,94165,94146,94169,93835,0.4243647060072848,320.1241916129512,41.86061230341195,293.2882264864233,313.0,0.9984318946277784
  1820A1025_scaffold_5,328711,160303,159657,160296,160265,160303,159613,0.4424371346762069,320.44356624642086,41.882241754677075,293.1935023050099,313.0,0.9983611885767244
  1820A1025_scaffold_2,428162,209675,208987,209651,209621,209675,208913,0.4234267318469059,328.2336473113151,41.89234768093478,293.3190318349827,313.0,0.9984327943894524
  1820A1025_scaffold_22,180103,87912,87523,87908,87880,87912,87488,0.4677632177632178,320.15531440531436,41.63107425607426,292.8853512603513,314.0,0.9982475272431892
  1820A1025_scaffold_1,461964,226192,225400,226155,226126,226192,225299,0.4299002617245526,336.90442190705244,41.908409669661175,293.1702359057792,313.0,0.9984114970717194

The following metrics are provided for all individual scaffolds, and for all scaffolds together (scaffold "all_scaffolds"). For the max insert cutoff, the median_insert for all_scaffolds is used

header line
  The header line (starting with #; not shown in the above table) describes the parameters that were used to filter the reads

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

This describes the SNPs that are detected in this mapping.

.. csv-table:: SNVs.tsv

  scaffold,position,refBase,A,C,T,G,conBase,varBase,allele_count,cryptic,baseCoverage,varFreqrefFreq
  1820A1025_scaffold_114,27192,A,13,2,0,0,A,C,2,False,15,0.13333333333333333,0.8666666666666667
  1820A1025_scaffold_114,27193,T,8,2,4,1,A,T,3,False,15,0.26666666666666666,0.5333333333333333
  1820A1025_scaffold_114,27194,A,11,4,0,0,A,C,2,False,15,0.26666666666666666,0.7333333333333333
  1820A1025_scaffold_114,27195,G,1,10,4,0,C,T,2,False,15,0.26666666666666666,0.6666666666666666
  1820A1025_scaffold_114,27196,A,0,10,0,0,C,A,1,False,10,0.0,1.0
  1820A1025_scaffold_106,11174,C,4,57,0,0,C,A,2,False,61,0.06557377049180327,0.9344262295081968
  1820A1025_scaffold_106,28450,A,37,2,0,0,A,C,2,True,39,0.05128205128205128,0.9487179487179488
  1820A1025_scaffold_106,28541,A,3,2,2,1,A,C,3,False,8,0.25,0.375
  1820A1025_scaffold_106,28542,G,0,4,0,4,C,G,2,False,8,0.5,0.5

See the :doc:`module_descriptions` for what constitutes a SNP (what makes it into this table)

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

This describes the linkage between pairs of SNPs in the mapping that are found on the same read pair at least min_snp times.

.. csv-table:: linkage.tsv

  r2,d_prime,r2_normalized,d_prime_normalized,total,countAB,countAb,countaB,countab,allele_A,allele_a,allele_B,allele_b,distance,position_A,position_B,scaffold
  1.000000000000002,1.0,0.9999999999999998,1.0,162,152,0,0,10,G,A,C,A,24,24679,24703,1820A1025_scaffold_115
  0.9999999999999986,1.0,0.9999999999999998,1.0,122,113,0,0,9,G,A,C,A,58,24679,24737,1820A1025_scaffold_115
  1.0000000000000004,1.0,126,118,0,0,8,G,A,C,G,59,24679,24738,1820A1025_scaffold_115
  0.9999999999999984,1.0,1.0,1.0,125,117,0,0,8,G,A,T,G,70,24679,24749,1820A1025_scaffold_115
  1.0,1.0,1.0,1.0,48,42,0,0,6,G,A,T,C,160,24679,24839,1820A1025_scaffold_115
  1.0000000000000016,1.0,0.9999999999999998,1.0,49,43,0,0,6,G,A,G,A,163,24679,24842,1820A1025_scaffold_115
  0.9999999999999998,1.0,0.9999999999999998,0.9999999999999998,46,40,0,0,6,G,A,GA,169,24679,24848,1820A1025_scaffold_115
  1.0000000000000009,1.0000000000000002,1.0,1.0,49,42,0,0,7,G,A,C,T,181,24679,24860,1820A1025_scaffold_115
  1.0,1.0,1.0,1.0,43,38,0,0,5,G,A,C,G,238,24679,24917,1820A1025_scaffold_115

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

inStrain compare
-------

A typical run of inStrain will yield the following files in the output folder:

inStrain profile_genes
-----------

A typical run of inStrain profile_genes will yield the following additional files in the output folder:

gene_info.tsv
+++++++++++

This describes some basic information about the genes being profiled

.. csv-table:: gene_info.tsv

  gene,scaffold,direction,partial,start,end,coverage,breadth,clonality,microdiversity,masked_breadthSNPs_per_bp,min_ANI
  JBBB007E_scaffold_233_1,JBBB007E_scaffold_233,-1,True,1,327,31.07645259938838,0.9908256880733946,0.9996904894417408,0.00030951055825922946,0.9357798165137616,0.0,0
  JBBB007E_scaffold_233_2,JBBB007E_scaffold_233,-1,False,344,1654,68.31273836765827,1.0,0.9995749550953246,0.0004250449046754312,1.0,0.0,0
  JBBB007E_scaffold_233_3,JBBB007E_scaffold_233,-1,False,1655,1840,57.55913978494624,1.0,0.9995984579286268,0.0004015420713732176,1.0,0.0,0
  JBBB007E_scaffold_233_4,JBBB007E_scaffold_233,-1,False,1824,2261,60.378995433789946,1.0,0.999627013456876,0.0003729865431241208,1.0,0.0,0
  JBBB007E_scaffold_233_5,JBBB007E_scaffold_233,-1,False,2254,2679,69.35915492957747,1.0,0.9996632362755252,0.00033676372447488667,1.0,0.0,0
  JBBB007E_scaffold_233_6,JBBB007E_scaffold_233,-1,False,2679,3026,49.05172413793103,1.0,0.9995279442304852,0.000472055769514812,1.0,0.0,0
  JBBB007E_scaffold_233_7,JBBB007E_scaffold_233,-1,False,3020,3421,53.14427860696517,1.0,0.9997191567029526,0.00028084329704736183,1.0,0.0,0
  JBBB007E_scaffold_233_8,JBBB007E_scaffold_233,-1,False,3431,3664,65.02991452991454,1.0,0.999618145899895,0.00038185410010505016,1.0,0.0,0
  JBBB007E_scaffold_233_9,JBBB007E_scaffold_233,-1,False,3664,4563,64.00333333333333,1.0,0.9994913809167014,0.0005086190832985782,1.0,0.0,0

SNP_mutation_types.tsv
+++++++++++++++

This describes whether SNPs are synonymous, nonsynonymous, or intergenic

.. csv-table:: SNP_mutation_types.tsv

  scaffold,position,refBase,A,C,T,G,conBase,varBase,baseCoverage,varFreq,refFreq,mutation_type,mutation,gene
  1820A1025_scaffold_5,96846,A,142,0,0,17,A,G,159,0.1069182389937107,0.8930817610062893,S,S:296,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96849,T,17,0,142,0,T,A,159,0.1069182389937107,0.8930817610062893,S,S:299,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96866,A,142,0,0,16,A,G,158,0.10126582278481013,0.8987341772151899,N,N:D316G,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96867,T,16,0,139,0,T,A,155,0.1032258064516129,0.8967741935483872,N,N:D317E,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96868,T,16,0,136,0,T,A,152,0.10526315789473684,0.8947368421052632,N,N:S318T,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96869,C,0,142,0,15,C,G,157,0.09554140127388536,0.9044585987261148,N,N:S319C,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96879,T,0,16,143,0,T,C,159,0.10062893081761007,0.89937106918239,S,S:329,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,96918,G,0,6,0,113,G,C,119,0.050420168067226885,0.9495798319327732,S,S:368,1820A1025_scaffold_5_96
  1820A1025_scaffold_5,345729,A,3,0,2,0,A,T,5,0.4,0.6,N,N:I756L,1820A1025_scaffold_5_314
  1820A1025_scaffold_5,345730,T,0,0,3,2,T,G,5,0.4,0.6,N,N:I757R,1820A1025_scaffold_5_314

inStrain genome_wide
------------

A typical run of inStrain genome_wide will yield the following additional files in the output folder:

genomeWide_scaffold_info.tsv
+++++++++++++

This is a genome-wide version of the scaffold report described above

.. csv-table:: genomeWide_scaffold_info.tsv

  genome,detected_scaffolds,true_scaffolds,true_length,SNPs,breadth,coverage,std_cov,mean_clonality,ANI,unmaskedBreadth,expected_breadth
  S2_002_005G1_phage_Clostridioides_difficile.fasta,1,1,21096,1,0.9993837694349641,61.11220136518772,11.758200608087787,0.9996280813990492,0.9999525098542053,0.9981513083048921,1.0
  S2_018_020G1_bacteria_Clostridioides_difficile.fasta,34,34,4075786,1235,0.9978202977290761,134.79707030742046,25.46058779187456,0.9994819824195843,0.9996962394194163,0.9975258759905451,1.0

genomeWide_read_report.tsv
++++++++++++

This is a genome-wide version of the read report described above

.. csv-table:: genomeWide_read_report.tsv

  genome,scaffolds,unfiltered_reads,unfiltered_pairs,pass_filter_cutoff,pass_max_insert,pass_min_insert,pass_min_mapq,filtered_pairs,mean_mistmaches,mean_insert_distance,mean_mapq_score,mean_pair_length,median_insert,mean_PID
  S2_002_005G1_phage_Clostridioides_difficile.fasta,1,10605,5062,5048,5062,5062,5062,5048,0.3832477281706835,312.3638877913868,1.3024496246542872,293.6845120505729,308.0,0.998581261373412
  S2_018_020G1_bacteria_Clostridioides_difficile.fasta,34,4453547,2163329,2149205,2163040,2162730,2163329,2148394,0.5636466689761853,321.3510672021471,41.47419579138972,293.33494491093336,312.5147058823529,0.9979527547934701

inStrain plot
------------

This is what those plots look like
