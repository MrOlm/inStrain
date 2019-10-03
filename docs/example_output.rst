Example Output
===================

InStrain produces a variety of output in the IS folder depending on which operations are run. Generally, output that is meant for human eyes to be easily interpretable is located in the ``output`` folder.

Typical output
-------

A typical run of inStrain will yield the following files in the output folder:

scaffold_info.tsv
+++++++++++++++++

This gives basic information about the scaffolds in your sample at the highest allowed level of read identity. Example below:

::
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
