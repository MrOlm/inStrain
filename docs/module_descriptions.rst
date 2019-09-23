Module Descriptions
===================

The functionality of inStrain is broken up into modules. To see a list of available modules, check the help::

 $ inStrain -h

               ...::: inStrain v0.8.0 :::...

 Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

 Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
 Example: inStrain profile -h

   profile         -> Calculate microdiversity metrics and SNVs from a mapping. Must run this first to perform most other operations
   compare         -> Compare multiple inStrain profiles on a microdiversity level. Calculates popANI, coverage_overlap, and other things
   profile_genes   -> Calculate gene-level metrics on an inStrain profile
   quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
   filter_reads    -> Commands related to filtering reads from .bam files
   other           -> Other miscellaneous operations

IS_profile
--------------

An IS_profile (inStrain profile) is created by running the `inStrain profile` command. It contains  all of the program's internal workings, cached data, and output is stored. Additional modules can then be run on an IS_profile (to analyze genes, compare profiles, etc.), and there is lots of nice cached data stored in it that can be accessed using python.

.. seealso::

  :doc:`example_output`
    for help finding where the output from your run is located in the IS_profile

  :doc:`advanced_use`
    for access to the raw internal data (which can be very useful)

profile
------

The most complex part of inStrain, profile has several steps.

First, all reads in the .bam file are filtered to only keep those that map with sufficient quality. Reads must be paired (all non-paired reads will be filtered) and in all cases filters are applied on the pair, not the individual read. Command line parameters can be adjusted to change the specifics, but in general:

 * Pairs must be mapped in the proper orientation. The minimum insert distance can be set with a command line parameter. The maximum insert distance is a multiple of the median insert distance. So if pairs have a median insert size of 500bp, by default all pairs with insert sizes over 1500bp will be excluded.

 * Pairs must have a minimum mapQ score. In this case, the read in the pair with the higher mapQ is used. MapQ scores are confusing and how they're calculated varies based on the mapping algorithm being used, but are meant to represent both the number of mismatches in the mapping and how unique that mapping is. With bowtie2, if the read maps equally well to two positions on the genome, its mapQ score will be set to 2.

 * Pairs must be above some minimum nucleotide identity (ANI) value. For example if reads in a pair are 100bp each, and each read has a single mismatch, the ANI of that pair would be 0.99

Next, using only read pairs that pass filters, a number of microdiveristy metrics are calculated on a scaffold-by-scaffold basis. This includes:

 * Calculate the coverage at each position along the scaffold

 * Calculate the clonality at each position along the scaffold in which the coverage is greater than the min_cov argument. The formula for calculating clonality is the sum of the frequency of each base squared - [(frequency of A)^2 + (frequency of C)^2 + (frequency of G)^2 + (frequency of T)^2 ]. This clonality definition is nice because it is not effected by coverage

 * Identify SNPs. The criteria for being called a SNP are 1) More than min_cov number of bases at that position, 2) More than min_freq percentage of reads that are a variant base, 3) The number of reads with the variant base is more than the null model for that coverage. The null model describes the probability that the number of true reads that support a variant base could be due to random mutation error, assuming Q30 score. The default false discovery rate with the null model is 1e-6 (one in a million)

 * Calculate linkage between SNPs on the same read pair. For each pair harboring a SNP, calculate the linkage of that SNP with other SNPs within that same pair. This is only done for pairs of SNPs that are both on at least MIN_SNP reads

 * Calculate scaffold-level properties. These include things like the overall coverage, breadth of coverage, average nucleotide identity (ANI) between the reads and the reference genome, and the expected breadth of coverage based on that true coverage.

Finally, this information is stored as an IS_profile object. This includes the locations of SNPs, the number of read pairs that passed filters (and other information) for each scaffold, the linkage between SNV pairs, ect.

.. seealso::

  :doc:`example_output`
    for help interpreting the output

  :doc:`advanced_use`
    for access to the raw internal data (which can be very useful)

  :doc:`choosing_parameters`
    for information about the pitfalls and other things to consider when running inStrain

To see the command-line options, check the help::

  $ inStrain profile -h
  usage: inStrain profile [-o OUTPUT] [-p PROCESSES] [-d] [-h]
                         [-l FILTER_CUTOFF] [--min_mapq MIN_MAPQ]
                         [--max_insert_relative MAX_INSERT_RELATIVE]
                         [--min_insert MIN_INSERT] [-c MIN_COV] [-f MIN_FREQ]
                         [-fdr FDR] [-s MIN_SNP]
                         [--min_fasta_reads MIN_FASTA_READS]
                         [--store_everything] [--skip_mm_profiling]
                         [--scaffolds_to_profile SCAFFOLDS_TO_PROFILE]
                         bam fasta

  REQUIRED:
   bam                   Sorted .bam file
   fasta                 Fasta file the bam is mapped to

  I/O PARAMETERS:
   -o OUTPUT, --output OUTPUT
                         Output prefix (default: inStrain)

  SYSTEM PARAMETERS:
   -p PROCESSES, --processes PROCESSES
                         Number of processes to use (default: 6)
   -d, --debug           Make extra debugging output (default: False)
   -h, --help            show this help message and exit

  READ FILTERING OPTIONS:
   -l FILTER_CUTOFF, --filter_cutoff FILTER_CUTOFF
                         Minimum percent identity of read pairs to consensus to
                         use the reads. Must be >, not >= (default: 0.95)
   --min_mapq MIN_MAPQ   Minimum mapq score of EITHER read in a pair to use
                         that pair. Must be >, not >= (default: -1)
   --max_insert_relative MAX_INSERT_RELATIVE
                         Multiplier to determine maximum insert size between
                         two reads - default is to use 3x median insert size.
                         Must be >, not >= (default: 3)
   --min_insert MIN_INSERT
                         Minimum insert size between two reads - default is 50
                         bp. If two reads are 50bp each and overlap completely,
                         their insert will be 50. Must be >, not >= (default:
                         50)

  VARIANT CALLING OPTIONS:
   -c MIN_COV, --min_cov MIN_COV
                         Minimum coverage to call an variant (default: 5)
   -f MIN_FREQ, --min_freq MIN_FREQ
                         Minimum SNP frequency to confirm a SNV (both this AND
                         the FDR snp count cutoff must be true to call a SNP).
                         (default: 0.05)
   -fdr FDR, --fdr FDR   SNP false discovery rate- based on simulation data
                         with a 0.1 percent error rate (Q30) (default: 1e-06)

  OTHER OPTIONS:
   -s MIN_SNP, --min_snp MIN_SNP
                         Absolute minimum number of reads connecting two SNPs
                         to calculate LD between them. (default: 20)
   --min_fasta_reads MIN_FASTA_READS
                         Minimum number of reads mapping to a scaffold to
                         proceed with profiling it (default: 0)
   --store_everything    Store intermediate dictionaries in the pickle file;
                         will result in significantly more RAM and disk usage
                         (default: False)
   --skip_mm_profiling   Dont perform analysis on an mm level; saves RAM and
                         time (default: False)
   --scaffolds_to_profile SCAFFOLDS_TO_PROFILE
                         Path to a file containing a list of scaffolds to
                         profile- if provided will ONLY profile those scaffolds
                         (default: None)
