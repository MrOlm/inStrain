Program documentation
===================

Preparing input
----------------

There are two simple inputs to *inStrain*: a *.fasta* file and a mapping file in *.bam* format. A third, a prodigal `.faa` file, can be used in later steps  Here we go over some considerations involved in choosing these inputs.

Preparing the .fasta file
+++++++++++++

A *.fasta* file contains the DNA sequences of the contigs that you map your reads to. Choosing what *.fasta* you will use (consensus / reference genomes) is extremely important and will affect the interpretation of your *inStrain* results. Below we describe the three most common strategies.

Please note that the *.fasta* file must always be the same as, or a subset of, the *.fasta* file used to create the *.bam* file, i.e. the *.fasta* file that reads were mapped to.

Using a collection of genomes (recommended)
```````````

The recommended workflow for running *inStrain*:

1. Assemble reads into contigs for each sample collected from the environment. Recommended software: IDBA_UD, MEGAHIT, metaSPADES.

2. Bin genomes out of each assembly using differential coverage binning. Recommended software: Bowtie2 (for mapping), MetaBAT, CONCOCT, DasTOOL (for binning).

3. Dereplicate the entire set of genomes across all samples from the environment at 97-99% identity, filter genome set draft quality genomes. Recommended software: dRep, checkM.

4. Create a bowtie2 index of the representative genomes from this dereplicated set and map reads to this set from each sample: Recommended software: Bowtie2

5. Profile the resulting mapping *.bam* files using inStrain.

6. Use *inStrain genome_wide* to calculate genome-level microdiveristy metrics for each originally binned genome.

The most important aspect of this workflow is to **map to many genomes at once**. Mapping to just one genome at a time is highly discouraged, because this encourages mismapped reads from other genomes to be recruited by this genome. By including many (dereplicated) genomes in your bowtie2 index, you will be able to far more accurately filter out mismapped reads and reduce false positive SNPs.

For more information on this, see :doc:`choosing_parameters`

Using a single genome .FASTA file
```````````

If your .fasta file is a single genome, the main consideration is that it should be a good representitive genome for some organism in your sample. Ideally, it was assembled directly from that sample, isolated from that sample, or you have some other evidence that this genome is highly representation of a species in that sample. Regardless, you should check your `inStrain plot` output and `scaffold_info.tsv` output file to be sure that your inStrain run had decent coverage and breadth of coverage of the genome that you use before attempting to interpret the results.

Remember, your *.fasta* file can be a subset of the *.fasta* file that was used to create the *.bam* file. You can create a BAM with all dereplicated genomes from your environment, but then just pass a *.fasta* file for only the genomes of particular interest. This approach is recommended as opposed to creating a BAM for just each genome, as it reduces mismapping.

Using a metagenomic assembly
```````````

You can also pass *inStrain* an entire metagenomic assembly from a sample, including either binned or unbinned contigs. In this case, the output inStrain profile will include population information for each contig in the set. To then break it down by microbial genome / species, You can use ``inStrain genome_wide`` including a scaffold to bin file to generate results by genome.

Preparing the .bam file
+++++++++++++

*inStrain* requires paired-end Illumina read sequencing. We recommend using Bowtie2 to map your reads to your genome.

Bowtie2 default parameters are what we use for mapping, but it may be worth playing around with them to see how different settings perform on your data. It is important to note that the `-X` flag (capital X) is the expected insert length and is by default `500`. In many cases (e.g., 2x250 bp or simply datasets with longer inserts) it may be worthwhile to increase this value up to `-X 1000` for passing to bowtie2.

Preparing the prodigal `.fna` genes file for gene-level profiling
+++++++++++++

You can run prodigal on your *.fasta* file to generate the *.fna* file with the gene-level information that `inStrain profile_genes` requires.

Example:

 $ prodigal -i assembly.fasta -d genes.fna

Module descriptions
----------------

The functionality of inStrain is broken up into modules. To see a list of available modules, check the help::

  $ inStrain -h

               ...::: inStrain v1.0.0 :::...

  Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

  Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
  Example: inStrain profile -h

   profile         -> Create an inStrain profile (microdiversity analysis) from a mapping.
   compare         -> Compare multiple inStrain profiles (popANI, coverage_overlap, etc.)
   profile_genes   -> Calculate gene-level metrics on an inStrain profile
   genome_wide     -> Calculate genome-level metrics on an inStrain profile
   quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
   filter_reads    -> Commands related to filtering reads from .bam files
   plot            -> Make figures from the results of "profile" or "compare"
   other           -> Other miscellaneous operations

IS_profile
+++++++++++++

An IS_profile (inStrain profile) is created by running the `inStrain profile` command. It contains  all of the program's internal workings, cached data, and output is stored. Additional modules can then be run on an IS_profile (to analyze genes, compare profiles, etc.), and there is lots of nice cached data stored in it that can be accessed using python.

.. seealso::

:doc:`example_output`
  For help finding where the output from your run is located in the IS_profile

:doc:`Advanced_use`
  For access to the raw internal data (which can be very useful)


profile
+++++++++++++

The most complex part of inStrain, and must be run before any other modules can be. The functionality of *profile* is broken into several steps.

First, all reads in the .bam file are filtered to only keep those that map with sufficient quality. Reads must be paired (all non-paired reads will be filtered out) and an additional set of filters are applied to the read pair (not the individual reads). Command line parameters can be adjusted to change the specifics, but in general:

 * Pairs must be mapped in the proper orientation with an expected insert size. The minimum insert distance can be set with a command line parameter. The maximum insert distance is a multiple of the median insert distance. So if pairs have a median insert size of 500bp, by default all pairs with insert sizes over 1500bp will be excluded.

 * Pairs must have a minimum mapQ score. MapQ scores are confusing and how they're calculated varies based on the mapping algorithm being used, but are meant to represent both the number of mismatches in the mapping and how unique that mapping is. With bowtie2, if the read maps equally well to two positions on the genome, its mapQ score will be set to 2. The read in the pair with the higher mapQ is used for the pair.

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
    For help interpreting the output

  :doc:`Advanced_use`
    For access to the raw internal data (which can be very useful)

  :doc:`choosing_parameters`
    For information about the pitfalls and other things to consider when running inStrain

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

compare
+++++++++++++

Compare provides the ability to compare two *IS_profile* folders (created by running *inStrain profile*). Both *IS_profile* objects must created based on mapping to the same *.bam* file for *compare* to work.

*inStrain compare* compares a set of different *IS_profile* folders (created by running *inStrain profile*). These *IS_profile* folders represent sets of different sample reads mapped to the same *.fasta* file. To use, we recommend assembly and binning of each sample, and then dereplication of genomes using the software dRep (https://drep.readthedocs.io/) at a high percent ANI, e.g. 96%-99%. Samples which contain multiple populations of the same dRep cluster (members of similar species or sub-species) can then be mapped back to the best genome from this dRep cluster, and then inStrain should be run on these dRep cluster genomes.

.. note::
  *inStrain* can only compare read profiles that have been mapped to the same .fasta file

Compare does pair-wise comparisons between each input *IS_profile*. For each pair, a series of steps are undertaken.

1. All positions in which both *IS_profile* objects have at least *min_cov* coverage (5x by default) are identified. This information can be stored in the output by using the flag *--store_coverage_overlap*, but due to it's size, it's not stored by default

2. Each position identified in step 1 is compared. If the flag *--compare_consensus_bases* is used, the consensus base at each position is compared. That means that if the position is 60% A 40% G in sample 1, and 40% A 60% G in sample 2, they will considered different. By default, however, this position would be considered the same. The way that is compares positions is by testing whether the consensus base in sample 1 is detected at all in sample 2 and vice-verse. Detection of an allele in a sample is based on that allele being above the set *-min_freq* and *-fdr*. All detected differences between each pair of samples can be reported if the flag *--store_mismatch_locations* is set.

3. The coverage overlap and the average nucleotide identify for each scaffold is reported. For details on how this is done, see :doc:`example_output`


To see the command-line options, check the help::

  $ inStrain compare -h
  usage: inStrain compare -i [INPUT [INPUT ...]] [-o OUTPUT] [-p PROCESSES] [-d]
                          [-h] [-c MIN_COV] [-f MIN_FREQ] [-fdr FDR]
                          [-s SCAFFOLDS] [--store_coverage_overlap]
                          [--store_mismatch_locations]
                          [--compare_consensus_bases]
                          [--include_self_comparisons] [--greedy_clustering]
                          [--g_ani G_ANI] [--g_cov G_COV] [--g_mm G_MM]

  REQUIRED:
    -i [INPUT [INPUT ...]], --input [INPUT [INPUT ...]]
                          A list of inStrain objects, all mapped to the same
                          .fasta file (default: None)
    -o OUTPUT, --output OUTPUT
                          Output prefix (default: instrainComparer)

  SYSTEM PARAMETERS:
    -p PROCESSES, --processes PROCESSES
                          Number of processes to use (default: 6)
    -d, --debug           Make extra debugging output (default: False)
    -h, --help            show this help message and exit

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
    -s SCAFFOLDS, --scaffolds SCAFFOLDS
                          Location to a list of scaffolds to compare. You can
                          also make this a .fasta file and it will load the
                          scaffold names (default: None)
    --store_coverage_overlap
                          Also store coverage overlap on an mm level (default:
                          False)
    --store_mismatch_locations
                          Store the locations of SNPs (default: False)
    --compare_consensus_bases
                          Only compare consensus bases; dont look for lower
                          frequency SNPs when calculating ANI (default: False)
    --include_self_comparisons
                          Also compare IS profiles against themself (default:
                          False)

  GREEDY CLUSTERING OPTIONS [THIS SECTION IS EXPERIMENTAL!]:
    --greedy_clustering   Dont do pair-wise comparisons, do greedy clustering to
                          only find the number of clsuters. If this is set, use
                          the parameters below as well (default: False)
    --g_ani G_ANI         ANI threshold for greedy clustering- put the fraction
                          not the percentage (e.g. 0.99, not 99) (default: 0.99)
    --g_cov G_COV         Alignment coverage for greedy clustering- put the
                          fraction not the percentage (e.g. 0.5, not 10)
                          (default: 0.99)
    --g_mm G_MM           Maximum read mismatch level (default: 100)

profile_genes
+++++++++++++

After running *inStrain profile* on a sample, you can calculate the coverage, microdiveristy, and SNP type for each gene. You do this by providing a file of gene calls. See doc:`example_output` for example results, and doc:`preparing_input` for information about creating the input file.

To see the command-line options, check the help::

  $ inStrain profile_genes -h
  usage: inStrain profile_genes -i IS -g GENE_FILE [-p PROCESSES] [-d] [-h]

  REQUIRED:
  -i IS, --IS IS        an inStrain profile object (default: None)
  -g GENE_FILE, --gene_file GENE_FILE
                        Path to prodigal .fna genes file. (default: None)

  SYSTEM PARAMETERS:
  -p PROCESSES, --processes PROCESSES
                        Number of processes to use (default: 6)
  -d, --debug           Make extra debugging output (default: False)
  -h, --help            show this help message and exit

genome_wide
+++++++++++++

After running *inStrain profile*, most results are presented on a scaffold-by-scaffold basis. To have the results summarized in a genome-by-genome way instead, you can use the module *inStrain genome_wide*. It is also required to run this module before making plots.

There are a number of ways of telling *inStrain* which scaffold belongs to which genome

1. Individual .fasta files. As recommended in :doc:`preparing_input`, if you want to run *inStrain* on multiple genomes in the same sample, you should first concatenate all of the individual genomes into a single *.fasta* file and map to that. To view the results of the individual genomes used to create the concatenated .fasta file, you can pass a list of the individual *.fasta files to *inStrain genome_wide*. (e.g. inStrain genome_wide -i inStrain_folder -s genome1.fasta genome2.fasta genome3.fasta)

2. Scaffold to bin file. This text file consists of two columns, with one column listing the scaffold name, and the second column listing the genome bin name. Columns should be separated by tabs.

3. Nothing. If all of your scaffolds belong to the same genome, by running *inStrain genome_wide* without any *-s* options it will summarize the results of all scaffolds together.

The flag `--mm_level` produces output for each mm. You probably don't want this. For information on what I mean by mm_level see :doc:`Advanced_use`, for information on the output see :doc:`example_output`

To see the command-line options, check the help::

  $ inStrain genome_wide -h
  usage: inStrain genome_wide -i IS [-s [STB [STB ...]]] [--mm_level]
                            [-p PROCESSES] [-d] [-h]

  REQUIRED:
  -i IS, --IS IS        an inStrain profile object (default: None)
  -s [STB [STB ...]], --stb [STB [STB ...]]
                        Scaffold to bin. This can be a file with each line
                        listing a scaffold and a bin name, tab-seperated. This
                        can also be a space-seperated list of .fasta files,
                        with one genome per .fasta file. If nothing is
                        provided, all scaffolds will be treated as belonging
                        to the same genome (default: [])
  --mm_level            Create files on the mm level (see documentation for
                        info) (default: False)

  SYSTEM PARAMETERS:
  -p PROCESSES, --processes PROCESSES
                        Number of processes to use (default: 6)
  -d, --debug           Make extra debugging output (default: False)
  -h, --help            show this help message and exit

quick_profile
+++++++++++++

This is a quirky module that is not really related to any of the others. It is used to quickly profile a *.bam* file to pull out scaffolds from genomes that are at a sufficient breadth.

To use it you must provide a *.bam* file, the *.fasta* file that you mapped to to generate the *.bam* file, and a *scaffold to bin* file (see above section for details). The *stringent_breadth_cutoff* removed scaffolds entirely which have less breath than this (used to make the program run faster and produce smaller output). All scaffolds from genomes with at least the *breadth_cutoff* are then written to a file. In this way, you can then choose to run inStrain profile only on scaffolds from genomes that known to be of sufficient breadth, speeding up the run and reducing RAM usage (though not by much).

To see the command-line options, check the help::

  $ inStrain quick_profile -h
  usage: inStrain quick_profile -b BAM -f FASTA -s STB [-o OUTPUT]
                              [-p PROCESSES] [-d] [-h]
                              [--breadth_cutoff BREADTH_CUTOFF]
                              [--stringent_breadth_cutoff STRINGENT_BREADTH_CUTOFF]

  REQUIRED:
  -b BAM, --bam BAM     A bam file to profile (default: None)
  -f FASTA, --fasta FASTA
                        The .fasta file to profile (default: None)
  -s STB, --stb STB     Scaffold to bin file for genome-wide coverage and
                        breadth (default: None)
  -o OUTPUT, --output OUTPUT
                        Output prefix (default: None)

  SYSTEM PARAMETERS:
  -p PROCESSES, --processes PROCESSES
                        Number of processes to use (default: 6)
  -d, --debug           Make extra debugging output (default: False)
  -h, --help            show this help message and exit

  OTHER OPTIONS:
  --breadth_cutoff BREADTH_CUTOFF
                        Minimum breadth to pull scaffolds (default: 0.5)
  --stringent_breadth_cutoff STRINGENT_BREADTH_CUTOFF
                        Minimum breadth to let scaffold into coverm raw
                        results (default: 0.01)

plot
+++++++++++++

This module produces plots based on the results of *inStrain profile* and *inStrain compare*. In both cases, before plots can be made, *inStrain genome_wide* must be run on the output folder first. In order to make plots 8 and 9, *inStrain profile_genes* must be run first as well.

The recommended way of running this module is with the default `-pl a`. It will just try and make all of the plots that it can, and will tell you about any plots that it fails to make.

See :doc:`example_output` for an example of the plots it can make.

To see the command-line options, check the help::

  $ inStrain plot -h
  usage: inStrain plot -i IS [-pl [PLOTS [PLOTS ...]]] [-p PROCESSES] [-d] [-h]

  REQUIRED:
    -i IS, --IS IS        an inStrain profile object (default: None)
    -pl [PLOTS [PLOTS ...]], --plots [PLOTS [PLOTS ...]]
                          Plots. Input 'all' or 'a' to plot all
                          1) Coverage and breadth vs. read mismatches
                          2) Genome-wide microdiversity metrics
                          3) Read-level ANI distribution
                          4) Major allele frequencies
                          5) Linkage decay
                          6) Read filtering plots
                          7) Scaffold inspection plot (large)
                          8) Linkage with SNP type (GENES REQUIRED)
                          9) Gene histograms (GENES REQUIRED)
                          10) Compare dendrograms (RUN ON COMPARE; NOT PROFILE)
                           (default: a)

  SYSTEM PARAMETERS:
    -p PROCESSES, --processes PROCESSES
                          Number of processes to use (default: 6)
    -d, --debug           Make extra debugging output (default: False)
    -h, --help            show this help message and exit

other
+++++++++++++

This module holds odds and ends functionalities. As of version 1.0.0, all it can do is convert old *IS_profile* objects (>v0.3.0) to newer versions (v0.8.0). As the code base around *inStrain* matures, we expect more functionalities to be included here.

To see the command-line options, check the help::

  $ inStrain other -h
  usage: inStrain other [-p PROCESSES] [-d] [-h] [--old_IS OLD_IS]

  SYSTEM PARAMETERS:
    -p PROCESSES, --processes PROCESSES
                          Number of processes to use (default: 6)
    -d, --debug           Make extra debugging output (default: False)
    -h, --help            show this help message and exit

  OTHER OPTIONS:
    --old_IS OLD_IS       Convert an old inStrain version object to the newer
                          version. (default: None)
