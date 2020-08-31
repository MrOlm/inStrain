User Manual
=============================

Generating inStrain input
----------------------------

There are two main inputs to inStrain: a :term:`fasta file` containing reference genome sequences, and a :term:`bam file` containing reads mapped to these sequences. Additionally and optionally, by providing a genes `.fna` file inStrain can calculate gene-level metrics, and by providing a scaffold-to-bin file inStrain can calculate metrics on a genome level. Here we go over some considerations involved in generating these inputs.

Preparing the .fasta file
+++++++++++++++++++++++++++++++

A :term:`fasta file` contains the DNA sequences of the contigs that you map your reads to. Choosing what :term:`fasta file` you will use (consensus / reference genomes) is important and will affect the interpretation of your *inStrain* results. Below we describe the three most common strategies.

Please note that the :term:`fasta file` provided to inStrain must always be the same as, or a subset of, the :term:`fasta file` used to create the :term:`bam file` (i.e. the :term:`fasta file` that reads were mapped to).

Using *de novo* assembled genomes (recommended)
`````````````````````````````````````````````````````

This strategy involves assembling genomes from the metagenomic samples that you'd like to profile. This is the recommended workflow for running *inStrain*:

1. Assemble reads into contigs for each sample collected from the environment. Recommended software: IDBA_UD, MEGAHIT, metaSPADES.

2. Bin genomes out of each assembly using differential coverage binning. Recommended software: Bowtie2 (for mapping), MetaBAT, CONCOCT, DasTOOL (for binning).

3. Dereplicate the entire set of genomes that you would like to profile (all genomes from all environments) at 97-99% identity, and filter out low quality genomes. Recommended software: dRep, checkM.

4. Create a bowtie2 index of the representative genomes from this dereplicated set and map reads to this set from each sample: Recommended software: Bowtie2

5. Profile the resulting mapping *.bam* files using inStrain.

6. Use *inStrain genome_wide* to calculate genome-level microdiveristy metrics for each originally binned genome.

An important aspect of this workflow is to **map to many genomes at once**. Mapping to just one genome at a time is highly discouraged, because this encourages :term:`mismapped reads<mismapped read>` from other genomes to be recruited by this genome. By including many (dereplicated) genomes in your bowtie2 index, you will be able to far more accurately filter out :term:`mismapped reads<mismapped read>` and reduce false positive SNPs.

Using a single genome .fasta file
``````````````````````````````````````
If your .fasta file is a single genome, the main consideration is that it should be a good representative genome for some organism in your sample. Ideally, it was assembled directly from that sample, isolated from that sample, or you have some other evidence that this genome is highly representation of a species in that sample. Regardless, you should check your `inStrain plot` output and `scaffold_info.tsv` output file to be sure that your inStrain run had decent coverage and breadth of coverage of the genome that you use before attempting to interpret the results.

Remember, your .fasta file can be a subset of the .fasta file that was used to create the .bam file. You can create a .bam with all dereplicated genomes from your environment, but then just pass a .fasta file for only the genomes of particular interest. This approach is recommended as opposed to creating a :term:`bam file` for just each genome, as it reduces :term:`mismapped reads<mismapped read>`

Using a metagenomic assembly
`````````````````````````````````
You can also pass inStrain an entire metagenomic assembly from a sample, including both binned and unbinned contigs. In this case, the output inStrain profile will include population information for each contig in the set. To  break it down by microbial genome / species, you can include a scaffold to bin file to generate results by genome.

Preparing the .bam file
++++++++++++++++++++++++++

InStrain is designed primarily for paired-end Illumina read sequencing, though un-paired reads can also be used by adjusting the run-time parameters. We recommend using the program Bowtie2 to map your reads to your genome.

Bowtie2 default parameters are what we use for mapping, but it may be worth playing around with them to see how different settings perform on your data. It is important to note that the ``-X`` flag (capital X) is the expected insert length and is by default ``500``. In many cases (e.g., 2x250 bp or simply datasets with longer inserts) it may be worthwhile to increase this value up to ``-X 1000`` for passing to Bowtie2. By default, if a read maps equally well to multiple genomes, Bowtie2 will pick one of the positions randomly and give the read a MAPQ score of 1. Thus, if you'd like to remove :term:`multi-mapped reads<multi-mapped read>`, you can set you minimum mapQ score to 2.

Other mapping software can also be used to generate .bam files for inStrain. However, some software (e.g. BBmap and SNAP) use the fasta file scaffold descriptions when generating the .bam files, which causes problems for inStrain. If using mapping software that does this, include the flag ``--use_full_fasta_header`` to let inStrain account for this.

.. note::
  If the reads that you'd like to run with inStrain are not working, please post an issue on GitHub. We're happy to upgrade inStrain to work with new mapping software and/or reads from different technologies.

Preparing the genes file
++++++++++++++++++++++++++

You can run prodigal on your :term:`fasta file` to generate an .fna file with the gene-level information. This .fna file can then be provided to inStrain profile to get gene-level characterizations.

Example::

 $ prodigal -i assembly.fasta -d genes.fna

Preparing a scaffold-to-bin file
++++++++++++++++++++++++++++++++++++++++++++++++++++

After running ``inStrain profile``, most results are presented on a scaffold-by-scaffold basis. There are a number of ways of telling *inStrain* which scaffold belongs to which genome, so that results can be analyzed on a genome-by-gene level as well.

1. Individual .fasta files. As recommended above, if you want to run *inStrain* on multiple genomes in the same sample, you should first concatenate all of the individual genomes into a single *.fasta* file and map to that. To view the results of the individual genomes used to create the concatenated .fasta file, you can pass a list of the individual .fasta files the ``-s`` arguement.

2. Scaffold-to-bin file. This is a text file consists of two columns, with one column listing the scaffold name, and the second column listing the genome bin name. Columns should be separated by tabs. The script `parse_stb.py <https://github.com/MrOlm/drep/blob/master/helper_scripts/parse_stb.py>`_  can help you create a scaffold-to-bin file from a list of individual .fasta files, or to split a concatenated .fasta file into individual genomes. The script comes packaged with the program `dRep <https://github.com/MrOlm/drep>`_, and can be installed with the command ``pip install drep``.

3. Nothing. If all of your scaffolds belong to the same genome, by running ``inStrain profile`` without any *-s* options it will summarize the results of all scaffolds together as if they all belong to the same genome.


Description of inStrain modules and arguments
----------------------------------------------

The functionality of inStrain is broken up into modules. To see a list of available modules, check the help::

    $ inStrain -h

                    ...::: inStrain v1.3.2 :::...

      Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

      Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
      Example: inStrain profile -h

      Workflows:
        profile         -> Create an inStrain profile (microdiversity analysis) from a mapping.
        compare         -> Compare multiple inStrain profiles (popANI, coverage_overlap, etc.)

      Single operations:
        profile_genes   -> Calculate gene-level metrics on an inStrain profile [DEPRECATED; USE profile INSTEAD]
        genome_wide     -> Calculate genome-level metrics on an inStrain profile
        quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
        filter_reads    -> Commands related to filtering reads from .bam files
        plot            -> Make figures from the results of "profile" or "compare"
        other           -> Other miscellaneous operations


profile
+++++++++++++

Module description
````````````````````

The most complex part of inStrain, and must be run before any other modules can be. The input is a :term:`fasta file` and a :term:`bam file`, and the output is an :term:`IS_profile<inStrain profile>`. The functionality of ``inStrain profile`` is broken into several steps.

First, all reads in the .bam file are filtered to only keep those that map with sufficient quality. All non-paired reads will be filtered out by default, and an additional set of filters are applied to each read pair (not the individual reads). Command line parameters can be adjusted to change the specifics, but in general:

* Pairs must be mapped in the proper orientation with an expected insert size. The minimum insert distance can be set with a command line parameter. The maximum insert distance is a multiple of the median insert distance. So if pairs have a median insert size of 500bp, by default all pairs with insert sizes over 1500bp will be excluded. For the max insert cutoff, the median_insert for all scaffolds is used.

* Pairs must have a minimum mapQ score. MapQ scores are confusing and how they're calculated varies based on the mapping algorithm being used, but are meant to represent both the number of mismatches in the mapping and how unique that mapping is. With bowtie2, if the read maps equally well to two positions on the genome (:term:`multi-mapped read`), its mapQ score will be set to 2. The read in the pair with the higher mapQ is used for the pair.

* Pairs must be above some minimum nucleotide identity (ANI) value. For example if reads in a pair are 100bp each, and each read has a single mismatch, the ANI of that pair would be 0.99

Next, using only read pairs that pass filters, a number of microdiveristy metrics are calculated on a scaffold-by-scaffold basis. This includes:

* Calculate the coverage at each position along the scaffold

* Calculate the :term:`nucleotide diversity` at each position along the scaffold in which the coverage is greater than the min_cov argument.

* Identify :term:`SNSs<SNS>` and :term:`SNVs<SNV>`. The criteria for being reported as a :term:`divergent site` are 1) More than min_cov number of bases at that position, 2) More than min_freq percentage of reads that are a variant base, 3) The number of reads with the variant base is more than the :term:`null model` for that coverage.

* Calculate :term:`linkage` between :term:`divergent sites<divergent site>` on the same read pair. For each pair harboring a :term:`divergent site`, calculate the linkage of that site with other :term:`divergent sites<divergent site>` within that same pair. This is only done for pairs of :term:`divergent sites<divergent site>` that are both on at least MIN_SNP reads

* Calculate scaffold-level properties. These include things like the overall coverage, breadth of coverage, average nucleotide identity (ANI) between the reads and the reference genome, and the expected breadth of coverage based on that true coverage.

Finally, this information is stored as an :term:`IS_profile<inStrain profile>` object. This includes the locations of :term:`divergent sites<divergent site>`, the number of read pairs that passed filters (and other information) for each scaffold, the linkage between SNV pairs, ect.

Module parameters
````````````````````
To see the command-line arguments for inStrain profile, check the help::

    $ inStrain profile -h
    usage: inStrain profile [-o OUTPUT] [--use_full_fasta_header] [-p PROCESSES]
                            [-d] [-h] [--version] [-l MIN_READ_ANI]
                            [--min_mapq MIN_MAPQ]
                            [--max_insert_relative MAX_INSERT_RELATIVE]
                            [--min_insert MIN_INSERT]
                            [--pairing_filter {paired_only,all_reads,non_discordant}]
                            [--priority_reads PRIORITY_READS]
                            [--detailed_mapping_info] [-c MIN_COV] [-f MIN_FREQ]
                            [-fdr FDR] [-g GENE_FILE] [-s [STB [STB ...]]]
                            [--mm_level] [--skip_mm_profiling] [--database_mode]
                            [--min_scaffold_reads MIN_SCAFFOLD_READS]
                            [--min_genome_coverage MIN_GENOME_COVERAGE]
                            [--min_snp MIN_SNP] [--store_everything]
                            [--scaffolds_to_profile SCAFFOLDS_TO_PROFILE]
                            [--rarefied_coverage RAREFIED_COVERAGE]
                            [--window_length WINDOW_LENGTH] [--skip_genome_wide]
                            [--skip_plot_generation]
                            bam fasta

    REQUIRED:
      bam                   Sorted .bam file
      fasta                 Fasta file the bam is mapped to

    I/O PARAMETERS:
      -o OUTPUT, --output OUTPUT
                            Output prefix (default: inStrain)
      --use_full_fasta_header
                            Instead of using the fasta ID (space in header before
                            space), use the full header. Needed for some mapping
                            tools (including bbMap) (default: False)

    SYSTEM PARAMETERS:
      -p PROCESSES, --processes PROCESSES
                            Number of processes to use (default: 6)
      -d, --debug           Make extra debugging output (default: False)
      -h, --help            show this help message and exit
      --version             show program's version number and exit

    READ FILTERING OPTIONS:
      -l MIN_READ_ANI, --min_read_ani MIN_READ_ANI
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
      --pairing_filter {paired_only,all_reads,non_discordant}
                            How should paired reads be handled?
                            paired_only = Only paired reads are retained
                            non_discordant = Keep all paired reads and singleton reads that map to a single scaffold
                            all_reads = Keep all reads regardless of pairing status (NOT RECOMMENDED; See documentation for deatils)
                             (default: paired_only)
      --priority_reads PRIORITY_READS
                            The location of a list of reads that should be
                            retained regardless of pairing status (for example
                            long reads or merged reads). This can be a .fastq file
                            or text file with list of read names (will assume file
                            is compressed if ends in .gz (default: None)

    READ OUTPUT OPTIONS:
      --detailed_mapping_info
                            Make a detailed read report indicating deatils about
                            each individual mapped read (default: False)

    VARIANT CALLING OPTIONS:
      -c MIN_COV, --min_cov MIN_COV
                            Minimum coverage to call an variant (default: 5)
      -f MIN_FREQ, --min_freq MIN_FREQ
                            Minimum SNP frequency to confirm a SNV (both this AND
                            the FDR snp count cutoff must be true to call a SNP).
                            (default: 0.05)
      -fdr FDR, --fdr FDR   SNP false discovery rate- based on simulation data
                            with a 0.1 percent error rate (Q30) (default: 1e-06)

    GENE PROFILING OPTIONS:
      -g GENE_FILE, --gene_file GENE_FILE
                            Path to prodigal .fna genes file. If file ends in .gb
                            or .gbk, will treat as a genbank file (EXPERIMENTAL;
                            the name of the gene must be in the gene qualifier)
                            (default: None)

    GENOME WIDE OPTIONS:
      -s [STB [STB ...]], --stb [STB [STB ...]]
                            Scaffold to bin. This can be a file with each line
                            listing a scaffold and a bin name, tab-seperated. This
                            can also be a space-seperated list of .fasta files,
                            with one genome per .fasta file. If nothing is
                            provided, all scaffolds will be treated as belonging
                            to the same genome (default: [])

    READ ANI OPTIONS:
      --mm_level            Create output files on the mm level (see documentation
                            for info) (default: False)
      --skip_mm_profiling   Dont perform analysis on an mm level; saves RAM and
                            time; impacts plots and raw_data (default: False)

    PROFILE OPTIONS:
      --database_mode       Set a number of parameters to values appropriate for
                            mapping to a large fasta file. Will set:
                            --min_read_ani 0.92 --skip_mm_profiling
                            --min_genome_coverage 1 (default: False)
      --min_scaffold_reads MIN_SCAFFOLD_READS
                            Minimum number of reads mapping to a scaffold to
                            proceed with profiling it (default: 1)
      --min_genome_coverage MIN_GENOME_COVERAGE
                            Minimum number of reads mapping to a genome to proceed
                            with profiling it. MUST profile .stb if this is set
                            (default: 0)
      --min_snp MIN_SNP     Absolute minimum number of reads connecting two SNPs
                            to calculate LD between them. (default: 20)
      --store_everything    Store intermediate dictionaries in the pickle file;
                            will result in significantly more RAM and disk usage
                            (default: False)
      --scaffolds_to_profile SCAFFOLDS_TO_PROFILE
                            Path to a file containing a list of scaffolds to
                            profile- if provided will ONLY profile those scaffolds
                            (default: None)
      --rarefied_coverage RAREFIED_COVERAGE
                            When calculating nucleotide diversity, also calculate
                            a rarefied version with this much coverage (default:
                            50)
      --window_length WINDOW_LENGTH
                            Break scaffolds into windows of this length when
                            profiling (default: 10000)

    OTHER  OPTIONS:
      --skip_genome_wide    Do not generate tables that consider groups of
                            scaffolds belonging to genomes (default: False)
      --skip_plot_generation
                            Do not make plots (default: False)


compare
+++++++++++++

Module description
````````````````````

Compare provides the ability to compare multiple :term:`inStrain profiles<inStrain profile>` (created by running ``inStrain profile``).

.. note::
  *inStrain* can only compare :term:`inStrain profiles<inStrain profile>`that have been mapped to the same .fasta file

``inStrain compare`` does pair-wise comparisons between each input :term:`inStrain profile<inStrain profile>`. For each pair, a series of steps are undertaken.

1. All positions in which both IS_profile objects have at least *min_cov* coverage (5x by default) are identified. This information can be stored in the output by using the flag *--store_coverage_overlap*, but due to it's size, it's not stored by default

2. Each position identified in step 1 is compared to calculate both :term:`conANI` and :term:`popANI`. The way that it compares positions is by testing whether the consensus base in sample 1 is detected at all in sample 2 and vice-verse. Detection of an allele in a sample is based on that allele being above the set *-min_freq* and *-fdr*. All detected differences between each pair of samples can be reported if the flag *--store_mismatch_locations* is set.

3. The coverage overlap and the average nucleotide identify for each scaffold is reported. For details on how this is done, see :doc:`example_output`

Module parameters
````````````````````
To see the command-line options, check the help::


    $ inStrain compare -h
    usage: inStrain compare -i [INPUT [INPUT ...]] [-o OUTPUT] [-p PROCESSES] [-d]
                            [-h] [--version] [-c MIN_COV] [-f MIN_FREQ] [-fdr FDR]
                            [-s SCAFFOLDS] [--store_coverage_overlap]
                            [--store_mismatch_locations]
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
      --version             show program's version number and exit

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

Other modules
++++++++++++++

The other modules are not commonly used, and mainly provide auxiliary functions or allow you run certain steps of ``profile`` after the fact. It is recommended to provide a genes file and/or a scaffold-to-bin file during ``inStrain profile`` rather than using ``profile_genes`` or ``genome_wide``, as it is more computationally efficient to do things that way.

profile_genes
````````````````````
After running *inStrain profile* on a sample, you can providing a file of gene calls to calculate gene-level metrics after the fact. This is less computationally efficient than providing the genes file to ``profile`` in the first place, however. See above for information about creating the input file.

To see the command-line options, check the help::

    $ inStrain profile_genes -h
    usage: inStrain profile_genes [-g GENE_FILE] -i IS [--store_everything]
                                  [-p PROCESSES] [-d] [-h] [--version]

    GENE PROFILING OPTIONS:
      -g GENE_FILE, --gene_file GENE_FILE
                            Path to prodigal .fna genes file. If file ends in .gb
                            or .gbk, will treat as a genbank file (EXPERIMENTAL;
                            the name of the gene must be in the gene qualifier)
                            (default: None)

    INPUT / OUTPUT:
      -i IS, --IS IS        an inStrain profile object (default: None)
      --store_everything    Store gene sequences in the IS object (default: False)

    SYSTEM PARAMETERS:
      -p PROCESSES, --processes PROCESSES
                            Number of processes to use (default: 6)
      -d, --debug           Make extra debugging output (default: False)
      -h, --help            show this help message and exit
      --version             show program's version number and exit

genome_wide
````````````````````
After running *inStrain profile* on a sample, you can provide information to analyze things on the genome level after the fact. This is less computationally efficient than providing the stb file to ``profile`` in the first place, however. See above for information about creating the input file.

To see the command-line options, check the help::

    $ inStrain genome_wide -h
    usage: inStrain genome_wide [-s [STB [STB ...]]] -i IS [--store_everything]
                                [--mm_level] [--skip_mm_profiling] [-p PROCESSES]
                                [-d] [-h] [--version]

    GENOME WIDE OPTIONS:
      -s [STB [STB ...]], --stb [STB [STB ...]]
                            Scaffold to bin. This can be a file with each line
                            listing a scaffold and a bin name, tab-seperated. This
                            can also be a space-seperated list of .fasta files,
                            with one genome per .fasta file. If nothing is
                            provided, all scaffolds will be treated as belonging
                            to the same genome (default: [])

    INPUT / OUTPUT:
      -i IS, --IS IS        an inStrain profile object (default: None)
      --store_everything    Store gene sequences in the IS object (default: False)

    READ ANI OPTIONS:
      --mm_level            Create output files on the mm level (see documentation
                            for info) (default: False)
      --skip_mm_profiling   Dont perform analysis on an mm level; saves RAM and
                            time; impacts plots and raw_data (default: False)

    SYSTEM PARAMETERS:
      -p PROCESSES, --processes PROCESSES
                            Number of processes to use (default: 6)
      -d, --debug           Make extra debugging output (default: False)
      -h, --help            show this help message and exit
      --version             show program's version number and exit

quick_profile
````````````````````

This is a quirky module that is not really related to any of the others. It is used to quickly profile a :term:`bam file` to pull out scaffolds from genomes that are at a sufficient breadth.

To use it you must provide a *.bam* file, the *.fasta* file that you mapped to to generate the *.bam* file, and a *scaffold to bin* file (see above section for details). The *stringent_breadth_cutoff* removed scaffolds entirely which have less breath than this (used to make the program run faster and produce smaller output). All scaffolds from genomes with at least the *breadth_cutoff* are then written to a file. In this way, you can then choose to run inStrain profile only on scaffolds from genomes that known to be of sufficient breadth, speeding up the run and reducing RAM usage (though not by much).

On the backend, this module is really just calling the program coverM `coverM<https://github.com/wwood/CoverM>`_

To see the command-line options, check the help::

    $ inStrain quick_profile -h
    usage: inStrain quick_profile [-p PROCESSES] [-d] [-h] [--version]
                                  [-s [STB [STB ...]]] [-o OUTPUT]
                                  [--breadth_cutoff BREADTH_CUTOFF]
                                  [--stringent_breadth_cutoff STRINGENT_BREADTH_CUTOFF]
                                  bam fasta

    REQUIRED:
      bam                   Sorted .bam file
      fasta                 Fasta file the bam is mapped to

    SYSTEM PARAMETERS:
      -p PROCESSES, --processes PROCESSES
                            Number of processes to use (default: 6)
      -d, --debug           Make extra debugging output (default: False)
      -h, --help            show this help message and exit
      --version             show program's version number and exit

    OTHER OPTIONS:
      -s [STB [STB ...]], --stb [STB [STB ...]]
                            Scaffold to bin. This can be a file with each line
                            listing a scaffold and a bin name, tab-seperated. This
                            can also be a space-seperated list of .fasta files,
                            with one genome per .fasta file. If nothing is
                            provided, all scaffolds will be treated as belonging
                            to the same genome (default: [])
      -o OUTPUT, --output OUTPUT
                            Output prefix (default: QuickProfile)
      --breadth_cutoff BREADTH_CUTOFF
                            Minimum genome breadth to pull scaffolds (default:
                            0.5)
      --stringent_breadth_cutoff STRINGENT_BREADTH_CUTOFF
                            Minimum breadth to let scaffold into coverm raw
                            results (done with greater than; NOT greater than or
                            equal to) (default: 0.0)

plot
````````````````````

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
````````````````````

This module holds odds and ends functionalities. As of version 1.3.1, all it can do is convert old *IS_profile* objects (>v0.3.0) to newer versions (v0.8.0). As the code base around *inStrain* matures, we expect more functionalities to be included here.

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
