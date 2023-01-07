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

4. Create a :term:`scaffold-to-bin file` from the genome set. Recommended software: `parse_stb.py <https://github.com/MrOlm/drep/blob/master/helper_scripts/parse_stb.py>`_

5. Create a bowtie2 index of the representative genomes from this dereplicated set and map reads to this set from each sample. Recommended software: Bowtie2

6. Profile the resulting mapping *.bam* files using inStrain to calculate genome-level :term:`microdiveristy` metrics for each originally binned genome.

An important aspect of this workflow is to **map to many genomes at once**. Mapping to just one genome at a time is highly discouraged, because this encourages :term:`mismapped reads<mismapped read>` from other genomes to be recruited by this genome. By including many (dereplicated) genomes in your bowtie2 index, you will be able to far more accurately filter out :term:`mismapped reads<mismapped read>` and reduce false positive SNPs. See :doc:`important_concepts` for more info.

For instructions on merging your genomes with a public database, see Tutorial #3 of :doc:`tutorial`.

Using a single genome .fasta file
``````````````````````````````````````
If your .fasta file is a single genome, the main consideration is that it should be a good representative genome for some organism in your sample. Ideally, it was assembled directly from that sample, isolated from that sample, or you have some other evidence that this genome is highly representative of a species in that sample. Regardless, you should check your `inStrain plot` output and `scaffold_info.tsv` output file to be sure that your inStrain run had decent coverage and breadth of coverage of the genome that you use before attempting to interpret the results.

Remember, your .fasta file can be a subset of the .fasta file that was used to create the .bam file. You can create a .bam with all dereplicated genomes from your environment, but then just pass a .fasta file for only the genomes of particular interest. This approach is recommended as opposed to creating a :term:`bam file` for just each genome, as it reduces :term:`mismapped reads<mismapped read>`

Using a metagenomic assembly
`````````````````````````````````
You can also pass inStrain an entire metagenomic assembly from a sample, including both binned and unbinned contigs. In this case, the output inStrain profile will include population information for each contig in the set. To  break it down by microbial genome / species, you can include a :term:`scaffold-to-bin file` to generate results by genome.

Preparing the .bam file
++++++++++++++++++++++++++

InStrain is designed primarily for paired-end Illumina read sequencing, though un-paired reads can also be used by adjusting the run-time parameters. We recommend using the program Bowtie2 to map your reads to your genome.

Bowtie2 default parameters are what we use for mapping, but it may be worth playing around with them to see how different settings perform on your data. It is important to note that the ``-X`` flag (capital X) is the expected insert length and is by default ``500``. In many cases (e.g., 2x250 bp or simply datasets with longer inserts) it may be worthwhile to increase this value up to ``-X 1000`` for passing to Bowtie2. By default, if a read maps equally well to multiple genomes, Bowtie2 will pick one of the positions randomly and give the read a MAPQ score of 1. Thus, if you'd like to remove :term:`multi-mapped reads<multi-mapped read>`, you can set the minimum mapQ score to 2.

Other mapping software can also be used to generate .bam files for inStrain. However, some software (e.g. BBmap and SNAP) use the fasta file scaffold descriptions when generating the .bam files, which causes problems for inStrain. If using mapping software that does this, include the flag ``--use_full_fasta_header`` to let inStrain account for this.

.. note::
  If the reads that you'd like to run with inStrain are not working, please post an issue on GitHub. We're happy to upgrade inStrain to work with new mapping software and/or reads from different technologies.

Preparing the genes file
++++++++++++++++++++++++++

You can run prodigal on your :term:`fasta file` to generate an .fna file with the gene-level information. This .fna file can then be provided to inStrain profile to get gene-level characterizations.

Example::

 $ prodigal -i assembly.fasta -d genes.fna -a genes.faa

Preparing a scaffold-to-bin file
++++++++++++++++++++++++++++++++++++++++++++++++++++

After running ``inStrain profile``, most results are presented on a scaffold-by-scaffold basis. There are a number of ways of telling *inStrain* which scaffold belongs to which genome, so that results can be analyzed on a genome-by-gene level as well.

1. Individual .fasta files. As recommended above, if you want to run *inStrain* on multiple genomes in the same sample, you should first concatenate all of the individual genomes into a single *.fasta* file and map to that. To view the results of the individual genomes used to create the concatenated .fasta file, you can pass a list of the individual .fasta files the ``-s`` argument.

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

The heart of inStrain. The input is a :term:`fasta file` and a :term:`bam file`, and the output is an :term:`IS_profile<inStrain profile>`. The functionality of ``inStrain profile`` is broken into several steps.

First, all reads in the .bam file are filtered to only keep those that map with sufficient quality. All non-paired reads will be filtered out by default, and an additional set of filters are applied to each read pair (not the individual reads). Command line parameters can be adjusted to change the specifics, but in general:

* Pairs must be mapped in the proper orientation with an expected insert size. The minimum insert distance can be set with a command line parameter. The maximum insert distance is a multiple of the median insert distance. So if pairs have a median insert size of 500bp, by default all pairs with insert sizes over 1500bp will be excluded. For the max insert cutoff, the median_insert for all scaffolds is used.

* Pairs must have a minimum mapQ score. MapQ scores are confusing and how they're calculated varies based on the mapping algorithm being used, but are meant to represent both the number of mismatches in the mapping and how unique that mapping is. With bowtie2, if the read maps equally well to two positions on the genome (:term:`multi-mapped read`), its mapQ score will be set to 2. The read in the pair with the higher mapQ is used for the pair.

* Pairs must be above some minimum nucleotide identity (ANI) value. For example if reads in a pair are 100bp each, and each read has a single mismatch, the ANI of that pair would be 0.99

Next, using only read pairs that pass filters, a number of :term:`microdiversity` metrics are calculated on a scaffold-by-scaffold basis. This includes:

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
  *inStrain* can only compare inStrain profiles that have been mapped to the same .fasta file

``inStrain compare`` does pairwise comparisons between each input :term:`inStrain profile<inStrain profile>`. For each pair, a series of steps are undertaken.

1. All positions in which both IS_profile objects have at least *min_cov* coverage (5x by default) are identified. This information can be stored in the output by using the flag *--store_coverage_overlap*, but due to it's size, it's not stored by default

2. Each position identified in step 1 is compared to calculate both :term:`conANI` and :term:`popANI`. The way that it compares positions is by testing whether the consensus base in sample 1 is detected at all in sample 2 and vice-versa. Detection of an allele in a sample is based on that allele being above the set *-min_freq* and *-fdr*. All detected differences between each pair of samples can be reported if the flag *--store_mismatch_locations* is set.

3. The coverage overlap and the average nucleotide identity for each scaffold is reported. For details on how this is done, see :doc:`example_output`

4. **New in v1.6** Tables that list the coverage and base-frequencies of each SNV in all samples can be generated using the `--bams` parameter within inStrain compare. For each inStrain profile provided with the `-i` parameter, the corresponding bam file must be provided with the `--bams` parameter. The same read filtering parameters used in the original profile command will be used when running this analysis. See section `SNV POOLING OPTIONS:` in the help below for full information about this option, and see :doc:`example_output` for the tables it creates.

Module parameters
````````````````````
To see the command-line options, check the help::

    $ inStrain compare -h
    usage: inStrain compare -i [INPUT [INPUT ...]] [-o OUTPUT] [-p PROCESSES] [-d]
                        [-h] [--version] [-s [STB [STB ...]]] [-c MIN_COV]
                        [-f MIN_FREQ] [-fdr FDR] [--database_mode]
                        [--breadth BREADTH] [-sc SCAFFOLDS] [--genome GENOME]
                        [--store_coverage_overlap]
                        [--store_mismatch_locations]
                        [--include_self_comparisons] [--skip_plot_generation]
                        [--group_length GROUP_LENGTH] [--force_compress]
                        [-ani ANI_THRESHOLD] [-cov COVERAGE_TRESHOLD]
                        [--clusterAlg {centroid,weighted,ward,single,complete,average,median}]
                        [-bams [BAMS [BAMS ...]]] [--skip_popANI]

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

    GENOME WIDE OPTIONS:
      -s [STB [STB ...]], --stb [STB [STB ...]]
                            Scaffold to bin. This can be a file with each line
                            listing a scaffold and a bin name, tab-seperated. This
                            can also be a space-seperated list of .fasta files,
                            with one genome per .fasta file. If nothing is
                            provided, all scaffolds will be treated as belonging
                            to the same genome (default: [])

    VARIANT CALLING OPTIONS:
      -c MIN_COV, --min_cov MIN_COV
                            Minimum coverage to call an variant (default: 5)
      -f MIN_FREQ, --min_freq MIN_FREQ
                            Minimum SNP frequency to confirm a SNV (both this AND
                            the FDR snp count cutoff must be true to call a SNP).
                            (default: 0.05)
      -fdr FDR, --fdr FDR   SNP false discovery rate- based on simulation data
                            with a 0.1 percent error rate (Q30) (default: 1e-06)

    DATABASE MODE PARAMETERS:
      --database_mode       Using the parameters below, automatically determine
                            which genomes are present in each Profile and only
                            compare scaffolds from those genomes. All profiles
                            must have run Profile with the same .stb (default:
                            False)
      --breadth BREADTH     Minimum breadth_minCov required to count a genome
                            present (default: 0.5)

    OTHER OPTIONS:
      -sc SCAFFOLDS, --scaffolds SCAFFOLDS
                            Location to a list of scaffolds to compare. You can
                            also make this a .fasta file and it will load the
                            scaffold names (default: None)
      --genome GENOME       Run scaffolds belonging to this single genome only.
                            Must provide an .stb file (default: None)
      --store_coverage_overlap
                            Also store coverage overlap on an mm level (default:
                            False)
      --store_mismatch_locations
                            Store the locations of SNPs (default: False)
      --include_self_comparisons
                            Also compare IS profiles against themself (default:
                            False)
      --skip_plot_generation
                            Dont create plots at the end of the run. (default:
                            False)
      --group_length GROUP_LENGTH
                            How many bp to compare simultaneously (higher will use
                            more RAM and run more quickly) (default: 10000000)
      --force_compress      Force compression of all output files (default: False)

    GENOME CLUSTERING OPTIONS:
      -ani ANI_THRESHOLD, --ani_threshold ANI_THRESHOLD
                            popANI threshold to cluster genomes at. Must provide
                            .stb file to do so (default: 0.99999)
      -cov COVERAGE_TRESHOLD, --coverage_treshold COVERAGE_TRESHOLD
                            Minimum percent_genome_compared for a genome
                            comparison to count; if below the popANI will be set
                            to 0. (default: 0.1)
      --clusterAlg {centroid,weighted,ward,single,complete,average,median}
                            Algorithm used to cluster genomes (passed to
                            scipy.cluster.hierarchy.linkage) (default: average)

    SNV POOLING OPTIONS:
      -bams [BAMS [BAMS ...]], --bams [BAMS [BAMS ...]]
                            Location of .bam files used during inStrain profile
                            commands; needed to pull low-frequency SNVs. MUST BE
                            IN SAME ORDER AS THE INPUT FILES (default: None)
      --skip_popANI         Only run SNV Pooling; skip other compare operations
                            (default: False)

Other modules
++++++++++++++

The other modules are not commonly used, and mainly provide auxiliary functions or allow you to run certain steps of ``profile`` after the fact. It is recommended to provide a genes file and/or a scaffold-to-bin file during ``inStrain profile`` rather than using ``profile_genes`` or ``genome_wide``, as it is more computationally efficient to do things that way.

quick_profile
````````````````````
This is a quirky module that is not really related to any of the others. It is used to quickly profile a :term:`bam file` to pull out scaffolds from genomes that are at a sufficient breadth. To use it you must provide a *.bam* file, the *.fasta* file that you mapped to to generate the *.bam* file, and a *scaffold to bin* file (see above section for details). On the backend this module is really just calling the program `coverM <https://github.com/wwood/CoverM>`_

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

This module holds odds and ends functionalities. As of version 1.4, all it can do is convert old *IS_profile* objects (>v0.3.0) to newer versions (v0.8.0) and create runtime summaries of complete inStrain runs. As the code base around *inStrain* matures, we expect more functionalities to be included here.

To see the command-line options, check the help::

    $  inStrain other -h
    usage: inStrain other [-p PROCESSES] [-d] [-h] [--version] [--old_IS OLD_IS]
                          [--run_statistics RUN_STATISTICS]

    SYSTEM PARAMETERS:
      -p PROCESSES, --processes PROCESSES
                            Number of processes to use (default: 6)
      -d, --debug           Make extra debugging output (default: False)
      -h, --help            show this help message and exit
      --version             show program's version number and exit

    OTHER OPTIONS:
      --old_IS OLD_IS       Convert an old inStrain version object to the newer
                            version. (default: None)
      --run_statistics RUN_STATISTICS
                            Generate runtime reports for an inStrain run.
                            (default: None)

Other related operations
----------------------------

The goal of this section is to describe how to perform other operations that are commonly part of an inStrain-based workflow.

Gene Annotation
+++++++++++++++++

Below are some potential ways of annotating genes for follow-up inStrain analysis. The input to all operations is an amino acid fasta file (`.faa`), which should match the `.fna` file you passed to inStrain (see :doc:`user_manual.rst#preparing-the-genes-file` for an example command)

**If you have some other annotation you like to use, please add to to this list by submitting a pull request on GitHub!** (https://github.com/MrOlm/inStrain/blob/master/docs/user_manual.rst)

KEGG Orthologies (KOs)
`````````````````````````
KOs can be annotated using KofamScan / KofamKOALA (https://www.genome.jp/tools/kofamkoala/)::

    # Download the database and executables
    wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

    # Unzip and untar
    gzip -d ko_list.gz
    tar xf profiles.tar.gz
    tar xf kofam_scan-1.3.0.tar.gz

    # Run kofamscan
    exec_annotation -p profiles -k ko_list --cpu 10 --tmp-dir ./tmp -o genes.faa.kofamscan genes.faa

The following python code parses the resulting table

.. code-block:: python

    import pandas as pd
    from collections import defaultdict

    def parse_kofamscan(floc):
        """
        v1.0: 1/6/2023
        
        Parse kofamscan results. Only save results where KO score > threshold
        
        Returns:
            Adb: DataFrame with KOfam results
        """
        table = defaultdict(list)
        
        with open(floc, 'r') as o:
            
            # This blockis for RAM efficiency
            while True:
                line = o.readline()
                if not line:
                    break
        
                line = line.strip()
                if line[0] == '#':
                    continue
                    
                lw = line.split()
                if lw[0] == '*':
                    del lw[0]
                    
                if lw[2] == '-':
                    lw[2] = 0
                    
                try:
                    if float(lw[3]) >= float(lw[2]):
                        g = lw[0]
                        k = lw[1]
                        
                        table['gene'].append(g)
                        table['KO'].append(k)
                        table['thrshld'].append(float(lw[2]))
                        table['score'].append(float(lw[3]))
                        table['e_value'].append(float(lw[4]))
                        table['KO_definition'].append(' '.join(lw[4:]))
                except:
                    print(line)
                    assert False   
        o.close()
        
        Adb = pd.DataFrame(table)
        return Adb

    floc = "genes.faa.kofamscan genes.faa"
    Adb = parse_kofamscan(floc)

Where Adb is a pandas DataFrame that looks like:

.. csv-table:: Adb

    gene,KO,thrshld,score,e_value,KO_definition
    AP010889.1_1,K02313,130.13,593.2,1.200000e-178,chromosomal replication initiator protein
    AP010889.1_2,K02338,52.73,345.9,1.300000e-103,DNA polymerase III subunit beta [EC:2.7.7.7]
    AP010889.1_2,K22359,0.00,12.5,3.000000e-02,alkene monooxygenase gamma subunit [EC:1.14.13...
    AP010889.1_3,K03629,115.43,397.5,1.600000e-119,DNA replication and repair protein RecF
    AP010889.1_5,K02470,946.10,986.6,1.500000e-297,DNA gyrase subunit B [EC:5.6.2.2]

Carbohydrate-Active enZYmes (CAZymes)
``````````````````````````````````````
CAZymes can be profiled using the HMMs provided by dbCAN, which are based on CAZyDB (http://www.cazy.org/)::

  # Download the HMMs and executables
  wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt
  wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/hmmscan-parser.sh

  # Prepare HMMs
  hmmpress dbCAN-HMMdb-V11.txt

  # Run (based on readme here - https://bcb.unl.edu/dbCAN2/download/Databases/V11/readme.txt)
  hmmscan --domtblout genes.faa_vs_dbCAN_v11.dm dbCAN-HMMdb-V11.txt genes.faa > /dev/null ; sh /hmmscan-parser.sh genes.faa_vs_dbCAN_v11.dm > genes.faa_vs_dbCAN_v11.dm.ps ; cat genes.faa_vs_dbCAN_v11.dm.ps | awk '$5<1e-15&&$10>0.35' > genes.faa_vs_dbCAN_v11.dm.ps.stringent

The following python code parses the resulting table

.. code-block:: python

  import pandas as pd
  from collections import defaultdict

  def parse_dbcan(floc):
      """
      v1.0 - 1/6/2023
      
      Parse dbcan2 results
      
      Returns:
          Cdb: DataFrame with dbCAN2 results
      """
      
      h = ['Family_HMM', 'HMM_length', 'gene', 'Query_length', 'E-value', 'HMM_start', 'HMM_end', 'Query_start', 'Query_end', 'Coverage']
      Zdb = pd.read_csv(floc, sep='\t', names=h)
      
      # Parse names
      def get_type(f):
          for start in ['PL', 'AA', 'GH', 'CBM', 'GT', 'CE']:
              if f.startswith(start):
                  return start
          if f in ['dockerin', 'SLH', 'cohesin']:
              return 'cellulosome'
          print(f)
          assert False

      def get_family(f):
          for start in ['PL', 'AA', 'GH', 'CBM', 'GT', 'CE']:
              if f.startswith(start):
                  if f == 'CBM35inCE17':
                      return 35
                  try:
                      return int(f.replace(start, '').split('_')[0])
                  except:
                      print(f)
                      assert False
          if f in ['dockerin', 'SLH', 'cohesin']:
              return f
          print(f)
          assert False

      def get_subfamily(f):
          if f.startswith('GT2_'):
              if f == 'GT2_Glycos_transf_2':
                  return 0
              else:
                  return f.split('_')[-1]
          if '_' in f:
              try:
                  return int(f.split('_')[1])
              except:
                  print(f)
                  assert False
          else:
              return 0

      t2n = {'GH':'glycoside hydrolases',
            'PL':'polysaccharide lyases',
            'GT':'glycosyltransferases',
            'CBM':'non-catalytic carbohydrate-binding modules',
            'AA':'auxiliary activities',
            'CE':'carbohydrate esterases',
            'cellulosome':'cellulosome'}    

      ZIdb = Zdb[['Family_HMM']].drop_duplicates()
      ZIdb['raw_family'] = [x.split('.')[0] for x in ZIdb['Family_HMM']]
      ZIdb['class'] = [get_type(f) for f in ZIdb['raw_family']]
      ZIdb['class_name'] = ZIdb['class'].map(t2n)
      ZIdb['family'] = [get_family(f) for f in ZIdb['raw_family']]
      ZIdb['subfamily'] = [get_subfamily(f) for f in ZIdb['raw_family']]

      ZIdb['CAZyme'] = [f"{c}{f}_{s}" for c, f, s in zip(ZIdb['class'], ZIdb['family'], ZIdb['subfamily'])]
      
      ZSdb = pd.merge(Zdb, ZIdb[['Family_HMM', 'class', 'family',
        'subfamily', 'CAZyme']], on='Family_HMM', how='left')
      
      # Reorder
      ZSdb = ZSdb[[
          'gene', 
          'CAZyme',
          'class',
          'family',
          'subfamily',
          'Family_HMM',
          'HMM_length',
          'Query_length',
          'E-value',
          'HMM_start',
          'HMM_end',
          'Query_start',
          'Query_end',
          'Coverage',
          ]]
      
      return ZSdb

  floc = "/LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/2022_Misame/gene_annotation/dbCAN/DeltaI_NewBifido.faa_vs_dbCAN_v11.dm.ps"
  Cdb = parse_dbcan(floc)

Where Cdb is a pandas DataFrame that looks like:

.. csv-table:: Cdb

  gene,CAZyme,class,family,subfamily,Family_HMM,HMM_length,Query_length,E-value,HMM_start,HMM_end,Query_start,Query_end,Coverage
  AP010888.1_103,GH51_0,GH,51,0,GH51.hmm,630,516,1.700000e-137,84,542,9,515,0.726984
  AP010888.1_107,GH13_18,GH,13,18,GH13_18.hmm,343,509,4.100000e-114,2,343,35,379,0.994169
  AP010888.1_113,GH13_30,GH,13,30,GH13_30.hmm,365,605,3.100000e-163,1,365,33,403,0.997260
  AP010888.1_115,GH77_0,GH,77,0,GH77.hmm,494,746,4.700000e-134,2,482,204,728,0.971660
  AP010888.1_120,GH31_0,GH,31,0,GH31.hmm,427,846,1.300000e-129,1,427,198,627,0.997658

Antibiotic Resistance Genes
``````````````````````````````````````
There are many, many different ways of identifying antibiotic resistance genes. The method below is based on identifying homologs to know antibiotic resistance genes using the CARD database (https://card.mcmaster.ca/download)::

  # Download and unzip database
  wget https://card.mcmaster.ca/download/0/broadstreet-v3.2.5.tar.bz2
  tar -xvjf broadstreet-v3.2.5.tar.bz2

  # Make a diamond database out of it
  diamond makedb --in protein_fasta_protein_homolog_model.fasta -d protein_fasta_protein_homolog_model.dmd --threads 6

  # Run
  diamond blastp -q genes.faa -d protein_fasta_protein_homolog_model.dmd -f 6 -e 0.0001 -k 1 -p 6 -o genes.faa_vs_CARD.dm

The following python code parses the resulting table

.. code-block:: python

  import pandas as pd
  from collections import defaultdict

  def parse_card(floc, jloc=None):
      """
      v1.0 - 1/6/2023
      
      Parse CARD 
      
      Returns:
          Rdb: DataFrame with CARD results
      """
      h = ['gene', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',
          'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',
          'extra']
      db = pd.read_csv(floc, sep='\t', names=h)
      del db['extra']
      
      db['protein_seq_accession'] =  [t.split('|')[1] for t in db['target']]
      db['ARO'] =  [t.split('|')[2].split(':')[-1] for t in db['target']]
      db['CARD_short_name'] =  [t.split('|')[3].split(':')[-1] for t in db['target']]
      
      # Reorder
      header = ['gene', 'CARD_short_name', 'ARO', 'target']
      db = db[header + [x for x in list(db.columns) if x not in header]]
      
      if jloc is None:
          return db
      
      # Parse more
      import json
      j = json.load(open(jloc))
      
      aro2name = {}
      aro2categories = {}

      for n, m2t in j.items():
          if type(m2t) != type({}):
              continue

          if 'ARO_description' in m2t: 
              aro2name[m2t['ARO_accession']] = m2t['ARO_description']
              
          if 'ARO_category' in m2t:
              cats = []
              for cat, c2t in m2t['ARO_category'].items():
                  if 'category_aro_accession' in c2t:
                      cats.append(c2t['category_aro_accession'])
              aro2categories[m2t['ARO_accession']] = '|'.join(cats)
              
              
      db['ARO_description'] = db['ARO'].map(aro2name)
      db['ARO_category_accessions'] = db['ARO'].map(aro2categories)
      
      header = ['gene', 'CARD_short_name', 'ARO', 'ARO_description', 'ARO_category_accessions', 'target']
      db = db[header + [x for x in list(db.columns) if x not in header]]
      
      return db
    
  floc = "/genes.faa_vs_CARD.dm"
  Rdb = parse_card(floc, jloc = "card.json")

Where Rdb is a pandas DataFrame that looks like:

.. csv-table:: Rdb

  gene,CARD_short_name,ARO,ARO_description,ARO_category_accessions,target,percentID,alignment_length,mm,gaps,querry_start,querry_end,target_start,target_end,e-value,bit_score,protein_seq_accession
  AP010889.1_21,macB,3000535,MacB is an ATP-binding cassette (ABC) transpor...,0010001|0000006|0000000|3000159|0010000,gb|AAV85982.1|ARO:3000535|macB,34.6,231,137,5,1,227,3,223,7.170000e-33,123.0,AAV85982.1
  AP010889.1_53,lin,3004651,Listeria monocytogenes EGD-e lin gene for linc...,3000221|0000046|0000017|0001004,gb|AEO25219.1|ARO:3004651|lin,21.4,415,260,14,170,560,80,452,1.560000e-07,51.2,AEO25219.1
  AP010889.1_106,patA,3000024,PatA is an ABC transporter of Streptococcus pn...,0010001|0000036|3000662|0000001|3000159|0010000,gb|AAK76137.1|ARO:3000024|patA,25.9,228,149,7,80,298,344,560,1.560000e-14,70.9,AAK76137.1
  AP010889.1_107,bcrA,3002987,bcrA is an ABC transporter found in Bacillus l...,0010001|0000041|3000629|3000630|3000631|300005...,gb|AAA99504.1|ARO:3002987|bcrA,28.3,212,149,2,5,216,4,212,7.430000e-26,99.8,AAA99504.1
  AP010889.1_129,Abau_AbaF,3004573,Expression of abaF in E. coli resulted in incr...,0010002|0000025|3007149|3000159|0010000,gb|ABO11759.2|ARO:3004573|Abau_AbaF,31.5,435,280,7,24,453,4,425,7.750000e-72,230.0,ABO11759.2


Human milk oligosaccharide (HMO) Utilization genes
````````````````````````````````````````````````````
This is pretty niche, but it's something I (Matt Olm) am interested in. So here is how it can be done!::

  # Download the Supplemental Table S4 from here: https://data.mendeley.com/datasets/gc4d9h4x67/2

  wget https://data.mendeley.com/public-files/datasets/gc4d9h4x67/files/565528fe-585a-4f71-bb84-9f76625a872b/file_downloaded -O humann2_HMO_annotation.csv

  # Download the reference genome from here: https://www.ncbi.nlm.nih.gov/nuccore/CP001095.1

  Click “Send to:" -> “Coding Sequences” -> Format: “FASTA Protein” -> Rename to "Bifidobacterium_longum_subsp_infantis_ATCC_15697.NCBI.faa"

  # Pull the HMO genes using pullseq

  pullseq -i Bifidobacterium_longum_subsp_infantis_ATCC_15697.NCBI.faa -n HMO_list > Blon_HMO_genes.faa

  # Search against them

  diamond makedb --in Blon_HMO_genes.faa -d /LAB_DATA/DATABASES/HMO_ID/Blon_HMO_genes.faa.dmd

  diamond blastp -q genes.faa -d Blon_HMO_genes.faa.dmd.dmnd  -f 6 -e 0.0001 -k 1 -p 6 -o genes.faa_vs_HMO.b6

The following python code parses the resulting table

.. code-block:: python
  def parse_HMOs(floc, iloc):
      """
      v1.0 - 1/6/2023
      
      Parse HMOs 
      
      Returns:
          Hdb: DataFrame with HMO results
      """
      
      Hdb = pd.read_csv(iloc, sep=';')
      Hdb['target'] = [x.replace('_cds_', '_prot_').replace('lcl.', 'lcl|').strip() for x in Hdb['HMOgenes']]
      Hdb = Hdb[['target', 'Blon', 'Cluster']]

      h = ['gene', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',
          'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',
          'extra']
      db = pd.read_csv(floc, sep='\t', names=h)
      del db['extra']
      
      # Filter a bit
      db = db[(db['percentID'] >= 50)]
      
      # Merge
      Hdb = pd.merge(db, Hdb, how='left')
      
      # Re-order
      header = ['gene', 'Blon', 'Cluster', 'target']
      Hdb = Hdb[header + [x for x in list(Hdb.columns) if x not in header]]
      
      return Hdb
      
  floc = '/LAB_DATA/CURRENT/CURRENT_Metagenomics_PROJECTS/2022_Misame/gene_annotation/HMO/DeltaI_NewBifido.faa_vs_HMO.b6'
  iloc = '/LAB_DATA/DATABASES/HMO_ID/humann2_HMO_annotation.csv'

  Hdb = parse_HMOs(floc, iloc)

Where Hdb is a pandas DataFrame that looks like:

.. csv-table:: Hdb

  gene,Blon,Cluster,target,percentID,alignment_length,mm,gaps,querry_start,querry_end,target_start,target_end,e-value,bit_score
  AP010889.1_103,Blon_0104,Urease,lcl|CP001095.1_prot_ACJ51233.1_97,99.8,433,1,0,1,433,1,433,1.570002e-318,852.0
  AP010889.1_104,Blon_0105,Urease,lcl|CP001095.1_prot_ACJ51234.1_98,100.0,294,0,0,1,294,1,294,1.660000e-195,530.0
  AP010889.1_105,Blon_0106,Urease,lcl|CP001095.1_prot_ACJ51235.1_99,100.0,371,0,0,1,371,1,371,1.930000e-267,718.0
  AP010889.1_106,Blon_0107,Urease,lcl|CP001095.1_prot_ACJ51236.1_100,100.0,257,0,0,49,305,14,270,3.000000e-186,506.0
  AP010889.1_107,Blon_0108,Urease,lcl|CP001095.1_prot_ACJ51237.1_101,100.0,235,0,0,1,235,2,236,2.300000e-166,451.0