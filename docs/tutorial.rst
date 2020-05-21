Tutorial
===================

The following tutorial goes through an example run of inStrain. You can follow along with your own data, or use a small set of reads that are included in the inStrain install for testing. They can be found in the folder `test/test_data/` of your install folder, or can be downloaded from the inStrain source code at `this link on GitHub
<https://github.com/MrOlm/inStrain/tree/master/test/test_data>`_. The only files that you'll need for this tutorial are forward and reverse metagenomic reads (`N5_271_010G1.R1.fastq.gz` and `N5_271_010G1.R2.fastq.gz`) and a .fasta file to map to (`N5_271_010G1_scaffold_min1000.fa`). In case you're curious, these metagenomic reads come from a premature infant fecal sample.

.. seealso::
  :doc:`overview`
    To get started using the program
  :doc:`program_documentation`
    For descriptions of what the modules can do and information on how to prepare data for inStrain
  :doc:`example_output`
    To view example output
  :doc:`Advanced_use`
    For detailed information on how to rationally adjust inStrain parameters

Preparing .bam and .fasta files
----------------

After downloading the genome file that you would like to profile (.fasta file) and at least one set of paired reads, the first thing to do is to map the reads to the .fasta file in order to generate a .bam file.

When this mapping is performed it is important that you map to all genomes simultaneously, so the first thing to do is to combine all of the genomes that you'd like to map into a single .fasta file. In our case our .fasta file already has all of the genomes that we'd like to profile within it, but if you did want to profile a number of different genomes, you could combine them using a command like this ::

 $  cat raw_data/S2_002_005G1_phage_Clostridioides_difficile.fasta raw_data/S2_018_020G1_bacteria_Clostridioides_difficile.fasta > allGenomes_v1.fasta

Next we must map our reads to this .fasta file to create .bam files. In this tutorial we will use the mapping program Bowtie 2 ::

 $ mkdir bt2

 $ bowtie2-build ~/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa bt2/N5_271_010G1_scaffold_min1000.fa

 $ bowtie2 -p 6 -x bt2/N5_271_010G1_scaffold_min1000.fa -1 ~/Programs/inStrain/test/test_data/N5_271_010G1.R1.fastq.gz -2 ~/Programs/inStrain/test/test_data/N5_271_010G1.R2.fastq.gz > N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam

At this point we  have generated a .sam file, the precursor to .bam files. Lets make sure it's there and not empty ::

 $ ls -lht

 total 34944
 -rw-r--r--  1 mattolm  staff    16M Jan 23 11:56 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam
 drwxr-xr-x  8 mattolm  staff   256B Jan 23 11:54 bt2/

Perfect. At this point we could convert the .sam file to a sorted and indexed .bam file, but since inStrain can do that for us automatically we won't bother.

Preparing genes file
---------------

If we want inStrain to do gene-level profiling we need to give it a list of genes to profile. **Note - this is an optional step that is not required for inStrain to work in general, but without this you will not get gene-level profiles**

We will profile our genes using the program prodigal, which can be run using the following example command ::

 $ prodigal -i ~/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa -d N5_271_010G1_scaffold_min1000.fa.genes.fna

Preparing for genome-level characterization
---------------

In the step above ("Preparing .bam and .fasta files"), we combined all of our genomes into a single .fasta file for mapping. However we likely want to profile the microdiversity of the individual genomes in that .fasta file. In order to do that we need to tell inStrain which scaffolds belong to which genomes.

There are two ways of providing this information. One is to give inStrain a list of the .fasta files that went into making the concatenated .fasta file. The other is to provide inStrain with a "scaffold to bin" file, which lists the genome assignment of each scaffold in a tab-delimited file. In this case we're going to use the scaffold to bin file provided by inStrain (called "N5_271_010G1.maxbin2.stb"). Here's what it looks like ::

  $ head ~/Programs/inStrain/test/test_data/N5_271_010G1.maxbin2.stb
  N5_271_010G1_scaffold_0 	 maxbin2.maxbin.001.fasta
  N5_271_010G1_scaffold_1 	 maxbin2.maxbin.001.fasta
  N5_271_010G1_scaffold_2 	 maxbin2.maxbin.001.fasta
  N5_271_010G1_scaffold_3 	 maxbin2.maxbin.001.fasta
  N5_271_010G1_scaffold_4 	 maxbin2.maxbin.001.fasta

Running inStrain profile
--------------

Now that we've gotten everything set up, it's time to run inStrain. To see all of the options, run ::

 $ inStrain -h

A long list of arguments and options will show up. For more details on what these do, see :doc:`program_documentation`. The **only** arguments that are absolutely required, however, are a .sam or .bam mapping file, and the .fasta file that the mapping file is mapped to.

.. note::
  In this case we're going to have inStrain profile the mapping, call genes, make the results genome wide, and plot the results all in one command. It is possible to do these all as separate steps, however, using the subcommands "inStrain profile", "inStrain profile_genes", "inStrain genome_wide", and "inStrain plot". See :doc:`program_documentation` for more information.

Using all of the files we generated above, here is going to be our inStrain command ::

 $ inStrain profile N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam ~/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa -o N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS -p 6 -g N5_271_010G1_scaffold_min1000.fa.genes.fna -s ~/Programs/inStrain/test/test_data/N5_271_010G1.maxbin2.stb

You should see the following as inStrain runs (should only take a few minutes) ::

  You gave me a sam- I'm going to make it a .bam now
  Converting N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam to N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam
  samtools view -S -b N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sam > N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam
  sorting N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam
  samtools sort N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam -o N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam
  Indexing N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam
  samtools index N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam.bai
  ***************************************************
      ..:: inStrain profile Step 1. Filter reads ::..
  ***************************************************

  Getting read pairs: 100%|██████████████████████████████████████████████████████████| 178/178 [00:00<00:00, 715.57it/s]
  Making read report
  /Users/mattolm/.pyenv2/versions/3.6.9/lib/python3.6/site-packages/numpy/core/fromnumeric.py:3335: RuntimeWarning: Mean of empty slice.
    out=out, **kwargs)
  /Users/mattolm/.pyenv2/versions/3.6.9/lib/python3.6/site-packages/numpy/core/_methods.py:161: RuntimeWarning: invalid value encountered in double_scalars
    ret = ret.dtype.type(ret / rcount)
  Filtering reads
  1,727 read pairs remain after filtering
  ***************************************************
  .:: inStrain profile Step 2. Profile scaffolds ::..
  ***************************************************

  Profiling scaffolds: 100%|████████████████████████████████████████████████████████████| 23/23 [00:06<00:00,  3.44it/s]
  Storing output
  ***************************************************
    .:: inStrain profile Step 3. Profile genes ::..
  ***************************************************

  20.67703568161025% of the input 1093 genes were marked as incomplete
  161 scaffolds with genes, 169 in the IS, 153 to compare
  Running gene-level calculations on scaffolds: 100%|█████████████████████████████████| 153/153 [00:18<00:00,  8.16it/s]
  ***************************************************
  .:: inStrain profile Step 4. Make genome-wide ::..
  ***************************************************

  Scaffold to bin was made using .stb file
  85.66% of scaffolds have a genome
  93.82% of scaffolds have a genome
  ***************************************************
   .:: inStrain profile Step 5. Generate plots ::..
  ***************************************************

  making plots 1, 2, 3, 4, 5, 6, 7, 8, 9
  85.66% of scaffolds have a genome
  Plotting plot 1
  Plotting plot 2
  85.66% of scaffolds have a genome
  Plotting plot 3
  57.37% of scaffolds have a genome
  Plotting plot 4
  97.33% of scaffolds have a genome
  Plotting plot 5
  93.82% of scaffolds have a genome
  Plotting plot 6
  Plotting plot 7
  97.33% of scaffolds have a genome
  Plotting plot 8
  93.96% of scaffolds have a genome
  Plotting plot 9
  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      ..:: inStrain profile finished ::..

  Output tables........ /Users/mattolm/Programs/testing_house/tutorial/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/output/
  Figures.............. /Users/mattolm/Programs/testing_house/tutorial/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/figures/

  See documentation for output descriptions - https://instrain.readthedocs.io/en/latest/

  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

The last note shows you where the plots and figures have been made. Here's a list of the files that you should see ::

  $ ls -lht N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/output/
  total 512
  -rw-r--r--  1 mattolm  staff   545B Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_genomeWide_mapping_info.tsv
  -rw-r--r--  1 mattolm  staff   602B Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_genomeWide_scaffold_info.tsv
  -rw-r--r--  1 mattolm  staff    25K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_SNP_mutation_types.tsv
  -rw-r--r--  1 mattolm  staff   125K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_gene_info.tsv
  -rw-r--r--  1 mattolm  staff    19K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_mapping_info.tsv
  -rw-r--r--  1 mattolm  staff    14K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_linkage.tsv
  -rw-r--r--  1 mattolm  staff    26K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_SNVs.tsv
  mattolm@Matts-MacBook-Pro-3:~/Programs/testing_house/tutorial$ caffold_min1000.fa-vs-N5_271_010G1.IS_scaffold_info.tsv

  $ ls -lht N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/figures
  total 7792
  -rw-r--r--  1 mattolm  staff   432K Jan 23 15:17 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_GeneHistogram_plot.pdf
  -rw-r--r--  1 mattolm  staff   422K Jan 23 15:17 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_LinkageDecay_types_plot.pdf
  -rw-r--r--  1 mattolm  staff   448K Jan 23 15:17 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_ScaffoldInspection_plot.pdf
  -rw-r--r--  1 mattolm  staff   419K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_ReadFiltering_plot.pdf
  -rw-r--r--  1 mattolm  staff   421K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_LinkageDecay_plot.pdf
  -rw-r--r--  1 mattolm  staff   420K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_MajorAllele_frequency_plot.pdf
  -rw-r--r--  1 mattolm  staff   419K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_readANI_distribution.pdf
  -rw-r--r--  1 mattolm  staff   443K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_genomeWide_microdiveristy_metrics.pdf
  -rw-r--r--  1 mattolm  staff   419K Jan 23 15:16 N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS_CoverageAndBreadth_vs_readMismatch.pdf

For help interpreting these output files, see :doc:`example_output`

inStrain compare
-----------

inStrain compare allows you to compare genomes that have been profiled by multiple mappings. To compare a genome in multiple samples, you must first map reads from multiple samples to the **same** .fasta file, then run run `inStrain profile on each mapping.

In this tutorial we profiled reads mapped to the .fasta file "N5_271_010G1_scaffold_min1000.fa". Provided in the inStrain test_data folder (<https://github.com/MrOlm/inStrain/tree/master/test/test_data>) is also a different set of reads mapped to the same .fasta file. We've also already run inStrain on this mapping for you! The resulting inStrain profile is the folder `N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.IS/`

To compare these inStrain profiles we will use the following command ::

 $ inStrain compare -i N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/ ~/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.IS/ -o N5_271_010G1_scaffold_min1000.fa.IS.COMPARE -p 6

  Loading N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/
  Loading /Users/mattolm/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.IS/
  Warning! Your inStrain folder is from version 1.0.0, while the installed version is 1.2.1.
  If you experience weird behavior, this might be why
  169 of 178 scaffolds are in at least 2 samples
  Profiling scaffolds: 100%|█████████████████████████████████████████████████████| 169/169 [00:22<00:00,  7.38it/s]

You should now have the following output file created ::

  $ ls -lht N5_271_010G1_scaffold_min1000.fa.IS.COMPARE/output/
  total 64
  -rw-r--r--  1 mattolm  staff    30K Jan 23 15:20 N5_271_010G1_scaffold_min1000.fa.IS.COMPARE_comparisonsTable.tsv

This file shows the comparison values between scaffolds, however. To make these on the genome level, we can run `inStrain genome_wide` ::

  $ inStrain genome_wide -i N5_271_010G1_scaffold_min1000.fa.IS.COMPARE/ -s ~/Programs/inStrain/test/test_data/N5_271_010G1.maxbin2.stb
  Scaffold to bin was made using .stb file
  89.62% of scaffolds have a genome

Now we should also have a table that compares these genomes on the genome level ::

  $ ls -lht N5_271_010G1_scaffold_min1000.fa.IS.COMPARE/output/
  total 72
  -rw-r--r--  1 mattolm  staff   556B Jan 23 15:23 N5_271_010G1_scaffold_min1000.fa.IS.COMPARE_genomeWide_compare.tsv
  -rw-r--r--  1 mattolm  staff    30K Jan 23 15:20 N5_271_010G1_scaffold_min1000.fa.IS.COMPARE_comparisonsTable.tsv

Finally, we can also plot these results using the `inStrain plot` function ::

  $ inStrain plot -i N5_271_010G1_scaffold_min1000.fa.IS.COMPARE/
  making plots 10
  89.62% of scaffolds have a genome
  Plotting plot 10
  Done!

This should make a figure in the figures folder ::

  $ ls -lht N5_271_010G1_scaffold_min1000.fa.IS.COMPARE/figures/
  total 936
  -rw-r--r--  1 mattolm  staff   419K Jan 23 15:25 N5_271_010G1_scaffold_min1000.fa.IS.COMPARE_inStrainCompare_dendrograms.pdf

As before, for help interpreting this output see :doc:`example_output` .
