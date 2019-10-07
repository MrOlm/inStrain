Preparing input
===================

There are two main inputs to inStrain- a .fasta file and a mapping file. This describes some considerations to keep in mind when you use them

Preparing the .fasta file
----------------

The .fasta file that you map to can be a single genome, a metagenomic assembly, or a collection of genomes concatenated together.

Using a single genome .fasta file
++++++++++++++

If your .fasta file is a single genome, the main consideration is that it should be a good representitive genome for some organism in your sample. You can confirm this by looking at the breadth of coverage (see see :doc:`example_output` for more info). You can download a previously sequenced reference genome from an online database for this task, or you can assemble the genome out of the reads yourself (see next section)

** Maybe show that plot here about how genome similarity effects inStrain performance? **

Using a metagenomic assembly
+++++++++++++++

** CC has some plots here about the effect of assembly **

You can use ``inStrain genome_wide`` to add the results of binning to this as well

Using a collection of genomes
+++++++++++++

If you want to run inStrain on many genomes, it may make sense to concatonate them into a single .fasta file and map reads to that altogether. That way each read can only map to one genome. You can then convert the files to an individual genome basis using the command ``inStrain genome_wide``

Preparing the .bam file
----------------

Required to use paired reads. We use bowtie2 to map, but as long as it generates a valid .sam file, and mapper should do. CC has a lot of thoughts about this I think

Other
-----------

Use prodigal if you want genes
