Overview
========

When you sequence any microbial genome(s), you are always sequencing a population of cells (with the exception of single cell genomics). This population may be a nearly clonal population in a culture flask, a nearly clonal population in the real world, or a highly heterogenous population in the real world, but there is always real biological genetic hetereogeneity within that population - every cell does not have the same genotype at every single position. 

 *inStrain* **is a program for measuring, comparing, and interrogating that population's genetic heterogeneity ("intraspecific"; we can call it "microdiversity") in and between samples.**
~~~~~~~~~~~~~~~~~~~~~~

The includes calculation of nucleotide diversity, SNPs (including non-synonymous and synonymous variants), coverage / breadth, and linkage disequilibrium in the contexts of genomes, contigs, and individual genes.

This also includes comparing the frequencies of fixed and segregating variants between sequenced populations at extremely high accuracies, out-performing other popular "strain-resolved" metagenomics programs.

The typical use-case is to generate a `.bam` file by mapping metagenomic reads to a bacterial genome that is present in the metagenomic sample, and using inStrain to characterize the microdiversity present.

Another common use-case is detailed strain comparisons that involves comparing the microdiversity of two populations to see the extent to which they overlap. This allows for the calculation of ANI values for extremely similar genomic populations (<99.999% average nucleotide identity).

To get started using it, see :doc:`quickstart`

For descriptions of what the modules can do, see :doc:`module_descriptions`

For an example of the output you can expect, see :doc:`example_output`


When should I use *inStrain*?
-----------------

inStrain is intended to be used as a genome-resolved approach. Genome-resolved meta(genomics) involves sequencing and then de novo assembly of the actual microbial genomes present in the sample(s) of interest. It is these microbial genomes, and not microbial genomes derived from reference databases, that we will then use as scaffolds to recruit reads from the original sample. 

inStrain can be run on individual microbial genomes assembled and binned from a metagenome, (recommended) or entire metagenomic assemblies at once. However, it is important to note that when run on entire metagenomic assemblies, the results must be interpreted in the context of each unique species in that community - this can be resolved with a scaffold to bin file after the fact. 

How does *inStrain* work?
---------------

The reasoning behind inStrain is that every read came from a single DNA molecule (and likely a single cell) in the original species population. The weighted consensus of these reads and therefore of these cells was assembled into contigs, and binned into genomes - but by returning to assess the variation in the reads that assembled into the contig, we can characterize the genetic diversity of the population that contributed to the assembled contig.

Steps:
<insert steps>

What are the metrics and terminology of *inStrain*?
--------------

Community: the collection of species in a metagenome, i.e. the species diversity of a microbiome.

Population: the collection of cells for each species in a metagenome, i.e. the genetic diversity of each species in a microbiome.

InStrain is for characterizing metagenomes at the population level, not at the community level. 

<microdiversity, SNPs, linkage, N vs S, FST, etc.>

How does the compare function work?
--------------

You're essentially looking for overlap in the microdiveristies. Maybe I should draw a figure here.
