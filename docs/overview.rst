Overview
========

When you sequence any microbial genome(s), you sequence a population of cells (with the exception of single cell genomics). This population may be a nearly clonal population in a culture flask, a nearly clonal population in the real world, or a highly heterogenous population in the real world, but there is always real biological genetic hetereogeneity within that population - every cell does not have the same genotype at every single position. 

 *inStrain* **is a program for measuring, comparing, and interrogating the genetic heterogeneity of microbial populations ("intraspecific"; we can call it "microdiversity") in and between metagenomic samples.**
~~~~~~~~~~~~~~~~~~~~~~

*inStrain* includes calculation of nucleotide diversity, calling SNPs (including non-synonymous and synonymous variants), reporting accurate coverage / breadth, and calculating linkage disequilibrium in the contexts of genomes, contigs, and individual genes.

*inStrain* also includes comparing the frequencies of fixed and segregating variants between sequenced populations with extremely high accuracy, out-performing other popular "strain-resolved" metagenomics programs.

The typical use-case is to generate a `.bam` file by mapping metagenomic reads to a bacterial genome that is present in the metagenomic sample, and using inStrain to characterize the microdiversity present.

Another common use-case is detailed strain comparisons that involves comparing the genetic diversity of two populations and calculating the extent to which they overlap. This allows for the calculation of population ANI values for extremely similar genomic populations (<99.999% average nucleotide identity).

To get started using it, see :doc:`quickstart`

For descriptions of what the modules can do, see :doc:`module_descriptions`

For an example of the output you can expect, see :doc:`example_output`


When should I use *inStrain*?
-----------------

inStrain is intended to be used as a genome-resolved metagenomics approach. Genome-resolved metagenomics involves sequencing and then de novo assembly of the actual microbial genomes present in the sample(s) of interest. It is these microbial genomes, and not microbial genomes derived from reference databases, that we will then use as scaffolds to recruit reads from the sample. 

We don't recommend using reference genomes for strain-resolved inferences in metagenomes. This is because reference databases have usually poorly sampled the true extent of microbial diversity below the species level across many environments. Using even partially inaccurate references can lead to inaccurate conclusions about the genetic variation within your samples. 

inStrain can be run on individual microbial genomes assembled and binned from a metagenome, (recommended) or entire metagenomic assemblies at once. However, it is important to note that when run on entire metagenomic assemblies, the results must be interpreted in the context of each unique species in that community - this can be resolved with a scaffold to bin file after the fact. 

When should I probably not use *inStrain*?
---------------

When you have not assembled genomes from the metagenomic samples you are interrogating; When breadth and coverage of the consensus genome are low (best above 10x;95%); when you wish to compare populations that are <95% ANI with each other; when you are interested in community composition, not intra-population diversity. 

How does *inStrain* work?
---------------

The reasoning behind inStrain is that every sequencing read derived from a single DNA molecule (and likely a single cell) in the original population of a given microbial species. The weighted consensus of these reads was what was assembled into contigs, and was binned into genomes - but by returning to assess the variation in the reads that assembled into the contig, we can characterize the genetic diversity of the population that contributed to the contigs and genomes.

The basic steps:

1. Map reads to a BAM

2. Stringently filter mapped reads and calculate coverage and breadth

3. Calculate nucleotide diversity and SNPs

4. Calculate SNP linkage

5. Optional: calculate gene statistics and SNP function

6. Optional: compare SNPs between samples.

What are the metrics and terminology of *inStrain*?
--------------

**Community**: the collection of species in a metagenome, i.e. the species diversity of a microbiome.

**Population**: the collection of cells for each species in a metagenome, i.e. the genetic diversity of each species in a microbiome.

*InStrain* is for characterizing metagenomes at the population level, not at the community level. 

**SNP**: A SNP is a Single Nucleotide Polymorphism, a genetic variant of a single nucleotide change that some percentage of the cells that comprise a species population. We identify and call SNPs using a simple model to distinguish them from errors, and more importantly in our experience, careful read mapping and filtering of 300 bp (2x150 bp paired reads carefully evaluated as a pair) to be assured that the variants (and the reads that contain them) are truly from the species being profiled, and not from another species in the metagenome (we call it 'mismapping' when this happens). Note that a SNP refers to genetic variation *within a read set*.

**Microdiversity**: We use the term microdiversity to refer to intraspecific genetic variation, i.e. the genetic variation between cells within a microbial species. To measure this, we calculate a per-site nucleotide diversity of all reads - thus this metric is slightly influenced by sequencing error, but within study error rates should be consistent, and this effect is extremely minor compared to the extent of biological variation observed within samples. The metric of nucleotide diversity (often referred to as 'pi' in the population genetics world) is from Nei and Li 1979, calculated per site and then averaged across all sites.

**refSNP**: a genetic difference between the consensus of a read set and a reference genome. This is in contrast to SNPs, which are variants within a population being studied - reference SNPs are differences between the population you are studying (your reads) and the genome that you are mapping to. If you are mapping to a genome that was assembled from that sample, there will be very few to no refSNPs, because the consensus of that genome was built from the consensus of the reads in that sample. However, refSNPs are useful to track and understand cross-mapping, and we also use the percentage of refSNPs per read pair to filter read mappings.

**popANI**: calculated by `inStrain compare` function between two different inStrain profiles. This reflects the percentage of *fixed sites* in sample A that are also *fixed* in sample B. Sites that are segregating in either sample are ignored. 

**N SNP**: a polymorphic variant that changes the amino acid code of the protein encoded by the gene in which it resides; non-synonymous. 

**S SNP**: a polymoprhic variant that does not change the amino acid code of the protein encoded by the gene in which it resides; synonymous.

How does the compare function work?
--------------

`inStrain compare` calculates popANI metrics comparing a set of different inStrain profiles (sets of different sample reads, mapped to the same consenuss / reference genome). To use, we recommend assembly and binning of each sample, and then dereplication of genomes using the software dRep (https://drep.readthedocs.io/) at a high percent ANI, e.g. 96%-99%. Samples which contain multiple populations of the same dRep cluster (members of similar species or sub-species) can then be mapped back to the best genome from this dRep cluster, and then inStrain should be run on these dRep cluster genomes. 

**To restate: inStrain can only compare read profiles that have been mapped to the same reference genome**.

The metric that `inStrain compare` returns is a popANI between all inStrain profiles in the set. This reflects the percentage of *fixed sites* in sample A that are also *fixed* in sample B. Sites that are segregating in either sample are ignored. We find that popANI is highly specific for tracking near identical microbial populations, which should have popANI values of almost exactly 100%.
