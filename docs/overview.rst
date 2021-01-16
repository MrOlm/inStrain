Glossary & FAQ
==========================

Glossary of terms used in inStrain
------------------------------------

.. note::
  This glossary is meant to give a conceptual overview of the terms used in inStrain. See :doc:`example_output` for explanations of specific output data.

.. glossary::

    ANI
      Average nucleotide identity. The average nucleotide distance between two genomes or .fasta files. If two genomes have a difference every 100 base-pairs, the ANI would be 99%

    conANI
      Consensus ANI - average nucleotide identity values calculated based on consensus sequences. This is commonly reported as "ANI" in other programs. Each position on the genome is represented by the most common allele (also referred to as the consensus allele), and minor alleles are ignored.

    popANI
      Population ANI - a new term to describe a unique type of ANI calculation performed by inStrain that considers both major and minor alleles. If two populations share any alleles at a loci, including minor alleles, it does not count as a difference when calculating popANI. It's easiest to describe with an example: consider a genomic position where the reference sequence is 'A' and 100 reads are mapped to the position. Of the 100 mapped reads, 60 have a 'C' and 40 have an 'A' at this position. In this example the reads share a minor allele with the reference genome at the position, but the consensus allele (most common allele) is different. Thus, this position **would** count as a difference in conANI calculations (because the consensus alleles are different) and **would not** count as a difference in popANI calculations (because the reference sequence is present as an allele in the reads). See :doc:`important_concepts` for examples.

    Representative genome
      Representative genomes (RGs) are genomes that are used to represent some taxa. For example you could have a series of representative genomes to represent each clade of E. coli (one genome for each clade), or you could have one representative genome for the entire species of E. coli (in that case it would be a Species Representative Genome (SRG)). The base unit of inStrain-based analysis is the representative genome, and they are usually generated using the program `dRep <https://drep.readthedocs.io/en/latest/>`_

    Species representative genome
      A Species Representative Genome (SRG) is a :term:`representative genome` that is used to represent an entire single microbial species.

    Genome database
      A collection of :term:`representative genome`s that are mapped to simultaneously (competitive mapping).

    nucleotide diversity
      A measurement of genetic diversity in a population (microdiversity). We measure nucleotide diversity using the method from Nei and Li 1979 (often referred to as 'pi' π in the population genetics world). InStrain calculates nucleotide diversity at every position along the genome, based on all reads, and averages values across genes / genomes. This metric is influenced by sequencing error, but within study error rates should be consistent and this effect is often minor compared to the extent of biological variation observed within samples. This metric is nice because it is not affected by coverage. The formula for calculating nucleotide diversity is the sum of the frequency of each base squared: 1 - [(frequency of A)^2 + (frequency of C)^2 + (frequency of G)^2 + (frequency of T)^2 ].

    microdiversity
      We use the term microdiversity to refer to intraspecific genetic variation, i.e. the genetic variation between cells within a microbial species.

    clonality
      The opposite of nucleotide diversity (1 - nucleotide diversity). A deprecated term used in older versions of the program.

    SNV
      Single nucleotide variant. A single nucleotide change that is present in a faction of a  population. Can also be described as a genomic loci with multiple alleles present. We identify and call SNVs using a simple model to distinguish them from errors, and more importantly in our experience, careful read mapping and filtering of paired reads to be assured that the variants (and the reads that contain them) are truly from the species being profiled, and not from another species in the metagenome (we call it 'mismapping' when this happens). Note that a SNV refers to genetic variation *within a read set*.

    SNS
      Single nucleotide substitution. A single nucleotide change that has a fixed difference between two populations. If the reference genome has a 'A' at some position, but all of the reads have a 'C' at that position, that would be a SNS (if half of the reads have an 'A' and half of the reads have a 'C', that would be an SNV).

    divergent site
      A position in the genome where either an SNV or SNS is present.

    SNP
      Single Nucleotide Polymorphism. In our experience this term means different things to different people, so we have tried to avoid using it entirely (instead referring to SNSs, SNVs, and divert sites).

    linkage
      A measure of how likely two divergent sites are to be inherited together. If two alleles are present on the same read, they are said to be "linked", meaning that they are found together on the same genome. Loci are said to be in "linkage disequilibrium" when the frequency of association of their different alleles is higher or lower than what would be expected if the loci were independent and associated randomly. In the context of microbial population genetics, linkage decay is often used as a way to detect recombination among members of a microbial population. InStrain uses the metrics r2 (r squared) and D' (D prime) to measure linkage.

    coverage
      A measure of sequencing depth. We calculate coverage as the average number of reads mapping to a region. If half the bases in a scaffold have 5 reads on them, and the other half have 10 reads, the coverage of the scaffold will be 7.5

    breadth
      A measure of how much of a region is covered by sequencing reads. Breadth is an important concept that is distinct from sequencing coverage, and gives you an approximation of how well the reference sequence you're using is represented by the reads. Calculated as the percentage of bases in a region that are covered by at least a single read. A breadth of 1 means that all bases in a region have at least one read covering them

    expected breadth
       The breadth that would be expected if reads are evenly distributed along the genome, given a specific coverage value. Based on the function breadth = 1 - e ^{-0.883  *  coverage}. This is useful to establish whether or not the scaffold is actually in the reads, or just a fraction of the scaffold. If your coverage is 10x, the expected breadth will be ~1. If your actual breadth is significantly lower then the expected breadth, this means that reads are mapping only to a specific region of your scaffold (transposon, prophage, etc.). See :doc:`important_concepts` for more info.

    relative abundance
      The percentage of total reads that map a particular entity. If a metagenome has 1,000,000 reads and 1,000 reads to a particular genome, that genome is at 0.1% relative abundance

    contig
      A contiguous sequence of DNA. Usually used as a reference sequence for mapping reads against. The terms contig and scaffold are used interchangeably by inStrain.

    scaffold
      A sequence of DNA that may have a string of "N"s in it representing a gap of unknown length. The terms contig and scaffold are used interchangeably by inStrain.

    iRep
      A measure of how fast a population was replicating at the time of DNA extraction. Based on comparing the sequencing coverage at the origin vs. terminus of replication, as described in `Brown et. al., Nature Biotechnology 2016 <http://dx.doi.org/10.1038/nbt.3704>`_

    mutation type
      Describes the impact of a nucleotide mutation on the amino acid sequence of the resulting protein. N = non-synonymous mutation (the encoded amino-acid changes due to the mutation). S = synonymous mutation (the encoded amino-acid does not change due to the mutation; should happen ~1/6 of the time by random chance due to codon redundancy). I = intergenic mutation. M = multi-allelic SNV with more than one change (rare).

    dN/dS
      A measure of whether the set of mutations in a gene are biased towards synonymous (S) or non-synonymous (N) mutations. dN/dS is calculated based on mutations relative to the reference genome. dN/dS > 1 means the bias is towards N mutations, indicating the gene is under active selection to mutate. dN/dS < 1 means the bias is towards S mutations, indicated the gene is under stabilizing selection to not mutate. dN/dS = 1 means that N and S mutations are at the rate expected by mutating positions randomly, potentially indicating the gene is non-functional.

    pN/pS
      Very similar to dN/dS, but calculated at positions with at least two alleles present rather than in relation to the reference genome.

    fasta file
      A file containing a DNA sequence. Details on this file format can be found on `wikipedia <https://en.wikipedia.org/wiki/FASTA_format>`_

    bam file
      A file containing metagenomic reads mapped to a DNA sequence. Very similar to a `.sam` file. Details can be found `online <https://samtools.github.io/hts-specs/SAMv1.pdf>`_

    scaffold-to-bin file
      A .text file with two columns separated by tabs, where the first column is the name of a scaffold and the second column is the name of the bin / genome the scaffold belongs to. Can be created using the script `parse_stb.py <https://github.com/MrOlm/drep/blob/master/helper_scripts/parse_stb.py>`_ that comes with the program ``dRep``  See :doc:`example_output` for more info

    genes file
      A file containing the nucleotide sequences of all genes to profile, as called by the program Prodigal. See :doc:`example_output` for more info

    mismapped read
      A read that is erroneously mapped to a genome. InStrain profiles a population by looking at the reads mapped to a genome. These reads are short, and sometimes reads that originated from one microbial population map to the representative genome of another (for example if they share homology). There are several techniques that can be used to reduce mismapping to the lowest extent possible.

    multi-mapped read
      A read that maps equally well to multiple different locations in the .fasta file. Most mapping software will randomly select one position to place multi-mapped reads. There are several techniques that can be used to reduce multi-mapped reads to the lowest extent possible, including increasing the minimum MAPQ cutoff to >2 (which will eliminate them entirely).

    inStrain profile
      An inStrain profile (aka IS_profile, IS, ISP) is created by running the ``inStrain profile`` command. It contains  all of the program's internal workings, cached data, and is where the output is stored. Additional commands can then be run on an IS_profile, for example to analyze genes, compare profiles, etc., and there is lots of nice cached data stored in it that can be accessed using python.

    null model
      The null model describes the probability that the number of true reads that support a variant base could be due to random mutation error, assuming Q30 score. The default false discovery rate with the null model is 1e-6 (one in a million).

    mm
      The maximum number of mismatches a read-pair can have to be considered in the metric being considered. Behind the scenes, inStrain actually calculates pretty much all metrics for every read pair mismatch level. That is, only including read pairs with 0 mismatches to the reference sequences, only including read pairs with >= 1 mis-match to the reference sequences, all the way up to the number of mismatches associated with the "PID" parameter. Most of the time when it then generates user-facing output, it uses the highest mm possible and deletes the column label. If you'd like access to information on the mm-level, see the section titled "Dealing with mm"

    mapQ score
      MapQ scores are a measure of how well a read aligns to a particular position. They are assigned to each read mapped by bowtie2, but the details of how they are generated are incredibly confusing (see the following `link <http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html>`_ for more information). MapQ scores of 0 and 1 have a special meaning: if a read maps equally well to multiple different locations on a .fasta file, it always gets a MapQ score of 0 or 1.


FAQ (Frequently asked questions)
---------------------------------------

How does inStrain compare to other bioinformatics tools for strains analysis?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

A major difference is inStrain's use of the popANI and conANI, which allow consideration of minor alleles when performing genomic comparisons. See :doc:`important_concepts` for more information.

What can inStrain do?
++++++++++++++++++++++++++++++

inStrain includes calculation of nucleotide diversity, calling SNPs (including non-synonymous and synonymous variants), reporting accurate coverage / breadth, and calculating linkage disequilibrium in the contexts of genomes, contigs, and individual genes.

inStrain also includes comparing the frequencies of fixed and segregating variants between sequenced populations with extremely high accuracy, out-performing other popular strain-resolved metagenomics programs.

The typical use-case is to generate a `.bam` file by mapping metagenomic reads to a bacterial genome that is present in the metagenomic sample, and using inStrain to characterize the microdiversity present.

Another common use-case is detailed strain comparisons that involve comparing the genetic diversity of two populations and calculating the extent to which they overlap. This allows for the calculation of population ANI values for extremely similar genomic populations (>99.999% average nucleotide identity).

.. seealso::
  :doc:`installation`
    To get started using the program
  :doc:`module_descriptions`
    For descriptions of what the modules can do
  :doc:`example_output`
    To view example output
  :doc:`user_manual`
    For information on how to prepare data for inStrain
  :doc:`important_concepts`
    For detailed information on how to make sure inStrain is running correctly

How does inStrain work?
++++++++++++++++++++++++++++++

The reasoning behind inStrain is that every sequencing read is derived from a single DNA molecule (and thus a single cell) in the original population of a given microbial species. During assembly, the consensus of these reads are assembled into contigs and these contigs are binned into genomes - but by returning to assess the variation in the reads that assembled into the contigs, we can characterize the genetic diversity of the population that contributed to the contigs and genomes.

The basic steps:

1. Map reads to a `.fasta` file to create a `.bam` file

2. Stringently filter mapped reads and calculate coverage and breadth

3. Calculate nucleotide diversity and SNVs

4. Calculate SNV linkage

5. Optional: calculate gene statistics and SNV function

6. Optional: compare SNVs between samples.

What is unique about the way that inStrain compares strains?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Most strain-resolved pipelines compare the dominant allele at each position. If you have two closely related strains A and B in sample 1, with B being at higher abundance, and two closely related strains A and C in sample 2, with C being at higher abundance, most strain comparison pipelines will in actuality compare strain B and C. This is because they work on the principle of finding the dominant strain in each sample and then comparing the dominant strains. InStrain, on the other hand, is able to identify the fact that A is present in both samples. This is because it doesn't just compare the dominant alleles, but compares all alleles in the two populations. See :doc:`module_descriptions` and :doc:`choosing_parameters` for more information.

What is a population?
++++++++++++++++++++++++++++++

To characterize intra-population genetic diversity, it stands to reason that you first require an adequate definition of "population". InStrain relies mainly on population definitions that are largely technically limited, but also coincide conveniently with possibly biological real microbial population constraints (see `Olm et. al. mSystems 2020 <https://msystems.asm.org/content/5/1/e00731-19>`_ and `Jain et. al. Nature Communications 2018 <https://www.nature.com/articles/s41467-018-07641-9>`_). Often, we dereplicate genomes from an environment at average nucleotide identities (ANI) from 95% to 99%, depending on the heterogeneity expected within each sample - lower ANIs might be preferred with more complex samples. We then assign reads to each genome's population by stringently requiring that combined read pairs for SNP calling be properly mapped pairs with an similarity to the consensus of at least 95% by default, so that the cell that the read pair came from was at least 95% similar to the average consensus genotype at that position. Within an environment, inStrain makes it possible to adjust these parameters as needed and builds plots which can be used to estimate the best cutoffs for each project.

What are inStrain's computational requirements?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The two computational resources to consider when running inStrain are the number of processes given (``-p``) and the amount of RAM on the computer (usually not adjustable unless using cloud-based computing). Using inStrain v1.3.3, running inStrain on a .bam file of moderate size (1 Gbp of less) will generally take less than an hour with 6 cores, and use about 8Gb of RAM. InStrain is designed to handle large .bam files as well. Running a huge .bam file (30 Gbp) with 32 cores, for example, will take ~2 hours and use about 128Gb of RAM. The more processes you give inStrain the longer it will run, but also the more RAM it will use. See :doc:'important_concepts` for information on reducing compute requirements.

What mapping software can be used to generate .bam files for inStrain?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Bowtie2 is a common one the works well, but any software that generates .bam files should work. Some mapping software modifies .fasta file headers during mapping (including the tool BBMap and SNAP). Include the flag ``--use_full_fasta_header`` when mapping with these programs to properly handle this.