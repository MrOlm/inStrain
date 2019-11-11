Preparing input
===================

There are two simple inputs to inStrain: a .FASTA file and a mapping file in .BAM format. A third, a prodigal `.FAA` file, can be used in later steps  Here we go over some considerations involved in choosing these inputs.

Preparing the .fasta file
----------------

Choosing what FASTA you will use (consensus / reference genomes) is extremely important and will affect the interpetation of your *inStrain* results. Below we describe the three most common strategies. 

Please note that the FASTA file must always be the same as, or a subset of, the FASTA file used to create the BAM, i.e. the FASTA file that reads were mapped to.

Using a collection of genomes (recommended)
+++++++++++++

The recommended workflow for running inStrain:

1. Assemble genomes out of each sample collected from the environment. Recommended software: IDBA_UD, MEGAHIT, metaSPADES.

2. Bin genomes out of each assembly using differential coverage binning. Recommended software: Bowtie2 (for mapping), MetaBAT, CONCOCT, DasTOOL.

3. Dereplicate the entire set of genomes across all samples from the environment at 97-99% identity, filter genome set draft quality genomes. Recommended software: dRep, checkM.

4. Create a bowtie2 index of the representative genomes from this dereplicated set and map reads to this set from each sample: Recommended software: Bowtie2

5. Profile the mapping BAMs using inStrain, passing FASTAs that are subsetted genomes of interest from the representative genome set.  

The most important aspect of this workflow is to **map to many genomes at once**. Mapping to just one genome at a time is highly discouraged, because this encourages mismapped reads from other genomes to be recruited by this genome. By including many (dereplicated) genomes in your bowtie2 index, you will be able to far more accurately filter out mismapped reads and reduce false positive SNPs.

Using a single genome .FASTA file
++++++++++++++

If your .fasta file is a single genome, the main consideration is that it should be a good representitive genome for some organism in your sample. Ideally, it was assembled directly from that sample, isolated from that sample, or you have some other evidence that this genome is highly representation of a species in that sample. Regardless, you should check your `inStrain plot` output and `scaffold_info.tsv` output file to be sure that your inStrain run had decent coverage and breadth of coverage of the genome that you use before attempting to interpret the results. 

Remember your FASTA only has to be a subset of the FASTA used to create a BAM. You can create a BAM with all dereplicated genomes from your environment, but then just pass a FASTA for only the genomes of particular interest. This approach is recommended as opposed to creating a BAM for just each genome, as it reduces mismapping.

Using a metagenomic assembly
+++++++++++++++

You can also pass *inStrain* an entire metagenomic assembly from a sample, including either binned or unbinned contigs. In this case, the output inStrain profile will include population information for each contig in the set. To then break it down by microbial genome / species, You can use ``inStrain genome_wide`` including a scaffold to bin file to generate results by genome.

Preparing the .bam file
----------------

*inStrain* requires paired-end Illumina read sequencing. We recommend using Bowtie2 to map your reads to your genome, and then convert the SAM created to an indexed, sorted bam. You can do this in samtools with the commands:

```
samtools view -S -b sample.sam > sample.bam

samtools sort sample.bam -o sample.sorted.bam
```

Preparing the prodigal `.faa` genes file for gene-level profiling
-----------

You can run prodigal on your FASTA file to generate the `.FAA` file that `inStrain profile_genes` requires.

Example:
`prodigal -i assembly.fasta -a proteins.faa`
