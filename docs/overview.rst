Overview
========

inStrain is a program that characterizes the microdiveristy present in a metagenomic sample. The includes calculation of microdiveristy, clonality, and other metrics, analysis of linkage, analysis of SNPs present, ect. The typical use-case is to generate a .bam file by mapping metagenomic reads to a bacterial genome that is present in the metagenomic sample, and using inStrain to characterize the microdiveristy present.

Another common use-case is detailed strain comparisons that involves comparing the microdiveristy of two populations to see the extent to which they overlap. This allows the calculation of ANI values for extremely similar genomes (<99.999% average nucleotide identity).

To get started using it, see :doc:`quickstart`

For descriptions of what the modules can do, see :doc:`module_descriptions`

For an example of the output you can expect, see :doc:`example_output`

What is microdiveristy?
-----------------

It's like a read cloud. Here is a paragraph about it or something

How does inStrain work?
---------------

It's based on the assumption that read read-pair comes from a single individual DNA molecule, and thus a single member of the microbial community.

How does the compare function work?
--------------

You're essentially looking for overlap in the microdiveristies. Maybe I should draw a figure here.
