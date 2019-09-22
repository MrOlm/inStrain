Quick Start
===========

The functionality of inStrain is broken up into modules. To see a list of the available modules, see :doc:`module_descriptions`.::

 $ inStrain -h

               ...::: inStrain v0.7.2 :::...

  Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

  Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
  Example: inStrain profile -h

   profile         -> Calculate microdiversity metrics and SNVs from a mapping. Must run this first to perform most other operations
   compare         -> Compare multiple inStrain profiles on a microdiversity level. Calculates popANI, coverage_overlap, and other things
   profile_genes   -> Calculate gene-level metrics on an inStrain profile
   plot            -> Make plots from an inStrain profile
   filter_reads    -> Commands related to filtering reads from .bam files
   other           -> Other miscellaneous operations (determine taxonomy / check dependencies)

profile
---------------

De-replication is the process of identifying groups of genomes that are the "same" in a genome set, and removing all but the "best" genome from each redundant set. How similar genomes need to be to be considered "same", how the "best" genome is chosen,  and other options can be adjusted (see :doc:`choosing_parameters`)

To de-replicate a set of genomes, run the following command::

 $ dRep dereplicate outout_directory -g path/to/genomes/*.fasta

This will automatically de-replicate the genome list and produce lots of information about it.

.. seealso::
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters


compare
-----------------

dRep is able to perform rapid genome comparisons for a group of genomes and visualize their relatedness. For example::

 $ dRep compare output_directory -g path/to/genomes/*.fasta

For help understanding the output, see :doc:`example_output`

To change the comparison parameters, see :doc:`choosing_parameters`

.. seealso::
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters
