Quick Start
===========

The functionality of inStrain is broken up into modules. To see a list of the available modules, see :doc:`module_descriptions`.::

 $ inStrain -h

                ...::: inStrain v0.8.0 :::...

  Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

  Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
  Example: inStrain profile -h

    profile         -> Calculate microdiversity metrics and SNVs from a mapping. Must run this first to perform most other operations
    compare         -> Compare multiple inStrain profiles on a microdiversity level. Calculates popANI, coverage_overlap, and other things
    profile_genes   -> Calculate gene-level metrics on an inStrain profile
    quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
    filter_reads    -> Commands related to filtering reads from .bam files
    other           -> Other miscellaneous operations

Below is a list of brief descriptions of each of the modules. For more information see :doc:`module_descriptions`, for help understanding the output see :doc:`example_output`, to change the comparison parameters see :doc:`choosing_parameters`

.. seealso::
  :doc:`module_descriptions`
    for more information on the modules
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters

profile
---------------

inStrain profile is the main part of this program. It takes .bam file, consisting of reads mapping to a .fasta file, and runs a series of steps to characterize the microdiversity present. Details on how to generate the mapping, how the profiling is done, explanations of the output, how to choose the parameters can be found at :doc:`choosing_parameters` and :doc:`module_descriptions`

To run inStrain on a mapping run the following command::

 $ inStrain profile .bam_file .fasta_file -o IS_output_name

compare
-----------------

inStrain is able to compare multiple read mappings to the same .fasta file. Each mapping file must first be make into an inStrain profile using the above command. The coverage overlap and popANI between all pairs is calculated::

 $ inStrain compare -i IS_output_1 IS_output_2 IS_output_3

profile_genes
-----------------

Once you've run `inStrain profile`, you can also calculate gene-wise clonality, coverage, SNV identitiy, ect. using this command. It relies on having gene calls in the `.fna` format from the program prodigal.::

 $ inStrain profile_genes -i IS_output -g called_genes.fna
