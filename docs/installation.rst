Installation and Quickstart
============

Installation
-----------

InStrain is written in python. There are a number of ways that is can be installed.

pip
+++++++++++++++++

To install inStrain using the PyPi python repository, simply run ::

$ pip install instrain

That's it!

Pip is a great package with many options to change the installation parameters in various ways. For details, see `pip documentation <https://packaging.python.org/installing/>`_

bioconda
+++++++++++++++++

To inStrain inStrain from `bioconda <https://anaconda.org/bioconda/instrain>`_, run ::

$ conda config --add channels bioconda; conda install instrain

from source
+++++++++++++++++

To install inStrain from the source code, run ::

  $ git clone https://github.com/MrOlm/instrain.git

  $ cd instrain

  $ pip install .

Dependencies
+++++++++++++++++

inStrain requires a few other programs to run. Not all dependencies are needed for all operations. There are a number of python
package dependencies, but those should install automatically when inStrain is installed using pip

**Essential**

* `samtools <http://www.htslib.org>`_ This is needed for pysam

**Optional**

* `coverM <https://github.com/wwood/CoverM>`_ This is needed for the quick_profile operation

* `Prodigal <https://github.com/hyattpd/Prodigal>`_ This is needed to profile on a gene by gene level

Docker image
-----------

A Docker image with inStrain and dependencies already installed is available on Docker Hub at `mattolm/instrain <https://hub.docker.com/repository/docker/mattolm/instrain>`_. This image also has a wrapper script in it to make it easier to use inStrain with AWS. See the `docker folder of the GitHub page <https://github.com/MrOlm/inStrain/tree/v1.3.0/docker>`_ for use instructions.

Quick Start
-----------

The functionality of inStrain is broken up into several core modules. For more details on these modules, see :doc:`module_descriptions`.::

  $ inStrain -h

                ...::: inStrain v1.0.0 :::...

  Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019

  Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
  Example: inStrain profile -h

    profile         -> Create an inStrain profile (microdiversity analysis) from a mapping.
    compare         -> Compare multiple inStrain profiles (popANI, coverage_overlap, etc.)
    profile_genes   -> Calculate gene-level metrics on an inStrain profile
    genome_wide     -> Calculate genome-level metrics on an inStrain profile
    quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
    filter_reads    -> Commands related to filtering reads from .bam files
    plot            -> Make figures from the results of "profile" or "compare"
    other           -> Other miscellaneous operations

Below is a list of brief descriptions of each of the modules. For more information see :doc:`module_descriptions`, for help understanding the output, see :doc:`example_output`, and to change the parameters see :doc:`choosing_parameters`

.. seealso::
  :doc:`module_descriptions`
    for more information on the modules
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters
  :doc:`preparing_input`
    for information on how to prepare data for inStrain

profile
+++++++++++++++++

inStrain profile is the main method of the program. It takes a `.fasta` file and a `.bam` file (consisting of reads mapping to the `.fasta` file) and runs a series of steps to characterize the microdiversity, SNPs, linkage, etc. Details on how to generate the mapping, how the profiling is done, explanations of the output, how to choose the parameters can be found at :doc:`preparing_input` and :doc:`module_descriptions`

To run inStrain on a mapping run the following command::

 $ inStrain profile .bam_file .fasta_file -o IS_output_name

compare
+++++++++++++++++

inStrain is able to compare multiple read mappings to the same .fasta file. Each mapping file must first be make into an inStrain profile using the above command. The coverage overlap and popANI between all pairs is calculated::

 $ inStrain compare -i IS_output_1 IS_output_2 IS_output_3

profile_genes
+++++++++++++++++

Once you've run `inStrain profile`, you can also calculate gene-wise microdiversity, coverage, and SNP functions using this command. It relies on having gene calls in the `.fna` format from the program prodigal::

 $ inStrain profile_genes -i IS_output -g called_genes.fna

genome_wide
+++++++++++++++++

This module is able to translate scaffold-level results to genome-level results. If the `.fasta` file you mapped to consists of a single genome, running this module on its own will average the results among all scaffolds. If the `.fasta` file you mapped to consists of several genomes, by providing a `scaffold to bin file` or a list of the individual `.fasta` files making up the combined `.fasta` file, you can get summary results for each individual genome. Running this module is also required before generating plots.

 $ inStrain genome_wide -i IS_output -s genome1.fasta genome2.fasta genome3.fasta

quick_profile
+++++++++++++++++

This auxiliary module  is merely a quick way to calculate the coverage and breadth using the blazingly fast program `coverM <https://github.com/wwood/CoverM>`_. This can be useful for quickly figuring out which scaffolds have any coverage, and then generating a list of these scaffolds to profile with inStrain profile, making it run faster::

 $ inStrain quick_profile -b .bam_file -f .fasta_file -s scaffold_to_bin_file -o output_name

filter_reads
+++++++++++++++++

This auxiliary module lets you do various tasks to filter and/or characterize a mapping file, and then generate a new mapping file with those filters applied::

 $ inStrain filter_reads .bam_file .fasta_file -g new_sam_file_location

plot
+++++++++++++++++

This method makes a number of plots from an inStrain object. It is required that you run `genome_wide` first before running this module::

 $ inStrain plot -i IS_output

other
+++++++++++++++++

This module lets you do random small things, like convert IS_profile objects that are in an old format to the newest format.
