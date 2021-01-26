Installation
============

Installation
-------------------

InStrain is written in python. There are a number of ways that it can be installed.

pip
+++++++++++++++++

To install inStrain using the PyPi python repository, simply run ::

$ pip install instrain

That's it!

Pip is a great package with many options to change the installation parameters in various ways. For details, see `pip documentation <https://packaging.python.org/installing/>`_

bioconda
+++++++++++++++++

To inStrain inStrain from `bioconda <https://anaconda.org/bioconda/instrain>`_, run ::

$ conda install -c conda-forge -c bioconda -c defaults instrain

From source
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

Pre-built genome database
--------------------------

An established set of public genomes can be downloaded for your inStrain analysis at the `following link - https://doi.org/10.5281/zenodo.4441269 <https://doi.org/10.5281/zenodo.4441269>`_. See Tutorial #2 in :doc:`tutorial` for usage instructions.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4441269.svg

   :target: https://doi.org/10.5281/zenodo.4441269

Docker image
-------------------

A Docker image with inStrain and dependencies already installed is available on Docker Hub at `mattolm/instrain <https://hub.docker.com/repository/docker/mattolm/instrain>`_. This image also has a wrapper script in it to make it easier to use inStrain with AWS. See the `docker folder of the GitHub page <https://github.com/MrOlm/inStrain/tree/v1.3.0/docker>`_ for use instructions.