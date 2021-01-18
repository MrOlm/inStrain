Acknowledgements
==========================

People
++++++++++++++++++++++++

InStrain was developed by `Matt Olm <mattolm@berkeley.edu>`_ and
`Alex Crits-Christoph <crits-christoph@berkeley.edu>`_ in the `Banfield Lab <https://geomicrobiology.berkeley.edu/>`_ at the University of California, Berkeley.

Special thanks to all those who have provided `feedback on GitHub <https://github.com/MrOlm/inStrain/issues>`_ and otherwise, especially early adopters Keith Bouma-Gregson, Ben Siranosian, Yue "Clare" Lou, Naïma Madi, and Antônio Pedro Camargo.

Many of the ideas in :doc:`important_concepts` were honed over many years during countless discussions with members of the `Banfield Lab <https://geomicrobiology.berkeley.edu/>`_ and the `Sonnenburg Lab <https://sonnenburglab.stanford.edu/>`_. Special thanks to Christopher Brown, Keith Bouma-Gregson, Yue "Clare" Lou, Spencer Diamond, Alex Thomas, Patrick West, Alex Jaffe, Bryan Merrill, Matt Carter, and Dylan Dahan.

Software
+++++++++++++++++++++++++

InStrain relies on several previously published programs and python modules to run - see `here <https://github.com/MrOlm/inStrain/blob/master/setup.py>`_ `and here <https://bioconda.github.io/recipes/instrain/README.html>`_ for a full list of dependencies. Of special importance are `samtools <http://www.htslib.org>`_ (the basis for parsing .bam files) and `coverM <https://github.com/wwood/CoverM>`_ (the heart of ``quick_profile``).

While not a direct dependency, the open-source program `anvi’o <http://merenlab.org/software/anvio/>`_ was used as significant inspiration for several implementation details, especially related to multiprocessing efficiency and memory management.

Citation
+++++++++++++++++++++++++

The manuscript describing inStrain is available in `Nature Biotechnology <https://doi.org/10.1038/s41587-020-00797-0>`_ and on `bioRxiv <https://www.biorxiv.org/content/10.1101/2020.01.22.915579v1>`_
and can be cited as follows::

    Olm, M.R., Crits-Christoph, A., Bouma-Gregson, K., Firek, B.A., Morowitz, M.J., Banfield, J.F., 2021. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. Nature Biotechnology. https://doi.org/10.1038/s41587-020-00797-0

