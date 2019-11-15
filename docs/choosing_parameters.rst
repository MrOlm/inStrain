Choosing Parameters
===================

There are a number of important considerations when running inStrain. Here is some theory and data about how to make inStrain work best

Reference genome selection
------------------

inStrain relies on mapping reads from a sample to a reference genome. How similar the reference genome is to the reads, and the minimum read ANI threshold that you set, are very important and will determine much of what you get out of inStrain.

Below are a series of plots made by introducing a known number of mutations into an E. coli genome, simulating reads from these mutated genomes with known ANI differences from the original reference genome, mapping the synthetic reads back to the original reference genome, and running inStrain.

.. figure:: images/Fig1.png
  :width: 400px
  :align: center

In the above plot, inStrain was run with a minimum read ANI (inStrain profile parameter -l or --filter_cutoff) of 0.99. The reported genome breadth is reported on the y-axis.
