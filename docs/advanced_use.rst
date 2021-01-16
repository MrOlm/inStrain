Raw data access and API
==============================

API for accessing raw data
----------------------------

inStrain stores much more data than is shown in the output folder. It is kept in the ``raw_data`` folder, and is mostly stored in compressed formats. This data can be easily accessed using python, as described below.

To access the data, you first make an SNVprofile object of the inStrain output profile, and then you access data from that object. For example, the following code accessed the raw SNP table ::

  import inStrain
  import inStain.SNVprofile

  IS = inStain.SNVprofile.SNVprofile(``/home/mattolm/inStrainOutputTest/``)
  raw_snps = IS.get('raw_snp_table')


You can use the example above (``IS.get()``) to access any of the raw data described in the following section. There are also another special things that are accessed in other ways, as described in the section "Accessing other data"

Basics of raw_data
++++++++++++++++++++

A typical run of inStrain will yield a folder titled "raw_data", with lots of individual files in it. The specifics of what files are in there depend on how inStrain was run, and whether or not additional commands were run as well (like profile_genes).

There will always be a file titled "attributes.tsv". This describes some basic information about each item in the raw data. Here's an example

.. csv-table:: attributes.tsv

    name,value,type,description
    location,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test,value,Location of SNVprofile object
    version,1.3.3,value,Version of inStrain
    Rdic,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/Rdic.json,dictionary,Read pair -> mismatches
    mapping_info,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/mapping_info.csv.gz,pandas,Report on reads
    fasta_loc,/Users/mattolm/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa,value,Location of .fasta file used during profile
    scaffold2length,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/scaffold2length.json,dictionary,Dictionary of scaffold 2 length
    object_type,profile,value,Type of SNVprofile (profile or compare)
    bam_loc,/Users/mattolm/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam,value,Location of .bam file
    scaffold_list,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/scaffold_list.txt,list,1d list of scaffolds that were profiled
    raw_linkage_table,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/raw_linkage_table.csv.gz,pandas,Raw table of linkage information
    raw_snp_table,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/raw_snp_table.csv.gz,pandas,Contains raw SNP information on a mm level
    cumulative_scaffold_table,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/cumulative_scaffold_table.csv.gz,pandas,Cumulative coverage on mm level. Formerly scaffoldTable.csv
    cumulative_snv_table,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/cumulative_snv_table.csv.gz,pandas,Cumulative SNP on mm level. Formerly snpLocations.pickle
    scaffold_2_mm_2_read_2_snvs,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/scaffold_2_mm_2_read_2_snvs.pickle,pickle,crazy nonsense needed for linkage
    genes_coverage,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/genes_coverage.csv.gz,pandas,Coverage of individual genes
    genes_clonality,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/genes_clonality.csv.gz,pandas,Clonality of individual genes
    genes_SNP_count,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/genes_SNP_count.csv.gz,pandas,SNP density and counts of individual genes
    SNP_mutation_types,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/SNP_mutation_types.csv.gz,pandas,The mutation types of SNPs
    covT,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/covT.hd5,special,Scaffold -> mm -> position based coverage
    clonT,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/clonT.hd5,special,Scaffold -> mm -> position based clonality
    clonTR,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/clonTR.hd5,special,Scaffold -> mm -> rarefied position based clonality
    genes_fileloc,/Users/mattolm/Programs/inStrain/test/test_data/N5_271_010G1_scaffold_min1000.fa.genes.fna,value,Location of genes file that was used to call genes
    genes_table,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/genes_table.csv.gz,pandas,Location of genes in the associated genes_file
    scaffold2bin,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/scaffold2bin.json,dictionary,Dictionary of scaffold 2 bin
    bin2length,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/bin2length.json,dictionary,Dictionary of bin 2 total length
    genome_level_info,/Users/mattolm/Programs/inStrain/test/test_backend/testdir/test/raw_data/genome_level_info.csv.gz,pandas,Table of genome-level information

This is what the columns correspond to:

name
  The name of the data. This is the name that you put into ``IS.get()`` to have inStrain retrieve the data for you. See the section "Accessing raw data" for an example.

value
  This lists the path to where the data is located within the raw_data folder. If the type of data is a value, than this just lists the value

type
  This describes how the data is stored. Value = the data is whatever is listed under value; list = a python list; numpy = a numpy array; dictionary = a python dictionary; pandas = a pandas dataframe; pickle = a piece of data that's stored as a python pickle object; special = a piece of data that is stored in a special way that inStrain knows how to de-compress

description
  A one-sentence description of what's in the data.

.. warning::

  Many of these pieces of raw data have the column "mm" in them, which means that things are calculated at every possible read mismatch level. This is often not what you want. See the section "Dealing with mm" for more information.

Accessing other data
++++++++++++++++++++++++++++++

In addition to the raw_data described above, there are a couple of other things that inStrain can make for you. You access these from methods that run on the IS object itself, instead of using the ``get`` method. For example::

  import inStrain
  import inStain.SNVprofile

  IS = inStain.SNVprofile.SNVprofile(``/home/mattolm/inStrainOutputTest/``)
  coverage_table = IS.get_raw_coverage_table()

The following methods work like that:

get_nonredundant_scaffold_table()
  Get a scaffold table with just one line per scaffold, not multiple mms

get_nonredundant_linkage_table()
  Get a linkage table with just one line per scaffold, not multiple mms

get_nonredundant_snv_table()
  Get a SNP table with just one line per scaffold, not multiple mms

get_clonality_table()
  Get a raw clonality table, listing the clonality of each position. Pass `nonredundant=False` to keep multiple mms

Dealing with "mm"
++++++++++++++++++++

Behind the scenes, inStrain actually calculates pretty much all metrics for every read pair mismatch level. That is, only including read pairs with 0 mis-match to the reference sequences, only including read pairs with >= 1 mis-match to the reference sequences, all the way up to the number of mismatches associated with the "PID" parameter.

For most of the output that inStrain makes in the output folder, it removes the "mm" column and just gives the results for the maximum number of mismatches. However, it's often helpful to explore other mismatches levels, to see how parameters vary with more or less stringent mappings. Much of the data stored in "read_data" is on the mismatch level. Here's an example of what the looks like (this is the cumulative_scaffold_table)::

  ,scaffold,length,breadth,coverage,coverage_median,coverage_std,bases_w_0_coverage,mean_clonality,median_clonality,unmaskedBreadth,SNPs,breadth_expected,ANI,mm
  0,N5_271_010G1_scaffold_102,1144,0.9353146853146853,5.106643356643357,5,2.932067325774674,74,1.0,1.0,0.6145104895104895,0,0.9889923642060382,1.0,0
  1,N5_271_010G1_scaffold_102,1144,0.9353146853146853,6.421328671328672,6,4.005996333777764,74,0.9992001028104149,1.0,0.6748251748251748,0,0.9965522492489882,1.0,1
  2,N5_271_010G1_scaffold_102,1144,0.9423076923076923,7.3627622377622375,7,4.2747074564903285,66,0.9993874800638958,1.0,0.7928321678321678,0,0.998498542620078,1.0,2
  3,N5_271_010G1_scaffold_102,1144,0.9423076923076923,7.859265734265734,8,4.748789115369562,66,0.9992251555869703,1.0,0.7928321678321678,0,0.9990314705263914,1.0,3
  4,N5_271_010G1_scaffold_102,1144,0.9423076923076923,8.017482517482517,8,4.952541407151938,66,0.9992251555869703,1.0,0.7928321678321678,0,0.9991577528529144,1.0,4
  5,N5_271_010G1_scaffold_102,1144,0.9458041958041958,8.271853146853147,8,4.9911156795536105,62,0.9992512780077317,1.0,0.8024475524475524,0,0.9993271891539499,1.0,7

As you can see, the same scaffold is shown multiple times, and the last column is ``mm``. At the row with mm = 0, you can see what the stats are when only considering reads that perfectly map to the reference sequence. As the mm goes higher, so do stats like coverage and breadth, as you now allow reads with more mismatches to count in the generation of these stats. In order to convert this files to what is provided in the output folder, the following code is run::

  import inStrain
  import inStain.SNVprofile

  IS = inStain.SNVprofile.SNVprofile(``/home/mattolm/inStrainOutputTest/``)
  scdb = IS.get('cumulative_scaffold_table')
  ScaffDb = scdb.sort_values('mm')\
              .drop_duplicates(subset=['scaffold'], keep='last')\
              .sort_index().drop(columns=['mm'])

The last line looks complicated, but it's very simple what is going on. First, you sort the database by ``mm``, with the lowest mms at the top. Next, for each scaffold, you only keep the row with the lowest mm. That's done using the ``drop_duplicates(subset=['scaffold'], keep='last')`` command. Finally, you re-sort the DataFrame to the original order, and remove the ``mm`` column. In the above example, this would mean that the only row that would survive would be where mm = 7, because that's the bottom row for that scaffold.

You can of course subset to any level of mismatch by modifying the above code slightly. For example, to generate this table only using reads with <=5 mismatches, you could use the following code::

  import inStrain
  import inStain.SNVprofile

  IS = inStain.SNVprofile.SNVprofile(``/home/mattolm/inStrainOutputTest/``)
  scdb = IS.get('cumulative_scaffold_table')
  scdb = scdb[scdb['mm'] <= 5]
  ScaffDb = scdb.sort_values('mm')\
              .drop_duplicates(subset=['scaffold'], keep='last')\
              .sort_index().drop(columns=['mm'])

.. warning::

  You usually do not want to subset these DataFrames using something like ``scdb = scdb[scdb['mm'] == 5]``. That's because if there are no reads that have 5 mismatches, as in the case above, you'll end up with an empty DataFrame. By using the drop_duplicates technique described above you avoid this problem, because in the cases where you don't have 5 mismatches, you just get the next-highest mm level (which is usually what you want)

A note for programmers
+++++++++++++++++++++++

If you'd like to edit inStrain to add functionality for your data, don't hesitate to reach out to the authors of this program for help. Additionally, please consider submitting a pull request on GitHub so that others can use your changes as well.
