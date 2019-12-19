# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [1.2.0] - 2019-12-18
- Move from 'ANI' in scaffold_profile to "conANI" and "popANI"
- Store information on lots more types of SNPs in scaffold_profile
- Remove snpsCounted
- Make it so profile can handle genes, genome_wide, and figure generation
- Make the output of profile a bit more pretty

## [1.1.3] - 2019-12-18
- Profile genes no longer drops allele_count

## [1.1.2] - 2019-12-04
- Change genome_wide to include sums of SNPs
- GeneProfile now reports reference SNPs as well in terms of N and S
- Many updates to plotting function (listed below)
- The debug option will now produce a stack track for failed plots
- Plot basename now included in figure creation
- You can now plot only genomes meeting a breadth requirement
- You can now specify which genomes to plot
- Fixed bug preventing plots 2 and 7 from being made in a lot of cases
- Plot 4 now plots only freq > 0.5 and a set 50 bins

## [1.1.1] - 2019-11-18
- Change genome_wide to account for ANI, popANI, and conANI

## [1.1.0] - 2019-11-18
- Big changes to inStrain compare (stores SNVprofiles as globals; calculates popANI and con_ANI always; big changes to how popANI is calculated)

## [1.0.2] - 2019-11-17
- Update to try and get the Null model properly installed

## [1.0.1] - 2019-11-16
- Plot 9 displays titles now
- Fixed some of the test data; was previously on a non-mm level
- Set a maximum figure size of 100 inches on the scaffold inspection plot
- Update the docs for example_output
- Edit the GitHub README

## [1.0.0] - 2019-11-08
- InStain only yells at you if the minor version is different
- Add microdiversity in addition to clonality for user-facing output
- Change the name "morphia" to "allele_count"
- Bugfix for geneProfile on lien 92 related to "if scaffold not in covTs"
- Move calculate_null.py into the helper_scripts
- Delete "R_scripts" and "notebooks" folders
- Delete combine_samples.py
- Change the name of combined_null1000000.txt to NullModel.txt
- Changes to internals of calculate_null.py; add parameters so that others could change it if they want to as well
- Change the default NullModel to go up to 10,000 coverage, with X bootstraps
- Add plots 8 and 9
- Small changes to ReadComparer to try and reduce RAM usage by decreasing the amount of SNPtable stored
- Plot 10 can now be made
- Add dRep to required list
- Change the null model to be mutlithreaded
- NullModel.txt has 1,000,000 bootstraps
- Fixed a bug in making plot2

## [0.8.11] - 2019-10-30
- Optimize GeneProfile; no need to load SNV table twice, and establish a global variable of SNV locations

## [0.8.10] - 2019-10-30
- The raw_snp_table is now stored correctly
- Added plots 3-5
- Added plots 6-7
- Fixed critical bug on GeneProfile coverage
- Made storage of counts table make sense; only when store everything

## [0.8.9] - 2019-10-26
- Added plot number 2

## [0.8.8] - 2019-10-26
- Numerous little changes:
- Requires pandas >=0.25
- Fixed runtime warnings in genomeWide
- When you're making readReport genome_wide, remove the "all_scaffolds" column. Otherwise the program thinks that one of the scaffolds isn't in the .stb file
- Plotting now uses fonttype 42 for editable text
- Compare by default now skips mismatch locations

## [0.8.7] - 2019-10-26
- Greedy clustering is implemented in an experimental way

## [0.8.6] - 2019-10-15
- Bugfix on GeneProfile

## [0.8.5] - 2019-10-15
- Update to genome_wide calculations on the mm-level

## [0.8.4] - 2019-10-12
- Add the "plot" module

## [0.8.3] - 2019-10-07
- Add the "genome_wide" module

## [0.8.2] - 2019-10-06
- Compare is now stored in an SNVprofile object as well
- The log of compare is stored in the log folder
- The log of profile is stored in the log folder

## [0.8.1] - 2019-10-05
- Add a README to the raw_data folder
- Add mean_microdiversity and median_microdiversity to the scaffold table in addition to clonality
- Scaffolds with no clonality are np.nan, not 0
- Change the name from ANI to popANI in ReadComparer
- GeneProfile now stores the SNP mutation types in the output
- The name of the IS is now stored in front of everything in the output folder
- GeneProfile now correctly identifies genes as being incomplete or not

## [0.8.0] - 2019-09-22
- Whole new way of running inStrain is born- everything under the same parser

## [0.7.2] - 2019-09-16
- Optimize GeneProfile (faster iteration of clonT and covT)

## [0.7.1] - 2019-09-16
- Lots more bug-fixes / speed increases to RC
- Make the default not to self-compare with RC
- RC can now read scaffold list from a .fasta file

## [0.7.0] - 2019-09-15
- Fix the generation of the null model (prior to this it had an off-by-one error) - THIS IS BIG AND IMPACTS THE RESULTS IN TERMS OF WHERE SNPs WILL BE CALLED
- Make GeneProfile save an output in the output folder
- Make RC look for SNPs in each others databases, not just compare consensus bases
- Allow explicit setting of the fdr for use in null model

## [0.6.7] - 2019-09-10
- Add that bandaid to inStrain in general- rerun scaffolds that failed to multiprocess

## [0.6.6] - 2019-09-09
- Make readComparer storage of coverage information optional (and default off)

## [0.6.5] - 2019-09-08
- RAM logging for ReadComparer

## [0.6.4] - 2019-09-06
- Add more plottingUtilities
- Fix load_scaff2pair2mm2SNPs
- Re-write GeneProfile.py to actually work
- Make readComparer also store coverage information

## [0.6.3] - 2019-08-27
- Change "percent_compared" to "percent_genome_compared"
- Add junk to genomeUtilities and plottingUtilities

## [0.6.2] - 2019-08-27
- Change the name "considered_bases" to "percent_compared" in readComparer

## [0.6.1] - 2019-08-27
- ReadComparer now works with v0.6 IS objects

## [0.6.0] - 2019-08-26
- Note: Im making this a new minor version more because I should have made the last one a minor version that that this deserves to be one
- Give SNVprofile the ability to only load some scaffolds
- Add readComparer ability to only select certain scaffolds
- Make readComparer do the better form of multiprocessing
- Actually a lot of readComparer changes

## [0.5.6] - 2019-08-22
- NOTE: THIS UPDATE DOES SEEM TO INCREASE RAM SPIKES (probably related to shrinking), BUT DECREASES RAM USAGE OVERALL
- Add increased logging to read filtering
- Remove non-zero values from base-wise stored things
- Fix getting clonalities from clonT dataframe (need to sort by mm first)
- Remove the need to get clonalities from clonT dataframe during initial profile
- Add some more attempted logging and saving of failed scaffolds

## [0.5.5] - 2019-08-20
- Add an awk filter to quickProfile

## [0.5.4] - 2019-08-20
- Bugfix with regards to how unmaskedBreadth and ANI were calculated in an mm-based way

## [0.5.3] - 2019-08-19
- Add the script "quickProfile"

## [0.5.2] - 2019-08-19
- Shrink base-wise vectors before saving
- Use compression on hdf5 files
- Optimize multiprocessing
- Add the ability to provide a list of .fasta files to consider

## [0.5.1] - 2019-08-15
- Improvements to logging of multiprocessing; done by default now

## [0.5.0] - 2019-08-14
- Make _get_paired_reads_ much more RAM efficient
- Store clonT and covT in hd5 format, which prevents that final RAM spike
- Store clonality as float32 instead of float64
- Debug now stores useful RAM information

## [0.4.11] - 2019-08-13
- Make _clonT_to_table_ much more RAM efficient
- Another bugFix for ReadComparer

## [0.4.10] - 2019-08-13
- Fixed a bug in ReadComparer related to translating conBase and refBase to ints

## [0.4.9] - 2019-08-12
- Added a try-catch around chunking in the main method
- Added length to the output of read comparer
- Changed some coverage reporting in readComparer - still not clear if its right
- Fixed a readComparer coverage bug

## [0.4.8] - 2019-08-08
### Changed
- Optimized ReadComparer

## [0.4.7] - 2019-06-21
### Added
- Added the `--skip_mm_profiling` command line option

## [0.4.6] - 2019-06-19
### Added
- Bugfix multiprocessing chunking; no longer
- readComparer is almost there, but not quite

## [0.4.5] - 2019-06-17
### Added
- Added chunking to the multiprocessing

## [0.4.4] - 2019-06-13
### Added
- Added the argument "min_scaffold_reads"

## [0.4.3] - 2019-06-13
### Added
- Debug option
- A bunch of un-usable code and tests related to readComparer. Need to continue work on it

### Changed
- pileup_counts is now a "store_everything" attribute and not stored by default

## [0.4.2] - 2019-06-13
### Changed
- Big changes to the SNV table
  - Logs "morphia"; the number of bases reaching the above criteria
  - Included bases with 1 morphia when its different from the reference base
  - Include the reference base at that position

## [0.4.1] - 2019-06-13
### Fixed
- gene_statistics works more natively with the new SNVprofile

## [0.4.0] - 2019-06-13
### Changed
- SNV profile is totally different now. Pretty sweet

## [0.3.2] - 2019-06-05
### Changed
- gene_statistics now doesn't look for partial genes
- gene_statistics now actually saves as the file name you want it to
- reorganized everything

### Added
- expected breadth on the scaffold level

## [0.3.1] - 2019-06-04
- `Gene_statistics.py` now working
- Fixed critical bug of allele_B = allele_b in the linkage table.

## [0.3.0] - 2019-05-15
- Now calculates D', and normalized versions of D' and R2
- min_snp is now
- Has the new gene_statistics.py and combine_samples.py scripts

## [0.2.8] - 2019-05-09
- fixed varbase = refbase error when there's a tie for most abundant variant
- made default output prefix be the fasta file prefix
- now generate scafold counts numpy array which show total variant counts for all positions by scaffold

## [0.2.7] - 2019-04-10
### Changed
- Change to how scaff2pair2info is handled

## [0.2.6] - 2019-04-09
### Changed
- Minor change to the pickle protocol and stuff

## [0.2.5] - 2019-04-07
### Changed
- Allow changing filtering criteria, and some basic testing to make sure it works
- Will convert .sams to .bams and index for you

## [0.2.4] - 2019-03-29
### Changed
- Add a little bit of logging and efficiency into read report

## [0.2.3] - 2019-03-26
### Changed
- multiprocessing of scaffolds re-factored to be better
- now produces readable dataframe outputs
- make a nice read report
- add a --read_report option for filter_reads

## [0.2.2] - 2019-03-26
### Changed
- paired read filtering is now multi-processed

## [0.2.1] - 2019-03-26
### Changed
- argparse now displays the version number and defaults
- runs Matts way by default

## [0.2.0] - 2019-03-26
### Fixed
- implemented logging module
- made read filtering into a saved table
- implemented and tested "filter_reads.get_paired_reads"
- multiprocessing of SNP calling implemented
- setup.py is made; can now install itself
- changed the quoting on CC's table outputs to default
- changed "level" to "filter_cutoff"
- changed "min_coverage" to "min_cov"
- account for read overlap when calculating read length (half code is there, but reverted)
- changed to "no_filter" stepper
- clonality is now called the same with CC and Matt's versions
- min_cov is called with >= now

## [0.1.2] - 2019-03-12
- Stepper is now all

## [0.1.1] - 2019-03-12
- Tests now work

## [0.0.4] - 2018-06-25
- Adding some of the graph stuff and functions. The Python
is really slow, in the future it should call node2vec C++ library.

## [0.0.2] - 2018-06-21
### Fixed
- Test suite now runs on it's own and makes sure it succeeds

## [0.0.2] - 2018-06-21
### Fixed
- Test suite produces a command
- Program doesn't crash (lotsa small bug fixes)

## [0.0.1] - 2018-06-21
### Added
- Changelog and versioning is born
- Test suite is born

### Fixed
- Dumb print error is fixed
- Import statement is pretty :)
