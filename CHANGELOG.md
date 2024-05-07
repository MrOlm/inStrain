# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [1.9.0] - 2024-05-07
- Updates from Jing Wang pull request #181 https://github.com/MrOlm/inStrain/pull/181
- Efficiency improvements for inStrain compare and Error message clarification

## [1.8.1] - 2024-04-12
- Minor update to _get_covt_keys method

## [1.8.0] - 2023-08-24
- If you don't have any genomes detected, crash more gracefully
- Make it so all the genome_info columns are the same whether or not you have a linkage table
- Add "maximum_reads" argument

## [1.7.6] - 2023-07-26
- Make polymorpher work when scaffold names are ints
- Fix plots 2 and 7 when there's no linkage reported

## [1.7.5] - 2023-05-16
- Silly fix of 1.7.4 again

## [1.7.4] - 2023-05-16
- Fix a division by 0 when you don't have any annos detected in a sample

## [1.7.3] - 2023-05-11
- Fix but in 1.7.2 that resulted in all 0s with parse_gene_annotations

## [1.7.2] - 2023-05-09
- Update to parse_gene_annotations that spells things out on the genome level

## [1.7.1] - 2023-02-24
- Bring things up-to-date and working on python 3.10
- Make profile_genes use "fork" instead of "spawn" to make it work with higher python versions (even though it's already deprecated)
- Mark specific tests to skip when samtools isn't installed
- Remove "pd.append" calls (to avoid lots of warnings)
- Make all tests work with python3.10
- Remove numba and scikit-learn as dependencies (make numba still work if it's installed)
- Lots of refactoring to get around pandas warnings

## [1.7.0] - 2023-02-20
- Update help; remove deprecated things
- Added "parse_gene_annotations" method
- Added auxiliary script "recluster_instrain_compare"

## [1.6.4] - 2023-01-06
- Update numpy requirement to avoid numba bug

## [1.6.3] - 2022-08-27
- Fix SNP pooling bug https://github.com/MrOlm/inStrain/pull/111

## [1.6.2] - 2022-06-27
- Add "dtype='float64'" to all pd.Series() calls to supress warnings

## [1.6.1] - 2022-06-16
- Change a pandas "append" statement that was making a lot of noise
- Fix a bug when scaffold names are numbere

## [1.6.0] - 2022-06-10
- Include SNV pooling options (big update with lots of changes)
- Include "cov1 = pd.Series(dtype='float64')" in readComparer.py to suppress a future warning
- Remove cryptic SNVs from the nonredundant_snv_table (https://github.com/MrOlm/inStrain/issues/102)
- Only produce the run report when run in debug mode
- Various minor stylistic tweaks in runtime reporting

## [1.5.7] - 2022-02-09
- Make error log more readable - https://github.com/MrOlm/inStrain/issues/79
- Correct the reporting for the amount of reads removed during filtering to account for different filtering schemes
- Adjust the docs to account fo adjusting the pairing filter (https://github.com/MrOlm/inStrain/issues/70)

## [1.5.6] - 2022-02-09
- Logging handling of broken genes
- Check for sambamba in rarefaction_curve.py

## [1.5.5] - 2021-10-11
- Actual traceback in the log of gene failure exception

## [1.5.4] - 2021-06-24
- Improve the speed of auxillary script "rarefaction_curve.py"
- Fix a bug in the Docker "run_instrain.py" script related to using compressed .fasta files

## [1.5.3] - 2021-03-26
- BugFix related to inStrain clustering (https://github.com/MrOlm/inStrain/issues/52)
- BugFix related to genes problems (https://github.com/MrOlm/inStrain/issues/53)
- BugFix related to logging problems; caused by having ";" in gene names (https://github.com/MrOlm/inStrain/issues/49) (I actually didn't fix this)
- Add "force_compress" option (https://github.com/MrOlm/inStrain/issues/54)

## [1.5.2] - 2021-03-02
- BugFix when calling mutation effects (N, S, etc.) that impacted mutations with no reads supporting anything but the consensus variant
- Added a note about the weird way that positions are calculated in the "mutation" column of SNVs.tsv

## [1.5.1] - 2021-02-18
- Minor updates to docs; include warnings about 0-based indexing and some more FAQs
- Sam to bam conversion now uses the correct number of processes
- Automatically compress big dataframes on storage
- Do gene profiling of cryptic SNPs as well
- Test whether genome distance matricies contain blanks before clustering

## [1.5.0] - 2021-02-11
- Updated internals for plotting
- Removed dRep as a dependency
- Include depenency check option and automatically log dependencies to log
- Fixed numerous plotting bugs
- Supress pandas not report chainedassignment warnings

## [1.4.1] - 2021-01-25
- Add the auxillary script "rarefaction_curve.py"
- Update the Docker to work again

## [1.4.0] - 2020-10-30
- Big refactor of the code for "inStrain compare"
- Add support for .stb files in inStrain compare
- Add a "database mode" for inStrain compare
- Add the option to give inStrain compare a single genome to compare
- Completely refactor testing to work with pytest
- "compare" now clusters genomes when an .stb is provided
- gene_info.tsv no longer reports genes with coverage == 0
- Refactor testing of Docker images as well
- Adjust the 

## [1.3.11] - 2020-10-21
- Change the dependencies to specify a bioconda version that has the alphabet thing (https://github.com/MrOlm/inStrain/issues/27)

## [1.3.10] - 2020-10-21
- Change the dependencies to not allow pandas version 1.1.3 (which causes the following bug: https://github.com/pandas-dev/pandas/issues/37094)

## [1.3.9] - 2020-09-29
- More multiprocessing patching for compare

## [1.3.8] - 2020-09-23
- More multiprocessing patching

## [1.3.7] - 2020-09-15
- More multiprocessing patching

## [1.3.6] - 2020-09-14
- Apply a mutliprocessing patch to help python version 3.6
- Allow the Docker to accept a FOF for inStrain compare input
- Have the Docker download a reduced IS profile for running compare**

## [1.3.5] - 2020-09-10
- Backwards compatibility when running inStrain compare on old inStrain profiles
- Print some more info before crashing due to no scaffolds detected
- Update the Docker to be able to handle compare

## [1.3.4] - 2020-08-31
- More updates to the docs
- Bugfix regarding .stb files when scaffolds are present in the .fasta file, but not the .bam file

## [1.3.3] - 2020-08-31
- More upgrades to the docs
- An additional refinement to creating profile_merge jobs in linear time even when the amount of genes loaded is huge

## [1.3.2] - 2020-08-28
- Fixed a bug in quick_profile where it didn't work without an .stb file
- Optimize gene initilization
- Put the split merging + gene profiling in groups
- Add logging of how long groups take to run
- Edited the Docker image to work with version 1.3
- Overhaul of the Glossary
- Update the Docker

## [1.3.1] - 2020-08-19
- Undid some numba that actually slowed things down
- Change the way iterate_commands works
- Make compare log it's multiprocessing efficiency
- Save all counts into counts_table.npz (thanks https://github.com/apcamargo)
- Avoid sorting pre-sorted BAMs and use multiple threads to index and sort BAMs (thanks https://github.com/apcamargo)
- Handle "--version" in argparse correctly
- Make profile properly handle profile_genes
- Add a "DEPRECATED" flag to standalone profile_genes module

## [1.3.0w] - 2020-08-13
- Significant refactoring of controller.py and profile
- Re-writting the test suite to be in multiple modules
- Add numba to filter_reads evaluate_pair method (and a few others; just playing around)
- Optimize compare a bit (load stuff up front)

## [1.3.0v] - 2020-08-10
- Change internal structure of test suite
- Delete N_sites and S_sites from gene_info table
- Add "class" to SNVs.tsv
- Add some basic checkpoint logging to Compare
- Add a little bit of documentation to log_utils
- Make "compare" multi-thread in the new way (with spawn)
- Tiny docs change

## [1.3.0u] - 2020-08-06
- Add Docker and conda installation instructions to the README
- Edit parse_stb to handle None the same as []
- Add the Docker image and associated files

## [1.3.0t / 1.3.0.dev3] - 2020-07-29
- Fix "call_con_snps" to account for cases where there's only a SNP in one sample,
  but it's not a consensus SNP
- Make the output of compare generated through SNVprofile (like profile does it)
- Make the SNP table produced by compare actually legable and made in the output

## [1.3.0s] - 2020-07-14
- Add "FailureScaffoldHeaderTesting" to readComparer to make sure it can catch exceptions
- Add "high_cov" testing to test over 10,000x coverage
- Fix bug in readComparer when coverage was over 10,000x

## [1.3.0r / v1.3.0.dev1 / v1.3.0.dev2] - 2020-06-17
- UPLOADED AS v1.3.0.dev1 TO PYPI
- Fix a bug that broke iRep when skipping mm
- Translate .fasta files into all uppercase on loading them
- Removed some of the uninteresting columns from the genome_info output
- Change name of argument "filter_cutoff" to "min_read_ani"
- Add database_mode

## [1.3.0q] - 2020-06-15
- Fix a bug in GW having to do with the new mapping_info

## [1.3.0p] - 2020-06-11
- Add a little bit more checkpoints
- Fix the "min_genome_coverage" thing; there were problems with the read report when a scaffold had 0 reads

## [1.3.0o] - 2020-06-10
- Modify how split profiling is done and modify Rdic. Now only the reads for each scaffold are passed to worker threads, through the queue. This should lead to increased run-times (though hopefully not too bad; it should only impact efficiency) and significantly decreased memory usage
- Remove "--scaffold_level_mapping_info" and have it always on
- Change how read filtering is done after mutli-processing to improve speed

## [1.3.0n] - 2020-06-08
- Deep copy when making the final Rdic object in filter_reads

## [1.3.0m] - 2020-06-05
- Add logging for spawning and terminating workers to check the impact on RAM
- No more global scaff2sequence; just send over with commands
- Explicitly call the garbage collector in profile

## [1.3.0l] - 2020-06-04
- Add more checkpoints to filter_reads
- Adjust "ReadGroupSize" to 5000
- iRep is only run on mm = 1
- re-write \_iRep_filter_windows
- change the null_model to a regular dictionary; if a coverage isn't in there, use null_model[-1] for the biggest coverage
- change profile to no longer use globals; spawn new processes instead of forking

## [1.3.0k] - 2020-06-04
- Add a sanity check to "prepare_bam_fie" for proper indexing
- Dont have get_paired_reads write to the log when multiprocessing
- Remove \_validate_splits (it takes too long)
- Add a checkpoint for loading the .fasta file
- Change filter_reads to the "spawn" rather than "fork" multiprocessing
- Change the way bam files are initilized in filter_reads

## [1.3.0j] - 2020-06-02
- Make a argparse group for mm_level; calc GW on the mm level
- Plots work again
- Delete a pbar update that was crashing gene profiling in single thread mode (I think?)
- Fix partial gene identification problem
- If .sorted.bam in the name of the .bam file, dont try and sort it
- There can only be a maximum of 1000 jobs when filtering reads; hopefully cuts way down on queue access

## [1.3.0i] - 2020-06-01
- Fix bug with small scaffolds and add test to catch it
- Change --min_fasta_reads to --min_scaffold_reads
- Add --min_genome_coverage to profile, along with tests to make sure it works
- Handle iRepErro and StbErrors now

## [1.3.0h] - 2020-05-22
- Add some logging to gene profile in single thread mode
- Add the ability to catch iRep failures

## [1.3.0g] - 2020-05-21
- Add a 5 second timeout for single thread runs
- Store fasta_loc during profile
- Store scaffold 2 length from profile, not from within the profiled splits
- iRep can now be calculated
- Replace _mm_counts_to_counts_shrunk_ with a version that preserves base order
- Add "object_type" to IS profile objects
- If interrupted during profile, kill all processes
- Output tables are all made using the "generate" command of SNVprofile
- Big changes to the actual output tables
- This version breaks lots of plots and the "compare" function! Watch out!

## [1.3.0f] - 2020-05-13
- dN/dS for genes

## [1.3.0e] - 2020-05-13
- Re-write of genes profiling section
- Suppress warnings made during plotting
- Add real logging report to profile_genes
- Add failure testing to profile_genes

## [1.3.0d] - 2020-05-08
- Updating failure report to be better (including testing for split profiling and merging failures)
- Updating plots to not make all those errors
- Removed "MM" from all figures, replace with ANI level
- Rdic is stored by default

## [1.3.0a-c] - 2020-05-08
- Paralellization of individual scaffolds is done now using windows
- This required various performance tweaks to get right

## [1.2.14] - 2020-05-04
- RAM log actually reports min RAM usage now
- Dont store mm_reads_to_snvs by default
- Checkpoints report RAM usage as well
- Edits to make quick_profile have a similar input structure as profile

## [1.2.13] - 2020-04-20
- Correct RAM log
- Make a ram profiling plot

## [1.2.12] - 2020-04-09
- Bug fixes with logging
- Bug a problem with refBase being N calling multi-allelic SNPs

## [1.2.11] - 2020-04-08
- Changes to speed up multiprocessing
- gene2sequence is now a global
- removes sending over the IS object and Gdb
- runs 6000 genes per parallelization minimium
- speeds up creating figures by cacheing
- log now reports year and seconds
- generates a report at the end with runtime notes

## [1.2.10] - 2020-04-03
- Change the way that profile_genes works on the backend to optimize speed
- Turn off complete scaffold-level read profiling by default; make it a command line option
- Other small speed improvements

## [1.2.9] - 2020-03-29
- Dont crash if you have no gene clonality or coverage

## [1.2.8] - 2020-03-23
- Fix loading genes from GenBank files
- Allow storing gene2sequence

## [1.2.7] - 2020-03-21
- Calculates rarefied microdiversity as well, now

## [1.2.6] - 2020-03-19
- Don't crash geneprofile when no SNPs are called

## [1.2.5] - 2020-03-17
- Allow genbank files for calling genes (FILE MUST END IN .gb OR .gbk)

## [1.2.4] - 2020-02-27
- No longer crash when "N" is in the reference sequence
- Add the ability to --use_full_fasta_header

## [1.2.3] - 2020-02-21
- Messed with the internals of filter_reads. Moved old methods into deprecated_filter_reads
- Changed the options of filter_reads to just do things that it can actually do at the moment
- Remove break in "get_paired_reads" when a scaffold isn't in the .bam
- Fixed bug in profile_genes resulting from having no SNPs called
- Allow changing the options for pairing filters
- Allow creation of detailed read report
- Allow specification of priority_reads that will not be filtered out

## [1.2.2] - 2020-01-23
- Have genome_wide report microdiversity if it can
- Tried to fix a bug where the "percent_compared" in "inStrain compare" is underreported when run in "genome_wide"

## [1.2.1] - 2020-01-21
- Fixed typo in spelling of "Reference_SNPs"
- Fixed bug reporting refFreq
- Add a lot of "fillna(0)"s to the method "genome_wide_si_2" to make the average coverage work out correctly

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
