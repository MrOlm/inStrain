# strains_analysis
Code for analyzing population genomics in genome-resolved metagenomes

Requires: pysam, tqdm, BioPython.

### Usage

```python strainRep2.py -s 5 -f 0.05 -l 0.96 sorted_indexed_bam.bam scaffolds.fasta -o  output --log  log.txt
```

`python strainRep2.py -h` actually is pretty helpful, that's all of the documentation.

output: 3 tables (and a big python object). Linkage table (showing snp linkage), frequency table (showing SNPs and their frequencies), and clonality table (showing the clonality and coverage of each position - from this gene clonality can be calculated and compared to the genome average) (edited)
`-s 5` requires 5 reads to confirm a SNP, you can adjust depending on your coverage. `-f` means minimum snp frequency of 5%, `-l 0.96` means that read pairs must be 96% ID to reference. the statistics reported in the log file are also super useful

### Read Comparer

Here is a description of the values:

| Value  | Description | Calculation | Examples |
| ------------- | ------------- | -------------  | -------------  |
| percent_compared  | The percentage of bases in the scaffolds that are covered by both. | len([x for x in overlap if x])/len(overlap) | When ANI is np.nan, this must be 0. When comparing to self, this will be equal to unmasked breadth. If both scaffolds have 0 coverage, this will be 0.
| compared_bases_count | the number of considered bases | len([x for x in overlap if x]) | snps = compared_bases_count - (float(ani) * compared_bases_count) |
| coverage_overlap | The percentage of bases that are either covered or not covered in both | cov = len(coveredInBoth) / len(coveredInEither) | When comparing to self this is always 1 (unless there are no compared_bases, then its 0); If both scaffolds have 0 coverage, this will be 0 |


### inStrain

Here is a description of the values:

| Value  | Description | Calculation | Examples |
| ------------- | ------------- | -------------  | -------------  |
| unmaskedBreadth  | The percentage of bases hitting the minimum coverage requirement | len(clons) / lengt | This is really just the percentage of bases that have a clonality value |
