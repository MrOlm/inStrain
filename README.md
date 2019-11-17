# inStrain

inStrain is python program for  analysis of co-occurring genome populations from metagenomes that allows highly accurate genome comparisons, analysis of coverage, microdiversity, and linkage, and sensitive SNP detection with gene localization and synonymous non-synonymous identification.

Manual, installation instructions, and expected output are at available at [ReadTheDocs](https://instrain.readthedocs.io/en/latest/)

## Installation with pip
```
$ pip install instrain
```

## Quick start

### Show program help and modules:
```
$ inStrain -h
```

### Microdiversity and SNP-calling pipeline:
```
$ inStrain profile mapping.bam genome_file.fasta -o inStrain_profile1
```

### Details strain-level comparison:
```
$ inStrain compare inStrain_profile1 inStrain_profile2
```
