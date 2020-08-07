# inStrain

inStrain is python program for  analysis of co-occurring genome populations from metagenomes that allows highly accurate genome comparisons, analysis of coverage, microdiversity, and linkage, and sensitive SNP detection with gene localization and synonymous non-synonymous identification.

Manual, installation instructions, and expected output are at available at [ReadTheDocs](https://instrain.readthedocs.io/en/latest/)

Preprint publication is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.22.915579v1)

## Installation options

### pip
```
$ pip install instrain
```

### bioconda
```
$ conda config --add channels bioconda; conda install instrain
```

### Docker
Docker image is available on Docker Hub at [mattolm/instrain](https://hub.docker.com/repository/docker/mattolm/instrain). See [docker/](docker/) for use instructions.

## Quick start

### Show program help and modules:
```
$ inStrain -h
```

### Microdiversity and SNP-calling pipeline:
```
$ inStrain profile mapping.bam genome_file.fasta -o inStrain_profile1
```

### Detailed strain-level comparison:
```
$ inStrain compare inStrain_profile1 inStrain_profile2
```
