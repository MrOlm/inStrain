# inStrain

[![Downloads](https://pepy.tech/badge/instrain)](https://pepy.tech/project/instrain)
[![Downloads](https://pepy.tech/badge/instrain/week)](https://pepy.tech/project/instrain)

inStrain is a python program for analysis of co-occurring genome populations from metagenomes that allows highly accurate genome comparisons, analysis of coverage, microdiversity, and linkage, and sensitive SNP detection with gene localization and synonymous non-synonymous identification.

Manual, installation instructions, and expected output are at available at [ReadTheDocs](https://instrain.readthedocs.io/en/latest/)

Publication is available in [Nature Biotechnology](https://doi.org/10.1038/s41587-020-00797-0) and on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.22.915579v1)

## Installation options

### pip
```
$ pip install instrain
```

### bioconda
```
$ conda install -c conda-forge -c bioconda -c defaults instrain
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
