# README for the `instrain` Docker Image

## Overview

A Docker image containing inStrain is available on [Docker Hub with name "mattolm/instrain:latest"](https://hub.docker.com/repository/docker/mattolm/instrain).

The Docker image was made by running the command `make` in this directory. This puts all files in this directory inside of the resulting Docker image, as well as the most recent version of `inStrain` available on [Bioconda](https://bioconda.github.io/recipes/instrain/README.html). If you want to make sure inStrain is up to date, or use a different version, you can run `pip install instrain --upgrade` within the Docker image.

**Before inStrain can be run in the Docker image it must be initialized with the command ./prepare.sh; conda activate work;**

After inStrain is initialized it can be run as normal by moving data into the Docker image, or by using the `run_instrain.py` helper script to download data from AWS S3.

## `run_instrain.py`

This helper script was designed following the example provided by [aws-batch-genomics](https://github.com/aws-samples/aws-batch-genomics/tree/v1.0.0/tools), and is most useful for use in [AWS Genomic Batch Workflows](https://aws.amazon.com/blogs/compute/building-high-throughput-genomic-batch-workflows-on-aws-batch-layer-part-3-of-4/). Below are the arguments that can be provided to this wrapper script, as well as an example command of how it can be run.

```
usage: run_instrain.py [-h] [--bam BAM] [--fasta FASTA] [--genes GENES]
                       [--stb STB] [--results_directory RESULTS_DIRECTORY]
                       [--wd_name WD_NAME] [--command COMMAND]
                       [--cmd_args CMD_ARGS] [--light_upload]
                       [--working_dir WORKING_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --working_dir WORKING_DIR

File paths:
  --bam BAM             s3 path to the .bam file
  --fasta FASTA         s3 path to the .fasta file
  --genes GENES         s3 path to the genes file
  --stb STB             s3 path to the stb file
  --results_directory RESULTS_DIRECTORY
                        s3 path to the folder to put the results in

Run command args:
  --wd_name WD_NAME     Name of the output directory
  --command COMMAND     The command that should go after the command inStrain
  --cmd_args CMD_ARGS   A string (as long as you want) that will be put after
                        the inStrain command, .bam, fasta, and output
                        directory
  --light_upload        By default it will upload the /raw_data folder; this
                        will make it not
```

### Example commands

```
./prepare.sh; conda activate work; ./run_instrain.py --bam {0} --fasta {1} --results_directory {2} --wd_name {3} --cmd_args='--skip_plot_generation'
```

```
"./prepare.sh; conda activate work; ./run_instrain.py --IS {0} {1} --results_directory {2} --wd_name {3} --command compare --cmd_args='--store_mismatch_locations' --scaffolds {4}
```