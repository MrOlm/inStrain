#!/bin/bash

echo "shebang present"

# Set up
BINS="s3://czbiohub-microbiome/Sonnenburg_Lab/hwastyk/201912-fiber-fermented-study/Pilot_8024_6_4_B6/08_BINNING/METABAT215_SUBJECTMAPPING_SCAFFOLDS_1500_MERGED_K77_METASPADES__Pilot_8024_6_4_B6.binset.txt"
FASTA="s3://czbiohub-microbiome/Sonnenburg_Lab/hwastyk/201912-fiber-fermented-study/Pilot_8024_6_4_B6/06_ASSEMBLY/SCAFFOLDS_1500_MERGED_K77_METASPADES__Pilot_8024_6_4_B6.fasta"

BINS_NAME=$( basename $BINS )
FASTA_NAME=$( basename $FASTA )
OUT_LOCATION="s3://czbiohub-microbiome/Sonnenburg_Lab/Software/docker_testing/s3_results/"
OUT_BASE="testeroini"

# Make dirs
mkdir input_data; mkdir output_data

echo "which python"
which python

echo "which pip"
which pip

echo "which conda"
which conda

# Download input_data
aws s3 cp $BINS input_data/
aws s3 cp $FASTA input_data/

# Run the process
python BioScripts/parse_stb.py -s input_data/$BINS_NAME -f input_data/$FASTA_NAME -o output_data/$OUT_BASE
gzip output_data/*

# Upload the results
aws s3 sync output_data/ $OUT_LOCATION
