#!/bin/bash -x

# The next two lines are necessary because there are two different /mnt locations
# Without this odd copy step, Snakemake fails (other things do too).
LOCAL=$(pwd)
cp -pr * $LOCAL/
cd $LOCAL
export PATH="/opt/conda/bin:${PATH}"
