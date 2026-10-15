#!/bin/bash
set -euo pipefail

# establish parameters?

# read from a text file all of the sample numbers to be run

############## Step 0: Gather Reference Genome if Needed ##############
# make sure that the tools for this script are enabled
#conda init
#conda deactivate
#conda activate 6701-rnaseq

# read all samples to be processed
#SRA_LIST=$1

# input index variable with param

############## Step 1: Perform Index Build ##############

# Fast works by not getting known splice sites and exon sites for indexing

# build index, but first check if the indexes are already present
build=0
for i in {1..8}; do
    if [ ! -f "${REF_DIR}/potato_dm_v6.1.${i}.ht2" ]; then
        build=1
    fi
done

if [ $build == 1 ]; then
    # build the index if the above condition is met
    echo
    echo "performing hisat2 build"
    echo
    hisat2-build -p 4 "${REF_DIR}/potato_dm.v6.1.fa" "${REF_DIR}/potato_dm_v6.1"

else
    echo "index has already been built for potato genome"

fi