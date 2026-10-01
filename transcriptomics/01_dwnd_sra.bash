#!/bin/bash
set -euo pipefail

# the main directory is assumed to be UMGC-Capstone2026
PROJ_DIR=$(pwd)
# this is the script that requests SRA samples from NCBI. It is run nested insided of 01_dwnd_sra_samples.bash
SRA=$1
#TODO: specify an output -O for SRA

# fasterq-dump comes with built ins to check that prefetched data has been downloaded already :)
RAW_DATA="${PROJ_DIR}/transcriptomics/data"
if [ ! -f "${SRA}.fastq" ]; then
    # prefetch condition, perhaps we already prefetched
    prefetch "${SRA}" -O "${RAW_DATA}"
    echo "prefetch complete for ${SRA}, moving onto fasterq-dump"
    fasterq-dump "${SRA}" -O "${RAW_DATA}/" --split-files --threads 2 # this download might be slow
fi