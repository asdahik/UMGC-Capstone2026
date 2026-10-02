#!/bin/bash
set -euo pipefail

# the main directory is assumed to be UMGC-Capstone2026
PROJ_DIR=$(pwd)
# this is the script that requests SRA samples from NCBI. It is run nested insided of 01_dwnd_sra_samples.bash
SRA=$1
#TODO: specify an output -O for SRA

# fasterq-dump comes with built ins to check that prefetched data has been downloaded already :)
RAW_DATA="${PROJ_DIR}/transcriptomics/data"
if [ ! -f "${RAW_DATA}/${SRA}.fastq.gz" ]; then

    # prefetch will check if the prefetched SRA file is present.
    # it does its own conditional check when searching for it so
    # no need to establish a conditional check like the one present
    prefetch "${SRA}" -O "${RAW_DATA}"
    echo "prefetch complete for ${SRA}, moving onto fasterq-dump"
    
    # Extract then compress
    if [ ! ]
    fasterq-dump "${SRA}" -O "${RAW_DATA}/" --split-files --threads 4 # this download might be slow
    #gzip the extracted file for now to save space, its close to the prefetch size
    gzip -v "${RAW_DATA}/${SRA}.fastq"
    # once the fasterq-dump process has completed, remove the uncompressed file
    rm -rf "${RAW_DATA}/${SRA}.fastq"
    # save the prefetched file just in case. Like sending it to a database potentially?
    #rm -rf "${RAW_DATA}/${SRA}"
    
    rm -rf fasterq.tmp.*

else
    echo "${SRA}.fastq is already extracted"
fi