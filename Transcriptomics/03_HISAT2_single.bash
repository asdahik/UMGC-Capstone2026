#!/bin/bash
set -eou pipefail

PROJ_DIR=$(pwd)

# this isnt SRA but PRJNA
PRJNA=$1
# set directory

RESULTS="${PROJ_DIR}/Transcriptomics/results/${PRJNA}"

# the fastp step probably covered this, but add the conditional again
if [ ! -d "${RESULTS}/bam" ]; then
    mkdir "${RESULTS}/bam"
fi

TRIMMED="${RESULTS}/trimmed"
echo "$TRIMMED"
zcat "${TRIMMED}/SRR3161990_trimmed.fastq.gz" | head -n 4 || true
IDX="${PROJ_DIR}/Transcriptomics/ref/potato_dm_v6.1"
#echo "$IDX"
RAW_DATA="${PROJ_DIR}/Transcriptomics/data/${PRJNA}_tmp"

# gather raw_data read

while read line; do 
    echo "$line"
    r1="${TRIMMED}/${line}_trimmed.fastq.gz"
    #echo "r1=$r1"
    ls -lh "$r1"
    bam="${RESULTS}/bam/${line}.sorted.bam"
    hisat2 -p 8 -x $IDX -U $r1 \
    | samtools sort -@ 4 -o "$bam" -

    samtools index "$bam"


done < "${PROJ_DIR}/Transcriptomics/${PRJNA}.txt"