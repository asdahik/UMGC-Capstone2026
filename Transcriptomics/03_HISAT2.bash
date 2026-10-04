#!/bin/bash
set -eou pipefail

PROJ_DIR=$(pwd)
RESULTS="${PROJ_DIR}/transcriptomics/results"

for r1 in ${RESULTS}/trimmed/*_1.fastq.gz; do
    r2=${r1/_1.trim/_2.trim}
    base=$(basename "$r1" _1.trim.fastq.gz)

    hisat2 -p 8 -x genome_index \
        -1 "$r1" -2 "$r2" \
        --rg-id "$base" --rg SM:"$base" \
        2> "logs/${base}.hisat2.log" \
        | samtools sort -@ 4 -o "${RESULTS}/bam/${base}.sorted.bam" -

    samtools index "${RESULTS}/bam/${base}.sorted.bam"
done