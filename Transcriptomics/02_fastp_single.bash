#!/bin/bash
set -eou pipefail

PROJ_DIR=$(pwd)

# this isnt SRA but PRJNA
PRJNA=$1
# set directory

RESULTS="${PROJ_DIR}/Transcriptomics/results"

if [ ! -d "${RESULTS}/${PRJNA}/trimmed" ]; then
    mkdir "${RESULTS}/${PRJNA}/trimmed"
fi
TRIMMED="${RESULTS}/${PRJNA}/trimmed"
PRJNA_DIR="${PROJ_DIR}/Transcriptomics/data/${PRJNA}_tmp"

# gather raw_data read


while read line; do 

    if [ ! -f "${TRIMMED}/${line}_trimmed.fastq.gz" ]; then
        echo "Performing trimming and filtering of ${line}..."
        # line is read from the text file
        fastp \
            -i "${PRJNA_DIR}/${line}.fastq.gz" \
            -o "${TRIMMED}/${line}_trimmed.fastq.gz" \
            --qualified_quality_phred 20 \
            --length_required 30 \
            --detect_adapter_for_pe \
            --thread 4 \
            --json "${TRIMMED}/${line}_trimmed.json"

        # quick comparison to visualize while in the command line
        echo "Before: $(( $(wc -l < "${PRJNA_DIR}/${line}.fastq.gz") / 4 )) read pairs"
        echo "After: $(( $(wc -l < "${TRIMMED}/${line}_trimmed.fastq.gz") / 4 )) read pairs"
    fi

done < "${PROJ_DIR}/Transcriptomics/${PRJNA}.txt"