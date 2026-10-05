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

RAW_DATA="${PROJ_DIR}/Transcriptomics/data"

# gather raw_data read

while read line; do 

    if [ ! -f "${TRIMMED}/${line}_trimmed_R2.fastq.gz" ]; then
        echo "Performing trimming and filtering of ${line}..."
        # line is read from the text file
        fastp \
            -i "${RAW_DATA}/${line}_R1.fastq.gz" -I "${RAW_DATA}/${line}_R2.fastq.gz" \
            -o "${TRIMMED}/${line}_R1trimmed.fastq.gz" -O "${TRIMMED}/${line}_R2trimmed.fastq.gz" \
            --qualified_quality_phred 20 \
            --length_required 100 \
            --detect_adapter_for_pe \
            --thread 4 \
            --json "${TRIMMED}/${line}_trimmed.json"

        # quick comparison to visualize while in the command line
        echo "Before: $(( $(wc -l < "${RAW_DATA}/${line}_R1.fastq.gz") / 4 )) read pairs"
        echo "After: $(( $(wc -l < "${TRIMMED}/${line}_R1trimmed.fastq.gz") / 4 )) read pairs"
    fi

done < "${PROJ_DIR}/Transcriptomics/${PRJNA}.txt"