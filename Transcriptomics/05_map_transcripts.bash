#!/bin/bash

PRJNA=$1
PROJ_DIR=$(pwd)/Transcriptomics
REF=${PROJ_DIR}/ref
RESULTS=${PROJ_DIR}/results/${PRJNA}

# run DESeq2
echo "running DESeq2 script"
Rscript ${PROJ_DIR}/05_DESeq2.r ${RESULTS}/${PRJNA}_counts.txt ${RESULTS}/${PRJNA}/SraRunTable.csv
python3 ${PROJ_DIR}/05_map_transcripts.py $RESULTS $REF

# once the transcripts are mapped, create a map key .txt
cd ${RESULTS}
awk '{print $1 "\t" $7}' deseq_results_annotated.csv > transcript_map.txt