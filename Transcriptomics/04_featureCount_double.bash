#!/bin/bash

set -eou pipefail

PRJNA=$1
#strandedness value dependent on the data. Look for between 60-80% read alignment
STRAND_VAL=$2

PROJ_DIR=$(pwd)/Transcriptomics
ANNOTATION=${PROJ_DIR}/ref/potato_genome_annotation.v6.1.gtf
RESULT=${PROJ_DIR}/results/${PRJNA}
featureCounts -T 8 \
    -p --countReadPairs -B -C \
     $STRAND_VAL \
    -t exon -g gene_id \
    -a $ANNOTATION \
    -o $RESULT/${PRJNA}_counts.txt \
    $RESULT/bam/*.bam