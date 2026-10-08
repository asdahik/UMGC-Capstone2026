#!/bin/bash

set -eou pipefail

PRJNA=$1

PROJ_DIR=$(pwd)/Transcriptomics
ANNOTATION=${PROJ_DIR}/ref/potato_genome_annotation.v6.1.gtf
RESULT=${PROJ_DIR}/results/${PRJNA}
featureCounts -T 8 \
    -s 0 \
    -t exon -g gene_id \
    -a $ANNOTATION \
    -o $RESULT/${PRJNA}_counts.txt \
    $RESULT/bam/*.bam