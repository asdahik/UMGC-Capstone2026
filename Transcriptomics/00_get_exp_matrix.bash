#!/bin/bash
PROJ_DIR=$(pwd)/Transcriptomics
PRJNA=$1
RESULTS=$PROJ_DIR/results

# Need to also download SRA Run metadata

cd $RESULTS
# additionally, before GSEA, download the SRA table from run selector
# this data is stored in results instead of ref like the SpudDB data
esearch -db sra -query "${PRJNA}" | efetch -format txt > SraRunTable.txt