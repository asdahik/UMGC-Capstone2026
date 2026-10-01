#!/bin/bash
set -euo pipefail

# the goal of this script is to extract all SRA samples for a project accenssion number (PRJNA###) given as a variable input
# using NCBI's entrez-direct and produces a list of SRA samples from that project accension number. 
# The SRA samples are then read by a for loop that uses NCBI's sratoolkit's fasterq-dump to then download each 
# SRA sample for data processing. 

##### NOTE: This can become very memory intensive and should only be executed per project. Raw data extracted from NCBI that is
##### no longer needed for data processing should be deleted.


PRJNA=$1
PROJ_DIR=$(pwd)

if [ ! -f "SRR_${PRJNA}.txt" ]; then

    # attempt to query the inputted variable above. Esearch may return some curl command fail error logs but that is alright
    esearch -db sra -query "${PRJNA}" | efetch -format runinfo | cut -d ',' -f 1 | grep SRR > "SRR_${PRJNA}.txt"

else
    echo "SRR_${PRJNA}.txt is already queried and extracted from Entrez"
fi

# using xargs to grab two samples at a time, can be modified to fit specifications and request access frequency permissions with NCBI

#TODO: adjust project directories for scripts
RAW_DATA_SCRIPT="${PROJ_DIR}/transcriptomics/01_dwnd_sra.bash"
chmod -x "${RAW_DATA_SCRIPT}"
# this runs the get raw data script
xargs -a "SRR_${PRJNA}.txt" -P 2 -I{} bash "$RAW_DATA_SCRIPT" {}