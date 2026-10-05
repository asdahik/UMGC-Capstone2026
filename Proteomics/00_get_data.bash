#!/bin/bash
set -eou pipefail

# the first variable that is past through along with the bash script is the PRIDE project number

# entry format using a system input variable by the user to call this script
#bash Proteomics/00_get_data.bash <PXD_Accession_Number>
# example used in data - PXD020259
PXD=$1
PROJ_DIR=$(pwd)

# create data and result directories
PROT_DATA=$PROJ_DIR/Proteomics/data
PROT_RESULT=$PROJ_DIR/Proteomics/results
#wget -r -np -nd -c -A "*.mgf,*.mgf.gz,*mgf.zip" \ 
#    -P ${PROJ_DIR}/data \
#    https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD020259/ 

# this will grab all mass spectrometry data from the project
# using pridepy to extract the files from PRIDE repo

#### Step 0: Generate a JSON of metadata from download link ####

if [ ! -f "${PROT_DATA}/${PXD}.json" ]; then

    echo "${PXD}.json not detected. Downloading metadata now from PRIDE using pridepy"
    pridepy stream-files-metadata -a "${PROT_DATA}/${PXD}" -o "${PROT_DATA}/${PXD}.json"

else

    echo "${PXD}.json is already downloaded. Proceeding to collecting MGF and MZID file names."

fi

# now generate txt files of collected MGF and MZID file names from the json

if [ ! -f "${PROT_DATA}/${PXD}_MGF.txt" ]; then

    # finding all files with .mgf
    grep -oE -e "[A-Z][0-9]+.mgf" \
        -e "[A-Z][0-9]+.mgf.gz" \
        -e "[A-Z][0-9]+.mgf.zip" "${PROT_DATA}/${PXD}.json" | sort -u > "${PROT_DATA}/${PXD}_MGF.txt" 
fi

if [ ! -f "${PROT_DATA}/${PXD}_MZID.txt" ]; then
    grep -oE -e "[0-9]+.mzid.gz" "${PROT_DATA}/${PXD}.json" | sort -u > "${PROT_DATA}/${PXD}_MZID.txt"
fi

# adding a time sleep function for readability
time sleep 5

#### Step 1: Download MGF and MZID data if not completed ####

# now that the PXD file is created, download the samples
if [ ! -d "${PROT_DATA}/${PXD}" ]; then

    echo "PXD directory not detected, extracting MGF and MZID data from PRIDE"

    # MGF filr collection step
    while read line; do
        pridepy download-file-by-name -a ${PXD} -f ${line} -o "${PROT_DATA}/${PXD}"
    done < "${PROT_DATA}/${PXD}_MGF.txt"
    
    # MZID file collection step
    while read line; do
        pridepy download-file-by-name -a ${PXD} -f ${line} -o "${PROT_DATA}/${PXD}"
    done < "${PROT_DATA}/${PXD}_MZID.txt"

    time sleep 5
fi

echo "All .mgf and .mzid files downloaded from PXD directory."

#### Step 2: Run Casanovo, generate a mztab file for PXD ####

# Run Casanovo? 
if [ ! -f "${PROT_RESULT}/${PXD}.mztab" ]; then
    # clean all of the downloaded files to fit the format by removing the header
    for f in "${PROT_DATA}/${PXD}"/*.mgf; do
        mkdir -p ${PROT_DATA}/tmp
        sed -n '/^BEGIN IONS/,$p' "$f" | sed 's/\r$//' > "${PROT_DATA}/tmp/$(basename "$f")"
    done
    echo "Proceeding to sequencing with Casanovo."
    casanovo sequence "${PROT_DATA}"/tmp/*.mgf --output_dir "${PROT_RESULT}" --output_root ${PXD}

    # remove tmp directory
    rm -rf ${PROT_DATA}/tmp/*
    rmdir ${PROT_DATA}/tmp
fi

#### Step 3: Generate a tsv table mapping each MZID file to the associated MGF file ####

for f in "${PROT_DATA}/${PXD}/*.mzid" "${PROT_DATA}/${PXD}/*.mzid.gz"; do
    echo $f
    [ -e "$f" ] || continue
    loc=$(zcat -f "$f" | grep -m1 -o '<SpectraData[^>]*location="[^"]*"' | sed 's/.*location="//; s/"$//')
    n=$(zcat -f "$f" | grep -c '<SpectraData ')
    echo -e "$f\t$(basename "$loc")\t$n"
done | tee "${PROT_RESULT}/${PXD}_mzid_to_mgf.tsv"
