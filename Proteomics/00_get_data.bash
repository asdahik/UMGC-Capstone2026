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

if [ ! -f "${PROT_DATA}/${PXD}.txt" ]; then

    echo "${PXD}.txt not detected. Downloading now from PRIDE using pridepy"
    cd $PROT_DATA
    pridepy stream-files-metadata -a ${PXD} -o ${PXD}.json
    # finding all files with .mgf
    grep -oE -e "[A-Z][0-9]+.mgf" \
        -e "[A-Z][0-9]+.mgf.gz" \
        -e "[A-Z][0-9]+.mgf.zip" ${PXD}.json | sort -u > ${PXD}.txt 
    cd ../..
else

    echo "${PXD}.txt is already downloaded"

fi

# adding a time sleep function for readability
time sleep 5

# now that the PXD file is created, download the samples
if [ ! -d "${PROT_DATA}/${PXD}" ]; then

    echo "PXD directory not detected, extracting data from PRIDE"
    while read line; do
        pridepy download-file-by-name -a ${PXD} -f ${line} -o "${PROT_DATA}/${PXD}"
    done < "${PROT_DATA}/${PXD}.txt"
else
    echo "PXD directory detected. Redownloading any missing files from ${PXD}.txt"

    while read line; do
        if [ ! -f "${PROT_DATA}/${PXD}/${line}" ]; then
            echo "${line} from ${PXD} was detected missing, redownloading"
            pridepy download-file-by-name -a ${PXD} -f ${line} -o "${PROT_DATA}/${PXD}" 
        fi
    done < "${PROT_DATA}/${PXD}.txt"
    # add a pause
    time sleep 5
fi

echo "All .mgf files downloaded from PXD directory. Proceeding to sequencing with Casanovo."


# clean all of the downloaded files to fit the format
for f in "${PROT_DATA}/${PXD}"/*.mgf; do
    mkdir -p ${PROT_DATA}/tmp
    sed -n '/^BEGIN IONS/,$p' "$f" | sed 's/\r$//' > "${PROT_DATA}/tmp/$(basename "$f")"
done
# Run Casanovo? 
casanovo sequence "${PROT_DATA}"/tmp/*.mgf --output_dir "${PROT_RESULT}" --output_root ${PXD}

# remove tmp directory
rm -rf ${PROT_DATA}/tmp/*
rmdir ${PROT_DATA}/tmp