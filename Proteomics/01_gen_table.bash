#!/bin/bash

PXD=$1
PROJ_DIR=$(pwd)

# create data and result directories
PROT_DATA=$PROJ_DIR/Proteomics/data
PROT_RESULT=$PROJ_DIR/Proteomics/results

if [ ! -f "${PROT_RESULT}/${PXD}_match.tsv" ]; then
    for f in "${PROT_DATA}/${PXD}"/*.mzid "${PROT_DATA}/${PXD}"/*.mzid.gz; do

        [ -e "$f" ] || continue
        MZID=$(grep -oE "[0-9]+\.mzid(\.gz)?$" <<< "$f")
        #echo $MZID
        loc=$(zcat -f "$f" | grep -m1 -o '<SpectraData[^>]*location="[^"]*"' | sed 's/.*location="//; s/"$//')
        echo -e "$MZID\t$(basename "$loc")"
    # TSV generated is a Peptide-Spectra Match     
    done | tee "${PROT_RESULT}/${PXD}_match.tsv"
fi

# run python script to generate mascot values
python3 $PROJ_DIR/Proteomics/02_build_psm.py "${PROT_RESULT}/${PXD}_PSM.tsv" "${PROT_DATA}/${PXD}"/*.mzid.gz