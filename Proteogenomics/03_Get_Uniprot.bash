#!/bin/bash

# take results from PXD...._MS_Table.tsv after build_psm.py

PROJ_DIR=$(pwd)

PXD=$1
PROT_RESULTS=${PROJ_DIR}/$2/${PXD}
MS_TABLE=$PROT_RESULTS/${PXD}_MS_table.tsv
# run python script ... try to move output results into PROT results

cd ${PROT_RESULTS}

# generate #{PXD}_Prot_Ascession.txt

awk 'NR>1 {print $3}' ${MS_TABLE} > ${PXD}_Prot_Ascession.txt

# this python script generates a {pdx}_proteins.out fasta and {pdx}_annotation.tsv
python3 ../Uniprot_REST.py ${PXD}_Prot_Ascession.txt ${PROT_RESULTS}

# default to the annotation table now to create an additional table from the following ms_peptides script
ANNOT=$PROT_RESULTS/${PXD}_MS_table.tsv

# generate a fasta file with the casanovo peptides from $PXD_MS_TABLE
PEP_OUT="${PROT_RESULTS}/${PXD}_peptides.fasta"

awk -F'\t' 'NR>1 { gsub(/^[ \t]+|[ \t]+$/, "", $3)
                   print ">" $1 " " $3 "\n" $2 }' "$MS_TABLE" > "$PEP_OUT"

# run this script to generate ms_Result figures
python3 ../ms_peptides.py --proteins ${PXD}_proteins.out --peptides ${PXD}_peptides.fasta --ms-tables ${PXD}_MS_table.tsv --outdir ${PXD}

# 