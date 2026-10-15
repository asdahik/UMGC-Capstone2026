#!/bin/bash
set -euo pipefail

# system variable for building index with splice sites and stuff
PROJ_DIR=$(pwd)/Transcriptomics # run within the directory

# create directory

if [ ! -d "${PROJ_DIR}/ref" ]; then
    mkdir "${PROJ_DIR}/ref"
fi
REF_DIR="${PROJ_DIR}/ref"
echo "${REF_DIR}"

if [ ! -d "${PROJ_DIR}/results" ]; then
    mkdir "${PROJ_DIR}/results"
fi

RESULTS="${PROJ_DIR}/results"


# Check if all data has been downloaded for this particular project set
# need to get reference genome. Using DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3 and DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

# check if genome assembly and annotation have been downloaded


# downlaod genome assembly
if [ ! -f "${REF_DIR}/potato_genome_assembly.v6.1.fa" ]; then

    # genome assembly
    wget https://spuddb.uga.edu/data/dm_v61/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

    # decompress genome
    gunzip DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

    # move the assembly to their respective folders
    mv DM_1-3_516_R44_potato_genome_assembly.v6.1.fa "${REF_DIR}/potato_genome_assembly.v6.1.fa"

else
    # TODO MAKE A MORE INFORMATIVE PRINT STATEMENT
    echo "potato genome assembly has already been downloaded"

fi

# download genome annotation
if [ ! -f "${REF_DIR}/potato_genome_annotation.v6.1.gtf" ]; then
    if [ ! -f "${REF_DIR}/potato_genome_annotation.v6.1.gff3" ]; then
        # genome annotation
        wget https://spuddb.uga.edu/data/dm_v61/DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz

        # decompress annotation
        gunzip DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz

        # move the annotation to their respective folders
        mv DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3 "${REF_DIR}/potato_genome_annotation.v6.1.gff3"
    fi
    # perform gffread to change gff to gtf for extract_splice_sites.py
    gffread "${REF_DIR}/potato_genome_annotation.v6.1.gff3" -T -o "${REF_DIR}/potato_genome_annotation.v6.1.gtf"
 
else
    # TODO MAKE A MORE INFORMATIVE PRINT STATEMENT
    echo "potato genome annotation has already been downloaded"
fi

# download GO term's from SpudDB

if [ ! -f "${REF_DIR}/potato_go_term2gene_all.tsv" ]; then
    wget https://spuddb.uga.edu/data/dm_v61/DM_1-3_516_R44_potato.v6.1.working_models.go.txt.gz

    gunzip DM_1-3_516_R44_potato.v6.1.working_models.go.txt.gz
    mv DM_1-3_516_R44_potato.v6.1.working_models.go.txt ${REF_DIR}/DM_1-3_516_R44_potato.v6.1.working_models.go.txt
    # retain only the gene ID, GO term and Aspect
    awk -F'\t' '$5 ~ /^GO:/ {g=$2; sub(/\.[0-9]+$/,"",g); print $5 "\t" g "\t" $9}' \
  ${REF_DIR}/DM_1-3_516_R44_potato.v6.1.working_models.go.txt | sort -u > ${REF_DIR}/potato_go_term2gene_all.tsv
fi

# split the GO terms by aspect
# BP == Biological Process
if [ ! -f "${REF_DIR}/term2gene_GO_BP.tsv" ]; then
    awk -F'\t' '$3=="P" {print $1 "\t" $2}' ${REF_DIR}/potato_go_term2gene_all.tsv > ${REF_DIR}/term2gene_GO_BP.tsv
fi

# F == Molecular function

if [ ! -f "${REF_DIR}/term2gene_GO_MF.tsv" ]; then
    awk -F'\t' '$3=="F" {print $1 "\t" $2}' ${REF_DIR}/potato_go_term2gene_all.tsv > ${REF_DIR}/term2gene_GO_MF.tsv
fi

# C == Cellular Component
if [ ! -f "${REF_DIR}/term2gene_GO_CC.tsv" ]; then
    awk -F'\t' '$3=="C" {print $1 "\t" $2}' ${REF_DIR}/potato_go_term2gene_all.tsv > ${REF_DIR}/term2gene_GO_CC.tsv
fi

# get gene expression matrix data as well

# download gene expression matrix from SpudDB

PROJ_DIR=$(pwd)/Transcriptomics
REF=$PROJ_DIR/ref

cd ${REF}
wget https://spuddb.uga.edu/data/gene_exp/DM_v6_1.working_models.gene_exp_matrix.txt.gz
gunzip DM_v6_1.working_models.gene_exp_matrix.txt.gz
mv DM_v6_1.working_models.gene_exp_matrix.txt potato_gene_expression_matrix.txt

awk -F'\t' '{print $1 "\t" $NF}' ${REF}/potato_gene_expression_matrix.txt > ${REF}/potato_transcript_key.tsv

