#!/bin/bash
set -euo pipefail

# establish parameters?

# read from a text file all of the sample numbers to be run

############## Step 0: Gather Reference Genome if Needed ##############
# make sure that the tools for this script are enabled
#conda init
#conda deactivate
#conda activate 6701-rnaseq

# read all samples to be processed
#SRA_LIST=$1

# input index variable with param

# system variable for building index with splice sites and stuff
PROJ_DIR=$(pwd) # run within the directory

# create directory
REF_DIR="${PROJ_DIR}/transcriptomics/ref"
if [ ! -d ${REF_DIR} ]; then
    mkdir ${REF_DIR}
fi

echo "${REF_DIR}"
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

############## Step 1: Perform Index Build ##############

# Fast works by not getting known splice sites and exon sites for indexing

# build index, but first check if the indexes are already present
build=0
for i in {1..8}; do
    if [ ! -f "${REF_DIR}/potato_dm_v6.1.${i}.ht2" ]; then
        build=1
    fi
done

if [ $build == 1 ]; then
    # build the index if the above condition is met
    echo
    echo "performing hisat2 build"
    echo
    hisat2-build -p 4 "${REF_DIR}/potato_genome_assembly.v6.1.fa" "${REF_DIR}/potato_dm_v6.1"

else
    echo "index has already been built for potato genome"

fi