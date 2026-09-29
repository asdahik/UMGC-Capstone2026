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

PROJ_DIR=$(pwd) # run within the directory

# create a data directory
if [ ! -d "data" ]; then
    bash ${PROJ_DIR}/scripts/00_make_dir.bash
fi
RAWDATA_DIR="${PROJ_DIR}/data/raw"
REF_DIR="${PROJ_DIR}/data/reference"
# Check if all data has been downloaded for this particular project set
# need to get reference genome. Using DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3 and DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

# check if genome assembly and annotation have been downloaded


# downlaod genome assembly
if [ ! -f "${PROJ_DIR}/data/ref/potato_genome_assembly.v6.1.fa" ]; then

    # genome assembly
    wget https://spuddb.uga.edu/data/dm_v61/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

    # decompress genome
    gunzip DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz

    # move the assembly to their respective folders
    mv DM_1-3_516_R44_potato_genome_assembly.v6.1.fa "${PROJ_DIR}/data/ref/potato_genome_assembly.v6.1.fa"

else
    # TODO MAKE A MORE INFORMATIVE PRINT STATEMENT
    echo "potato genome assembly has already been downloaded"

fi

# download genome annotation
if [ ! -f "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gtf" ]; then
    if [ ! -f "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gff3" ]; then
        # genome annotation
        wget https://spuddb.uga.edu/data/dm_v61/DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz

        # decompress annotation
        gunzip DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3.gz

        # move the annotation to their respective folders
        mv DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3 "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gff3"
    fi
    # perform gffread to change gff to gtf for extract_splice_sites.py
    gffread -T -o "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gff3" "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gff3"

else
    # TODO MAKE A MORE INFORMATIVE PRINT STATEMENT
    echo "potato genome annotation has already been downloaded"
fi

############## Step 1: Perform Index Build ##############

# we have splice sites and exon sites mapped but keep them separate from index build for now
if [ ! -f "${PROJ_DIR}/data/ref/potato_dm_v6.1.ss" ]; then
    
    # python script from HISAT2 for getting splice sites .ss
    echo
    echo "performing HISAT2 index building: generating splice sites"
    extract_splice_sites.py --verbose "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gtf" > potato_dm_v6.1.ss
    mv potato_dm_v6.1.ss "${PROJ_DIR}/data/ref/potato_dm_v6.1.ss"

else
    echo
    echo "HISAT2 splice sites already generated!" 

fi


if [ ! -f "${PROJ_DIR}/data/ref/potato_dm_v6.1.exon" ]; then
    
    # python script from HISAT2 for getting exon sites .exon
    echo
    echo "performing HISAT2 index building: generating exon sites"
    extract_exons.py --verbose "${PROJ_DIR}/data/ref/potato_genome_annotation.v6.1.gtf" > potato_dm_v6.1.exon
    mv potato_dm_v6.1.exon "${PROJ_DIR}/data/ref/potato_dm_v6.1.exon"

else
    echo
    echo "HISAT2 exon sites already generated!" 
fi

# build index, but first check if the indexes are already present
build=0
for i in {1..8}; do
    if [ ! -f "${PROJ_DIR}/data/ref/potato_dm_v6.1.${i}.ht2" ]; then
        build=1
    fi
done

if [ $build == 1 ]; then
    # build the index if the above condition is met
    echo
    echo "performing hisat2 build"
    echo
    hisat2-build -p 4 data/ref/potato_genome_assembly.v6.1.fa \
                data/ref/potato_dm_v6.1

else
    echo "index has already been built for potato genome"

fi
# begin cleaning using fastq with fastqc and then aggregate with multiqc?
# trim data with fastp 

# once all data has been processed, run HISAT2
# make sure that reference genome is specified, better to download

