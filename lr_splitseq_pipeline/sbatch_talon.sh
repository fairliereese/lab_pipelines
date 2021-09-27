#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N TALON
#$ -m ea
#$ -cwd
#$ -j y

opref=$1
sample=$2
sam=${opref}_merged_primers.sam
gtf=~/mortazavi_lab/ref/gencode.vM21/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf

talon_initialize_database \
    --f ${gtf} \
    --g mm10 \
    --a gencode_vM21 \
    --l 0 \
    --idprefix ENCODEM \
    --5p 500 \
    --3p 300 \
    --o ${opref}

printf "${sample},SequelII,${sam}" > ${opref}_config.csv

talon \
    --f ${opref}_config.csv \
    --cb \
    --db ${opref}.db \
    --build mm10 \
    -t 16 \
    --o ${opref}
