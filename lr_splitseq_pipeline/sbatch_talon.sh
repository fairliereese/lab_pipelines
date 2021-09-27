#!/bin/bash
#SBATCH --job-name=talon
#SBATCH -n 16
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A
#SBATCH -e processing_tables/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=32G
#SBATCH --mail-user=freese@uci.edu

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
