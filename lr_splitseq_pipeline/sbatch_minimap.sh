#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A
#SBATCH -e processing_tables/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

fastq="${opref}_demux.fastq"
genome=~/mortazavi_lab/ref/mm10/mm10.fa
sam=${opref}_mapped.sam
log=${opref}_minimap.log


module load minimap2

minimap2  \
    -t 32 \
    -ax splice:hq \
    -uf \
    --MD \
    $ref \
    $fastq \
    > $sam \
    2> $log
