#!/bin/bash
#SBATCH --job-name=demux
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A
#SBATCH -e processing_tables/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=128G
#SBATCH --mail-user=freese@uci.edu

opref=$1
fastq=${opref}.fastq
python ~/mortazavi_lab/bin/LR-splitpipe/LR-splitpipe/demultiplex.py \
    -f ${fastq} \
    -o ${opref} \
    -t 32 \
    -rc 0
