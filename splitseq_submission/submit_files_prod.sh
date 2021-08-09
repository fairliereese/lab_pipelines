#!/bin/bash
#SBATCH -A SEYEDAM_LAB
#SBATCH -o submit_files.o%A
#SBATCH -e submit_files.e%A
#SBATCH --time=48:00:00
#SBATCH --job-name submit_files 
#SBATCH --partition=standard

# use the following to submit files to the production server from a tsv

source ~/.bash_profile
conda activate encode_submissions

sub_tab=$1

eu_register.py -m prod -p file -w -i $1
