#!/bin/bash
#SBATCH -A SEYEDAM_LAB
#SBATCH -o submit_prod.o%A
#SBATCH -e submit_prod.e%A
#SBATCH --time=48:00:00
#SBATCH --job-name submit_files
#SBATCH --partition=standard

# use the following to submit files to the production server from a tsv

source ~/.bash_profile
conda activate encode_submissions

opref=$1

eu_register.py -m prod -p biosample -w -i ${opref}_biosample.tsv
eu_register.py -m prod -p experiment -w -i ${opref}_experiment.tsv
eu_register.py -m prod -p library -w -i ${opref}_library.tsv
eu_register.py -m prod -p replicate -w -i ${opref}_rep.tsv
eu_register.py -m prod -p file -w -i ${opref}_file.tsv
