#!/bin/bash
#SBATCH -A SEYEDAM_LAB
#SBATCH -o submit_dev.o%A
#SBATCH -e submit_dev.e%A
#SBATCH --time=8:00:00
#SBATCH --job-name submit_files
#SBATCH --partition=standard

set -e

# use the following to submit files to the production server from a tsv

source ~/.bash_profile
conda activate encode_submissions

opref=$1
biosamp=$2
file=$3

# 1 = submit biosample
# 0 = do not submit biosample
if [ $biosamp -eq 1 ]
then
  eu_register.py -m dev -p biosample -w -i ${opref}_biosample.tsv
fi

eu_register.py -m dev -p experiment -w -i ${opref}_experiment.tsv
eu_register.py -m dev -p library -w -i ${opref}_library.tsv
eu_register.py -m dev -p replicate -w -i ${opref}_rep.tsv

# 1 = submit files
# 0 = do not submit files
if [ $file -eq 1 ]
then
  eu_register.py -m dev -p file -w -i ${opref}_r1_fastq.tsv
  eu_register.py -m dev -p file -w -i ${opref}_r2_fastq.tsv
fi
