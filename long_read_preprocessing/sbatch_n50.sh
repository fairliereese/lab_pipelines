#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A_%a
#SBATCH -e processing_tables/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=7-0
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

module load minimap2/2.17
module load samtools

set -x
set -e

ifile=$1

# extract PBID
i=$SLURM_ARRAY_TASK_ID
pb_id=`head -${i} $ifile | tail -1 | cut -f1`
smrt_cell=`head -${i} $ifile | tail -1 | cut -f4`

# get flnc post-Minimap reads for each data directory
dir=${pb_dir}Minimap/${smrt_cell}01/
files=($dir/mapped_flnc.sam)
sam=${files[0]}

bash calc_n50.sh \
  $sam \
  $pb_id
