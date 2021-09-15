#!/bin/bash
#SBATCH --job-name=concat_fastqs
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A_%a
#SBATCH -e processing_tables/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=7-0
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

# input tab-separated file with PB ID \t file location at GHTF \t md5sum
ifile=$1
opref=$2

while read sample
do
  # extract PBID
  i=$SLURM_ARRAY_TASK_ID
  pb_id=`cut -f1 ${sample}`
  smrt_cell=`cut -f4 ${sample}`

  # make directories
  pb_dir=~/pacbio/$pb_id/

  # get flnc post-Refine reads for each data directory
  dir=${pb_dir}Refine/${smrt_cell}01/
  files=($dir/flnc.fastq)
  fastq=${files[0]}
  echo $fastq
done < ${ifile}
