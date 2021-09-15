#!/bin/bash
#SBATCH --job-name=demux
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

# extract PBID
i=$SLURM_ARRAY_TASK_ID
pb_id=`head -${i} $ifile | tail -1 | cut -f1`
smrt_cell=`head -${i} $ifile | tail -1 | cut -f4`

# make directories
pb_dir=~/pacbio/$pb_id/

# get flnc post-Refine reads for each data directory
# for dir in ${pb_dir}Refine/*01/
# do
dir=${pb_dir}Refine/${smrt_cell}01/
files=($dir/flnc.fastq)
fastq=${files[0]}
name=$(basename "$dir" | cut -f1 -d"_")

fq=$1
i_bcs=$2
opref=$3
d=$4

python ~/mortazavi_lab/bin/lab_pipelines/bin/pacbio-splitpip/demultiplex.py \
    -f ${fastq} \
    -o ${opref}${name} \
    -t 32 \
    -rc 500
