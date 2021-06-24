#!/bin/sh
#SBATCH -A SEYEDAM_LAB
#SBATCH --output=talon_label.out
#SBATCH --error=talon_label.err
#SBATCH --time=10:00:00
#SBATCH -J talon_label
#SBATCH --mail-type=START,END
#SBATCH --partition=standard

# input arguments
samples=$1
genome=$2

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

sam=${p}_clean.sam
talon_label_reads \
    --f $sam \
    --g $genome \
    --tmpDir ${p}_tmp/ \
    --ar 20  \
    --deleteTmp  \
    --o $p
