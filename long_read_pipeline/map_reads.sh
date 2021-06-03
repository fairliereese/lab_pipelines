#!/bin/sh
#SBATCH -A SEYEDAM_LAB
#SBATCH --cpus-per-task 16
#SBATCH --output=minimap.out
#SBATCH --error=minimap.err
#SBATCH --time=06:00:00
#SBATCH -J minimap
#SBATCH --mail-type=START,END
#SBATCH --partition=standard

module load samtools
module load minimap2/2.17

# input arguments
samples=$1
genome=$2

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

# run minimap on this sample
fastq=${p}.fastq
sam=${p}_mapped.sam
log=${p}_minimap.log
minimap2 \
    -t 16 \
    -ax splice:hq \
    -uf \
    --MD \
    --secondary=no \
    $genome \
    $fastq > \
    $sam 2> \
    $log
