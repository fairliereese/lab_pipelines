#!/bin/sh
#SBATCH -A SEYEDAM_LAB
#SBATCH --cpus-per-task 16
#SBATCH --output=tc.out
#SBATCH --error=tc.err
#SBATCH --time=12:00:00
#SBATCH -J tc
#SBATCH --mail-type=START,END
#SBATCH --partition=standard

module load samtools

# input arguments
samples=$1
tc_path=$2

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

# input and output filenames
sam=${p}_mapped.sam
bam=${p}_mapped.bam
sort_bam=${p}_sorted.bam
sort_sam=${p}_sorted.sam

# first sort the sam file
samtools view -Sb $sam > $bam
samtools sort $bam > $sort_bam
samtols view -h $sort_bam > $sort_sam

# run TranscriptClean
python ${tc_path}/TranscriptClean.py \
   --sam $sort_sam \
   --genome $genome \
   -t 16 \
   --canonOnly \
   --deleteTmp \
   --outprefix $p
