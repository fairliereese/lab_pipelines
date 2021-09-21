#!/bin/bash
#SBATCH --job-name=tc
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A
#SBATCH -e processing_tables/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

module load samtools

opref=$1
sam=${opref}_mapped.sam
sam_noscaff=${opref}_sorted_no_scaff.sam

# references
tc_path=/dfs6/pub/freese/mortazavi_lab/bin/TranscriptClean/
genome=/dfs6/pub/freese/mortazavi_lab/ref/mm10/mm10.fa
sjs=/data/users/freese/mortazavi_lab/ref/mm10/mm10_SJs.tsv


# remove reads that mapped to scaffold chromosomes and sort
grep -v '^@' $sam | awk ""'{if($3 !~ "_") print $0}'" "| \
          cat <(samtools view -H $sam) - | samtools view -bS | \
          samtools sort | samtools view -h > $sam_noscaff

# run tc
time python ${tc_path}TranscriptClean.py \
    --sam $sam_noscaff \
    --genome $genome \
    --spliceJns $sjs \
    -t 4 \
    --canonOnly \
    --primaryOnly \
    --deleteTmp \
    --outprefix ${opref}
