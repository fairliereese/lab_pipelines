#!/bin/bash
#SBATCH --job-name=splitfastq    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-28               ## number of tasks to launch (number of samples)
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=fastqsplit.out ## output log file
#SBATCH --error=fastqsplit.err ## error log file

source activate seqtk
module load ucsc-tools
cd /share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/

sample=`head -n $SLURM_ARRAY_TASK_ID barcodes/prefixes.txt | tail -n 1`

# Grep unfiltered (no UMI or gene filter) cell IDs for each sample from barcoded trimmed fastq to make new sample-level fastqs
grep fastq/single_cells_barcoded_head.fastq -F -A 3 -f barcodes/${sample}.txt | grep -v '^--$' > fastq/${sample}.fastq

# Strip the barcodes
python3 barcodes/strip_sc_barcode.py --input fastq/${sample}.fastq --output fastq/${sample}_fixed.fastq

rm fastq/${sample}.fastq

# Get names from sample level fastqs
awk 'NR%4==1 {print substr($1,2)}' fastq/${sample}_fixed.fastq > fastq/${sample}.names.txt

rm fastq/${sample}_fixed.fastq

# Subset the raw, untrimmed fastq
seqtk subseq fastq/Adr_12k_R1.fastq fastq/${sample}.names.txt  > fastq/${sample}.fastq

# Check
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' fastq/${sample}.fastq
validateFiles -chromInfo=barcodes/mm10.chrom.sizes -type=fastq fastq/${sample}.fastq

# Zip
gzip fastq/${sample}.fastq
