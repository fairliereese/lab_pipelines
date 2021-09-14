#!/bin/bash
#$ -q som,bio,free64
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -ckpt restart
#$ -N trim

adapters=/data/users/freese/mortazavi_lab/ref/nextera_adapters/adapters.fa
combined_fastqs=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/combined_fastqs/
trimmed_fastqs=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/trimmed_fastqs/

# sample 1
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE ${combined_fastqs}ATAC1_S1_R1.fastq ${combined_fastqs}ATAC1_S1_R2.fastq ${trimmed_fastqs}pe_ATAC1_S1_R1.fastq ${trimmed_fastqs}se_ATAC1_S1_R1_.fastq.gz ${trimmed_fastqs}pe_ATAC1_S1_R2.fastq ${trimmed_fastqs}se_ATAC1_S1_R2.fastq.gz ILLUMINACLIP:${adapters}:2:30:8:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:17 MINLEN:30

# sample 2
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE ${combined_fastqs}ATAC2_S2_R1.fastq ${combined_fastqs}ATAC2_S2_R2.fastq ${trimmed_fastqs}ATAC2_S2_R1_pe.fastq ${trimmed_fastqs}ATAC2_S2_R1_se.fastq.gz ${trimmed_fastqs}ATAC2_S2_R2_pe.fastq ${trimmed_fastqs}ATAC2_S2_R2_se.fastq.gz ILLUMINACLIP:${adapters}:2:30:8:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:17 MINLEN:30

