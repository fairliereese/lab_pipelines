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
#$ -N mapsample1

module load bowtie2
module load samtools

bowtie_ind=/data/users/freese/mortazavi_lab/ref/mm10/bowtie2/mm10
odir=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/
read1=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/trimmed_fastqs/pe_ATAC1_S1_R1.fastq
read2=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/trimmed_fastqs/pe_ATAC1_S1_R2.fastq

# sample 1
bowtie2 \
        -X 2000 \
	-p 16 \
        -x $bowtie_ind \
        -1 $read1 \
        -2 $read2 \
         > ${odir}ATAC1_S1.sam 
cat ${odir}ATAC1_S1.sam | samtools view -Sb - \
	| samtools sort - > ${odir}ATAC1_S1.bam
# get uniquely-mapped reads
samtools view -b -q 30 ${odir}ATAC1_S1.bam > ATAC1_S1_unique.bam	
