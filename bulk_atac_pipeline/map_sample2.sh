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
#$ -N mapsample2

module load bowtie2
module load samtools

bowtie_ind=/data/users/freese/mortazavi_lab/ref/mm10/bowtie2/mm10
odir=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/
read1=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/trimmed_fastqs/ATAC2_S2_R1_pe.fastq
read2=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/trimmed_fastqs/ATAC2_S2_R2_pe.fastq

# sample 1
bowtie2 \
        -X 2000 \
		-p 16 \
        -x $bowtie_ind \
        -1 $read1 \
        -2 $read2 \
         > ${odir}ATAC2_S2.sam 
cat ${odir}ATAC2_S2.sam | samtools view -u - \
	| samtools sort - > ${odir}ATAC2_S2.bam
# get uniquely-mapped reads
samtools view -b -q 30 ${odir}ATAC2_S2.bam > ATAC2_S2_unique.bam
