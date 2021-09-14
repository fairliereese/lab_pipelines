#!/bin/bash

module load samtools

bamfile=$1
prefix=`printf "$bamfile" | cut -d. -f1`
set -e 

# remove mitochondrial reads
printf "Removing mitochondrial reads for ${bamfile}...\n"
readnum=($(wc -l $bamfile))
printf "reads before removing chrM reads: $readnum\n" > ${prefix}mit_reads
samtools view -h $bamfile \
	| cut -f3 \
	| grep -v "chrM" \
	| samtools view -b \
	> temp.bam
filtnum=($(wc -l temp.bam))
printf "reads after removing chrM reads: $filtnum\n" >> ${prefix}mit_reads
#pct=$((filtnum/readnum))
#pct=$((pct*100))
#printf "percent mit. reads: $pct"
printf "\n"

# remove PCR duplicates
printf "Removing PCR duplicates for ${bamfile}...\n"
module load picard-tools/1.96
java -Xmx2g -jar /data/apps/picard-tools/1.96/MarkDuplicates.jar \
	INPUT=temp.bam \
	OUTPUT=${prefix}_filtered.bam \
	METRICS_File= ${prefix}_picard_metrics \
	REMOVE_DUPLICATES=true

# remove temp file
rm temp.bam

# index
samtools index ${prefix}_filtered.bam

