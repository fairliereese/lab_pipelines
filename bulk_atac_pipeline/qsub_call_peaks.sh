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
#$ -N call_peaks

module load macs2 

set -e

bamfile=$1

prefix=$(printf "$bamfile" | cut -d. -f1)
bname=$(basename $prefix)
odir=$(dirname $bamfile)

macs2 callpeak \
	-t $bamfile \
	-f BED \
	-n $prefix \
	--outdir "peaks" \
	--nolambda \
	--keep-dup all \
	--call-summits 
	
