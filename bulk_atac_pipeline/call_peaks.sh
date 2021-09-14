#!/bin/bash

module load macs2 

set -e

bamfile=$1

prefix=$(printf "$bamfile" | cut -d. -f1)
bname=$(basename $prefix)
odir=$(dirname $bamfile)

macs2 callpeak \
	-t $bamfile \
	-n $bname \
	--outdir "peaks" \
	--nolambda \
	--keep-dup all \
	--call-summits 
	
