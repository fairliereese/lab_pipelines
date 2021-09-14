#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N demultiplex
#$ -m ea
#$ -cwd
#$ -j y

fq=$1
i_bcs=$2
opref=$3
d=$4

python ${d}demultiplex.py \
    -f $fq \
    -o $opref \
    -t 16 \
    -i_file $i_bcs \
    -rc 500