#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N TALON
#$ -m ea
#$ -cwd
#$ -j y

d=
f=$1
sample=$2
gtf=$3
oprefix=$4

talon_initialize_database \
    --f ${gtf} \
    --g mm10 \
    --a gencode_vM21 \
    --l 0 \
    --idprefix ENCODEM \
    --5p 500 \
    --3p 300 \
    --o ${oprefix}

printf "${sample},SequelII,${f}" > ${oprefix}_config.csv

talon \
    --f ${oprefix}_config.csv \
    --cb \
    --db ${oprefix}.db \
    --build mm10 \
    -t 16 \
    --o ${oprefix}