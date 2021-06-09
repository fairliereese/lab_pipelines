#!/bin/sh
#SBATCH -A SEYEDAM_LAB
#SBATCH --cpus-per-task 64
#SBATCH --output=minimap.out
#SBATCH --error=minimap.err
#SBATCH --time=24:00:00
#SBATCH -J minimap
#SBATCH --mail-type=START,END
#SBATCH --partition=standard

config=$1

annot=$2
annot_name=$3

genome_name=$4

p_dir=talon/
mkdir -p ${p_dir}
talon_initialize_database \
    --f $annot \
    --g $genome_name \
    --a $annot_name \
    --l 0 \
    --idprefix ENCODEH \
    --5p 500 \
    --3p 300 \
    --o ${p_dir}pgp1

talon \
    --f $config \
    --db ${p_dir}pgp1.db \
    --build $genome_name \
    --t 64 \
    --o ${p_dir}pgp1
