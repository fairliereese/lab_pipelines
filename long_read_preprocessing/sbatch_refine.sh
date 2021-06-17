#!/bin/bash
#SBATCH --job-name=refine
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A_%a
#SBATCH -e processing_tables/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=7-0
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

module load bioconda/4.8.3
module load samtools

set -x
set -e

# input tab-separated file with PB ID \t file location at GHTF \t md5sum
ifile=$1

# extract PBID
i=$SLURM_ARRAY_TASK_ID
pb_id=`head -${i} $ifile | tail -1 | cut -f1`

# make directories
pb_dir=~/pacbio/$pb_id/
refine_dir=${pb_dir}Refine/

mkdir -p $refine_dir
cd $refine_dir

# get fl post-Lima reads for each data directory and run Lima
for dir in ${pb_dir}Lima/*01/
do
  files=($dir/fl.bam)
  bam=${files[0]}
  name=$(basename "$dir" | cut -f1 -d"_")
  data_dir=${refine_dir}${name}
  mkdir -p $data_dir

  isoseq3 refine \
    ${bam}  \
    /share/crsp/lab/seyedam/share/PACBIO/scripts/June2021/PB_adapters.fasta\
    ${data_dir}/flnc.bam \
    --num-threads 32 \
    --require-polya

  # create fastq as well
  samtools bam2fq ${data_dir}/flnc.bam > ${data_dir}/flnc.fastq

  echo "Finished Refine for $pb_id, $name"
  n_reads=`samtools view -c ${data_dir}/flnc.bam`
  echo "$n_reads after Refine"
done
