#!/bin/bash
#SBATCH --job-name=lima
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
smrt_cell=`head -${i} $ifile | tail -1 | cut -f4`

# make directories
pb_dir=~/pacbio/$pb_id/
lima_dir=${pb_dir}Lima/

mkdir -p $lima_dir
cd $lima_dir

# get ccs reads for each data directory and run Lima
# for dir in ${pb_dir}CCS/*01/
# do
dir=${pb_dir}CCS/${smrt_cell}01/
files=($dir/ccs.bam)
bam=${files[0]}
name=$(basename "$dir" | cut -f1 -d"_")
data_dir=${lima_dir}${smrt_cell}01/
mkdir -p $data_dir

lima \
  ${bam} \
  /share/crsp/lab/seyedam/share/PACBIO/scripts/June2021/PB_adapters.fasta \
  ${data_dir}/fl.bam \
  --isoseq \
  --num-threads 32 \
  --min-score 0 \
  --min-end-score 0 \
  --min-signal-increase 10 \
  --min-score-lead 0

module unload bioconda/4.8.3

# because lima names the files weird, rename it
mv ${data_dir}/fl.primer_5p--primer_3p.bam ${data_dir}/fl.bam

echo "Finished Lima for $pb_id, $name"
n_reads=`samtools view -c ${data_dir}/fl.bam`
echo "$n_reads after Lima"
# done
