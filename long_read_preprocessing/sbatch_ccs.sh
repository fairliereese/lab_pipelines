#!/bin/bash
#SBATCH --job-name=ccs
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
ccs_dir=${pb_dir}CCS/

mkdir -p $ccs_dir
cd $ccs_dir

# get subreads for each data directory and run CCS
for dir in ${pb_dir}/*_data
do
  files=($dir/*subreads.bam)
  subreads=${files[0]}
  name=$(basename "$dir" | cut -f1 -d"_")
  data_dir=${ccs_dir}${name}
  mkdir -p $data_dir

  ccs \
    --skip-polish \
    --min-length=10 \
    --min-passes=3 \
    --min-rq=0.9 \
    --min-snr=2.5 \
    --report-file ${data_dir}/ccs_report.txt $subreads ${data_dir}/ccs.bam

    echo "Finished CCS for $pb_id, $name"
    n_reads=`samtools view -c ${data_dir}/ccs.bam`
    echo "$n_reads after CCS"
done
