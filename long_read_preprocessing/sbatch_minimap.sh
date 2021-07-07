#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing_tables/%x.o%A_%a
#SBATCH -e processing_tables/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=7-0
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

module load minimap2/2.17
module load samtools

set -x
set -e

# input tab-separated file with PB ID \t file location at GHTF \t md5sum
ifile=$1
organism=$2 # human or mouse

if [ "${organism}" == "human" ]; then
	genome=~/mortazavi_lab/ref/hg38/hg38.fa
elif [ "${organism}" == "mouse" ]; then
	genome=~/mortazavi_lab/ref/mm10/mm10.fa
fi

# extract PBID
i=$SLURM_ARRAY_TASK_ID
pb_id=`head -${i} $ifile | tail -1 | cut -f1`

# make directories
pb_dir=~/pacbio/$pb_id/
minimap_dir=${pb_dir}Minimap/

mkdir -p $minimap_dir
cd $minimap_dir

# get flnc post-Refine reads for each data directory and run Minimap
for dir in ${pb_dir}Refine/*01/
do
  files=($dir/flnc.fastq)
  fastq=${files[0]}
  name=$(basename "$dir" | cut -f1 -d"_")
  data_dir=${minimap_dir}${name}
  mkdir -p $data_dir

  minimap2 \
    -t 10 \
    -ax splice:hq \
    -uf \
    --MD \
    --secondary=no \
    ${genome}~/PACBIO/Sequel2_scripts/mm10_SIRV_ERCC_set4.fasta \
    ${fastq} \
    > ${data_dir}/mapped_flnc.sam \
    2> ${data_dir}/mapped_flnc.log

  echo "Finished Minimap for $pb_id, $name"
  n_reads=`samtools view -c ${data_dir}/mapped_flnc.sam`
  echo "$n_reads after Minimap"
done
