#!/bin/bash
#SBATCH --job-name=dl_pb
#SBATCH -A SEYEDAM_LAB
#SBATCH -o dl_files/%x.o%A_%a
#SBATCH -e dl_files/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mail-user=freese@uci.edu

set -x
set -e

# input tab-separated file with PB ID \t file location at GHTF \t md5sum
ifile=$1

# extract PBID, file location, and md5sum
i=$SLURM_ARRAY_TASK_ID
pb_id=`head -${i} $ifile | tail -1 | cut -f1`
link=`head -${i} $ifile | tail -1 | cut -f2`
md5sum=`head -${i} $ifile | tail -1 | cut -f3`

# make directories and download
mkdir -p ~/pacbio/$pb_id/A01_data/
cd ~/pacbio/$pb_id/A01_data/
wget $link

# check md5sum
test_md5sum=`md5sum *.subreads.bam | cut -d' ' -f1`

if [ "${test_md5sum}" == "${md5sum}" ]; then
	echo "md5sums match"
else
	echo "MD5SUM CONTENT ERROR"
	echo "md5sum from GHTF: $md5sum"
	echo "md5sum of downloaded file: $test_md5sum"
fi
