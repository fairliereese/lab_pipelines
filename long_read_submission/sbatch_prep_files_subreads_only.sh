#!/bin/sh
#SBATCH -A SEYEDAM_LAB
#SBATCH -o %x.o%A_%a
#SBATCH -e %x.e%A_%a
#SBATCH --time=06:00:00
#SBATCH --job-name prep_files
#SBATCH --partition=standard

file=$1

set -x
set -e

# use the following to batch run prep files, for both subreads and flnc fastqs'
pb_ids=`head -${SLURM_ARRAY_TASK_ID} $file  | tail -1 | cut -f1`

echo "PB ids"
echo $pb_ids
name=`head -${SLURM_ARRAY_TASK_ID} $file  | tail -1 | cut -f2`
echo "name"
echo $name

pb_ids=$(echo $pb_ids | tr ',' '\n')
echo $pb_ids

# # fastqs
# mkdir -p files/$name/raw_fastqs/
# cd files/$name/raw_fastqs/
#
# fastq=flnc.fastq
#
# # concatenate flnc.fastq files
# for pb_id in $pb_ids
# do
# 	d=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/Refine/
# 	bam_files=${d}*0*/flnc*.bam
# 	for bam_file in $bam_files
# 	do
# 		bam_file=`ls $bam_file`
# 		path=$(dirname "${bam_file}")
# 		if test -f "${path}/flnc.fastq"; then
# 			echo "flnc.fastq already exists, won't run bam2fq"
# 		else
# 			module load samtools
# 			samtools bam2fq $bam_file > ${path}/flnc.fastq
# 		fi
# 		cat ${path}/flnc.fastq >> $fastq
# 	done
# done
#
# # gzip and compute md5sum
# gzip $fastq
# md5sum flnc.fastq.gz > md5sum
#
# cd ../../../

mkdir -p files/$name/subreads/
cd files/$name/subreads/

# subreads
touch md5sum
for pb_id in $pb_ids
do
	subreads=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/*01_data/*subreads.bam
	# i had to use the following line for weird dir structure for PB197_run2
	# subreads=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/*0*_dat*/*subreads.bam
	for f in $subreads
	do
		ln -s ${f} .
	done
done

cd ../../../
