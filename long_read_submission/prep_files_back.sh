# use the following to prepare both subreads and fastqs for submission
# called from sbatch_prep_files.sh

set -x
set -e 

pb_ids=$1
pb_ids=$(echo $pb_ids | tr ',' '\n')
echo $pb_ids

name=$2


# fastqs
mkdir -p files/$name/raw_fastqs/
cd files/$name/raw_fastqs/

fastq=flnc.fastq

# concatenate flnc.fastq files 
for pb_id in $pb_ids
do
	d=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/Refine/
	bam_files=${d}*0*/flnc*.bam
	for bam_file in $bam_files
	do
		bam_file=`ls $bam_file`
		path=$(dirname "${bam_file}")
		if test -f "${path}/flnc.fastq"; then
			echo "flnc.fastq already exists, won't run bam2fq"
		else
			module load samtools
			samtools bam2fq $bam_file > ${path}/flnc.fastq
		fi
		cat ${path}/flnc.fastq >> $fastq
	done
done

# # copy the file to the submission dir
# cp /share/crsp/lab/seyedam/share/PACBIO/$pb_id/*Refine/*01/flnc*.fastq .

# gzip and compute md5sum
gzip $fastq
md5sum flnc.fastq.gz > md5sum

cd ../../../

mkdir -p files/$name/subreads/
cd files/$name/subreads/

# subreads
touch md5sum
for pb_id in $pb_ids
do
	#subreads=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/*01_data/*subreads.bam
	# i had to use the following line for weird dir structure for PB197_run2
	subreads=/share/crsp/lab/seyedam/share/PACBIO/$pb_id/*0*_dat*/*subreads.bam
	for f in $subreads
	do
		# gzip $f
		ln -s ${f} .
		# md5sum ${f} >> md5sum
	done
done

cd ../../../
