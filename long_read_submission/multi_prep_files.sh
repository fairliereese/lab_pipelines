set -x
set -e


pb_ids=$1
pb_ids=$(echo $pb_ids | tr ',' '\n')
echo $pb_ids

name=$2

mkdir -p files/$name/raw_fastqs/
cd files/$name/raw_fastqs/
fastq=${name}.fastq

# concatenate flnc.fastq files 
for pb_id in $pb_ids
do
	cat /share/crsp/lab/seyedam/share/PACBIO/$pb_id/*Refine/*01/flnc*.fastq >> $fastq
done

# gzip and compute md5sum
gzip *fastq
md5sum *fastq.gz > md5sum

cd ../../../
