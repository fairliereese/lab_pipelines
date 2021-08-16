f=$1 # sam /bam file
id=$2 # id of dataset

module load samtools

samtools view -hF 256 $f | awk '{print length($10)}' > ${i}_read_lengths.txt

python calc_n50.py \
  --f ${i}_read_lengths.txt \
  --id ${id}

rm ${i}_read_lengths.txt
