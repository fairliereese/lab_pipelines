#!/bin/sh

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_R1.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_R2.fastq.gz

done;

