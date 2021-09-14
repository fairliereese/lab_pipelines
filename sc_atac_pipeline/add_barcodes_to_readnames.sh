# from https://github.com/r3fang/SnapATAC/wiki/FAQs#cellranger_output
bam=$1
prefix=${bam%.bam}
echo $bam
echo $prefix

# extract the header file
samtools view ${bam} -H > ${prefix}.header.sam

# create a bam file with the barcode embedded into the read name
cat <( cat ${prefix}.header.sam ) \
<( samtools view ${bam} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^RG:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["RG"], $0 }' ) \
| samtools view -bS - > ${prefix}.snap.bam

# sort 
samtools sort -n -@ 10 -m 1G ${prefix}.snap.bam -o ${prefix}.snap.nsrt.bam

