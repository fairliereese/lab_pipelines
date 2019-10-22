# from https://github.com/r3fang/SnapATAC/wiki/FAQs#cellranger_output
bam=$1
prefix=${bam%.bam}

snaptools snap-pre  \
	--input-file=${bam}  \
	--output-snap=${bam}.snap  \
	--genome-name=mm10  \
	--genome-size=mm10.chrom.sizes  \
	--min-mapq=30  \
	--min-flen=50  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=20000  \
	--min-cov=500  \
	--verbose=True