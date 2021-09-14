# from https://github.com/r3fang/SnapATAC/wiki/FAQs#cellranger_output
snap=$1
binsize=$2
prefix=${snap%.snap}
echo $snap
echo $prefix
echo $binsize

snaptools snap-add-bmat  \
	--snap-file=${snap}  \
	--bin-size-list ${binsize}  \
	--verbose=True