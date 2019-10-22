# from https://github.com/r3fang/SnapATAC/wiki/FAQs#cellranger_output
snap=$1
prefix=${snap%.snap}
echo $snap
echo $prefix

snaptools snap-add-bmat  \
	--snap-file=${snap}  \
	--bin-size-list 50  \
	--verbose=True