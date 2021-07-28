files=$1

cat $files | while read p
do
	echo $p
	for f in ~/pacbio/${p}/*01_data/*subreads*
	do
		echo $f
	done
done
