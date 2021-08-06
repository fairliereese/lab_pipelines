files=$1

dos2unix $files

cat $files | while read p
do
	echo $p
	for f in ~/pacbio/${p}/*01_data/*subreads*
	do
		echo "Removing $f...\n"
		rm $f
	done
done
