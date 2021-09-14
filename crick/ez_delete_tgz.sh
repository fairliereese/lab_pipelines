fname=$1
runs_dir=/runs/
backup_dir=/big_backup2/

size_runs=`ls -lat ${runs_dir}${fname} | cut -d' ' -f5`
size_backup=`ls -lat ${backup_dir}${fname} | cut -d' ' -f5`

if [ "${size_runs}" == "${size_backup}" ]; then
	echo "Files are the same size in ${runs_dir} and ${backup_dir}, deleting from ${runs_dir}"
	rm -f ${runs_dir}${fname}
else
	echo "Files are not the same size. Not deleting"
	echo "Size of ${runs_dir}${fname}: $size_runs"
	echo "Size of ${backup_dir}${fname}: $size_backup"
fi
