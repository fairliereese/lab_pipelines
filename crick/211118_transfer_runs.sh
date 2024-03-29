# places to look for runs
tenner=/big_backup2/

# places to move runs
dir1=/run/media/nextseq2000/9ac0024d-1923-4e2b-a6b8-43e69a092d1c/

# tenner runs
tenner_runs=211121_transfer_runs.txt

# log file
log=211121_log.txt
touch $log
rm $log
touch $log

runs=$tenner_runs
dos2unix $runs
backup_dirs=( $dir1 )
dest_dir=$tenner
while read run
do
  echo "Run number ${run}"

  # loop through different backup dirs
  for d in "${backup_dirs[@]}"
  do
    echo "Looking in ${d}..."
    for f in ${d}*NS500169_0${run}_*
    do
      if test -f "$f"; then
        echo "Copying file ${f}..."
        cp $f $dest_dir
        fname=$(basename $f)
        old_size=`ls -lat ${f} | cut -d' ' -f5`
        new_size=`ls -lat ${dest_dir}${fname} | cut -d' ' -f5`

        if [ "${old_size}" == "${new_size}" ]; then
        	echo "Files are the same size in ${d} and ${dest_dir}, deleting from ${runs_dir}"
        	rm -f ${f}
        else
        	echo "Files are not the same size. Not deleting" >> $log
        	echo "Size of ${f}: $old_size" >> $log
        	echo "Size of ${dest_dir}${fname}: $new_size" >> $log
          printf "${run}\n" >> $log
        fi
      fi
    done
  done
  printf "\n"
done < $runs
