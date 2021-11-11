# places to look for runs
dir1=/run/media/nextseq2000/Elements/
dir2=/big_backup2/

# places to move runs
modelad=/run/media/nextseq2000/Elements/
tenner=/run/media/nextseq2000/tenner/

# tenner runs
tenner_runs=211110_tenner_files.tsv

# model ad runs
modelad_runs=211110_modelad_files.tsv

# log file
log=211110_log.txt
touch $log
rm $log
touch $log

# first tenner lab runs
echo "T E N N E R"
runs=$tenner_runs
dos2unix $runs
backup_dirs=( $dir1 $dir2 )
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

# then model ad runs
echo "M O D E L A D"
runs=$modelad_runs
dos2unix $runs
backup_dirs=( $dir1 $dir2 )
dest_dir=$modelad
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
