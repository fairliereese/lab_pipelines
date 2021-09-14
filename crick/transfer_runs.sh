d1=/big_backup2/
d2=/big_backup/
d3=/run/media/nargesr/Elements/
model_ad=/run/media/nextseq2000/MODEL-AD/
runs=/runs/210816_transfer_runs.tsv
runs_with_errors=/runs/210816_transfer_runs_error.tsv

touch $runs_with_errors

dos2unix $runs


backup_dirs=( $d1 $d2 $d3 )
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
        cp $f $model_ad
        fname=$(basename $f)
        old_size=`ls -lat ${f} | cut -d' ' -f5`
        new_size=`ls -lat ${model_ad}${fname} | cut -d' ' -f5`

        if [ "${old_size}" == "${new_size}" ]; then
        	echo "Files are the same size in ${d1} and ${model_ad}, deleting from ${runs_dir}"
        	rm -f ${f}
        else
        	echo "Files are not the same size. Not deleting" >> $runs_with_errors
        	echo "Size of ${f}: $old_size" >> $runs_with_errors
        	echo "Size of ${model_ad}${fname}: $new_size" >> $runs_with_errors
          printf "${run}\n" >> $runs_with_errors
        fi
      fi
    done
  done

  printf "\n"

done < $runs
