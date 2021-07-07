files=$1

# fix carriage returns
# https://its.ucsc.edu/unix-timeshare/tutorials/clean-ctrl-m.html
# sed -e "s///" $files > temp
# mv temp $files

n=`wc -l $files | cut -d' ' -f1`
sbatch --array=1-${n} sbatch_download.sh $files
