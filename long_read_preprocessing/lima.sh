files=$1

n=`wc -l $files | cut -d' ' -f1`
sbatch --array=1-${n} sbatch_lima.sh $files
# sbatch --array=1-1 sbatch_lima.sh $files
