files=$1

n=`wc -l $files | cut -d' ' -f1`
sbatch --array=1-${n} lima_sbatch.sh $files
# sbatch --array=1-1 lima_sbatch.sh $files
