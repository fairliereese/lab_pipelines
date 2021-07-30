files=$1
organism=$2

# please deal with the damn line endings
dos2unix $files

n=`wc -l $files | cut -d' ' -f1`
#sbatch --array=1-${n} sbatch_minimap.sh $files $organism
sbatch --array=1-1 sbatch_minimap.sh $files $organism
