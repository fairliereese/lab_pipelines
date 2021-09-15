files=$1
opref=$2

# please deal with the damn line endings
dos2unix $files

sbatch sbatch_concat_fastqs.sh $files $opref
