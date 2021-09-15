files=$1

# please deal with the damn line endings
dos2unix $files

sbatch sbatch_concat_fastqs.sh $files
