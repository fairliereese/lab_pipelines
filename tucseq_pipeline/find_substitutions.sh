sam=$1
opref=$2

# parse CIGAR and MD tags to return alignments of each read

# the following block is if the sam file has a header
# module load samtools
# samtools view $sam | sam2pairwise > ${opref}_temp.sam

# the following line is if the sam file does not have a header
cat $sam | sam2pairwise > ${opref}_temp.sam


# parse alignments to find counts of each type of substitution
# in each read
python find_substitutions.py \
    -sam ${opref}_temp.sam \
    -opref $opref

# clean up temp files
rm ${opref}_temp.sam
rm ${opref}_mistmatches.tsv
