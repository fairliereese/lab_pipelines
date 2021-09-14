 docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-filter-beads:v1.0.0 \
 -i /data/alignments/ \
 -o /data/bead_filtration/ \
 -r mm10 \
 --name bead_filtration \
 --cpus="10"
