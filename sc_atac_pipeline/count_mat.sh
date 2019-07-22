docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-chromvar:v1.0.0 \
 -d /data/cells_filtered \
 -p /data/peaks \
 -o /data/count_matrix \
 -r mm10 \
 --name count_matrix
