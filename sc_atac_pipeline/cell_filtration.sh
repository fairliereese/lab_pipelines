docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-cell-filter:v1.0.0 \
 -i /data/deconvoluted_data/ \
 -o /data/cells_filtered/ \
 -r mm10 \
 --name cell_filtration
