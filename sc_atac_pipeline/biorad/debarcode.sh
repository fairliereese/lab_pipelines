docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-debarcode-dbg:v1.0.0 \
 -i /data/ \
 -o /data/debarcoded_reads \
 --name debarcode
