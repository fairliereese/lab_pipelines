docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-report:v1.0.0 \
 -i /data/ \
 -o /data/report
