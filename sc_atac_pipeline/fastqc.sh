docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-fastqc:v1.0.0 \
 -i /data/ \
 -o /data/fastqc_results \
 --name fastqc
