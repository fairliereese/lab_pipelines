docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-macs2:v1.0.0 \
 -i /data/deconvoluted_data/ \
 -r mm10 \
 -o /data/peaks \
 --name peak_calling
