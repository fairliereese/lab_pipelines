docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-qc:v1.0.0 \
 -r mm10 \
 -d /data/deconvoluted_data/ \
 -p /data/peaks \
 -o /data/atac_qc \
 --name atac_qc
