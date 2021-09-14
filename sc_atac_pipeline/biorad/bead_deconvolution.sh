docker run --rm -v $1:/data/ \
 bioraddbg/atac-seq-deconvolute:v1.0.0 \
 -i /data/alignments/ \
 -f /data/bead_filtration/ \
 -r mm10 \
 -o /data/deconvoluted_data \
 --name deconv_${2}
