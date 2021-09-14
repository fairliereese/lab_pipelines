docker run --rm -v $1:/data/ \
 -v /home/freese/mortazavi_lab/ref/genomes/mm10/bwa/:/genome/ \
 bioraddbg/atac-seq-bwa:v1.0.0 \
 -i /data/debarcoded_reads/ \
 -o /data/alignments/ \
 -r /genome/ \
 --name align
