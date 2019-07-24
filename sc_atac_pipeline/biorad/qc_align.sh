docker run --rm -v $1:/data/ \
 -v /home/freese/mortazavi_lab/ref/mm10/:/genome/ \
 bioraddbg/atac-seq-alignment-qc:v1.0.0 \
 -i /data/alignments/ \
 -r /genome/mm10.fa \
 -o /data/alignment_qc \
 --name alignment_qc
