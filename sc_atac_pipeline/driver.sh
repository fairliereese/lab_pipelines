#!/bin/bash

set -e

dpath=$1
bpath=$2

# fastqc
bash ${bpath}fastqc.sh $dpath

# debarcoding
bash ${bpath}debarcode.sh $dpath

# alignment
bash ${bpath}align.sh $dpath

# alignment qc
bash ${bpath}qc_align.sh $dpath

# bead filtration
bash ${bpath}bead_filtration.sh $dpath

# bead deconvolution
bash ${bpath}bead_deconvolution.sh $dpath mb2_12

# cell filtration
bash ${bpath}cell_filtration.sh $dpath 

# peak calling
bash ${bpath}peak_calling.sh $dpath

# atac qc
bash ${bpath}atac_qc.sh $dpath

# # counts matrix (docker image currently doesn't work!)
# bash ${bpath}count_mat.sh $dpath

# report generation 
bash ${bpath}gen_report.sh $dpath


