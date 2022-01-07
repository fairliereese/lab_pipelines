
fname=sample_bc1_table.tsv
fastq=/share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/analysis_2k/process/single_cells_barcoded_head.fastq

# loop through the different samples
head single_cells_barcoded_head.fastq | grep ^@ | cut -c18-25
head -40000 single_cells_barcoded_head.fastq | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print substr($0,49,50)}' | head #read1
head -40000 single_cells_barcoded_head.fastq | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print substr($0,49,38) 2 substr($0,88,12)}' | head # Read2

head -40000 single_cells_barcoded_head.fastq | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print substr($0,18,8)}' | head

# from python test
head -40000 single_cells_barcoded_head.fastq | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print substr($0,49,37)}' | head #read1
head -40000 single_cells_barcoded_head.fastq | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print substr($0,49,38) 2 substr($0,88,12)}' # read 2

head -40000 ${sp_fastq} | grep ^@ | awk '{if (substr($0,18,8)=="CTGTCCCG"||substr($0,18,8)=="CTAAGGGA"||substr($0,18,8)=="TTACCTCG"||substr($0,18,8)=="GGGTAGCG") print $0}' | awk '{split($0,a,"_"); print a[4]}' | awk '{split($0,a," "); print a[1]}'| head #read1


d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
a_d=/share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/
opref=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/adrenal_deep/fastq/adrenal_deep
cell_meta=${a_d}analysis_2k/all-well/DGE_unfiltered/cell_metadata.csv
sp_fastq=${a_d}analysis_2k/process/single_cells_barcoded_head.fastq
f1=${a_d}fastq/Adr_2k_retry_R1.fastq
f2=${a_d}fastq/Adr_2k_retry_R2.fastq
python ${d}demux_fastqs.py \
  -o ${opref} \
  -cell_meta ${cell_meta} \
  -sp_fastq ${sp_fastq} \
  -f1 ${f1} \
  -f2 ${f2} \
  -t 4
