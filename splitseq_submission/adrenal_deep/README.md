```bash
sample='adrenal_deep'
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/deep/
cd ${fastq_dir}
mkdir submission
cd submission
ln -s ../A_*fastq.gz .
```

File name format <tissue>_<age>_<sex>_<replicate>
* A = adrenal
* G = gastroc
* H = heart
* HC = hippocampus
* C = cortex

```bash
# hpc
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=${fastq_dir}/submission/
python ${d}make_submission_spreadsheets.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}${sample} \
  --deep \
  -lib_meta=${meta_dir}/${sample}_metadata.tsv
```

Submit files

```bash
# first submit to test
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_dev.sh ${fastq_dir}${sample}_sr 1

# then submit to prod
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_prod.sh ${fastq_dir}${sample}_sr 1
```
