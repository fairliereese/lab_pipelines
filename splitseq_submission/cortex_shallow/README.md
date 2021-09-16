```bash
sample='cortex_shallow'
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/cortex/fastq/shallow/
cd ${fastq_dir}
mkdir submission
cd submission
ln -s ../C_*fastq.gz .
```

File name format <tissue>_<age>_<sex>_<replicate>
* A = adrenal
* G = gastroc
* H = hippocampus
* C = cortex

```bash
# hpc
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=${fastq_dir}/submission/
python ${d}make_submission_spreadsheets.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}${sample} \
  --shallow \
  -lib_meta=${meta_dir}/${sample}_metadata.tsv

# local
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/test/
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=/Users/fairliereese/Documents/programming/mortazavi_lab/bin/lab_pipelines/splitseq_submission/test/
python ${d}make_submission_spreadsheets.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}/test \
  --deep \
  -lib_meta=${meta_dir}/adrenal_metadata.tsv
```

Submit files

```bash
# first submit to test
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_dev.sh ${fastq_dir}/submission/${sample}_sr

# then submit to prod
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_prod.sh ${fastq_dir}/submission/${sample}_sr
```
