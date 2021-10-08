```bash
sample='hippocampus_deep'
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/hippocampus/fastq/deep/
cd ${fastq_dir}
mkdir submission
cd submission
ln -s ../HC_*fastq.gz .
```

File name format <tissue>_<age>_<sex>_<replicate>
* A = adrenal
* G = gastroc
* H = hippocampus
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
sbatch ${d}submit_dev.sh ${fastq_dir}${sample}_sr 0 1

# then submit to prod
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_prod.sh ${fastq_dir}${sample}_sr 1 1
```

Add R2 reads
```bash
# hpc
sample='hippocampus_deep'
tissue=hippocampus
depth=deep

fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/hippocampus/fastq/deep
cd ${fastq_dir}
meta_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${sample}
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=${fastq_dir}/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/

cd ${fastq_dir}
ln -s ../HC_*R2*fastq.gz .

python ${d}make_index_reads.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}${sample} \
  --deep \
  -lib_meta=${meta_dir}/${sample}_metadata.tsv

conda activate encode_submissions
eu_register.py -m dev -p file -i ${fastq_dir}/${tissue}_${depth}_sr_r2_file.tsv

eu_register.py -m prod -p file -i ${fastq_dir}/${tissue}_${depth}_sr_r2_file.tsv
```
