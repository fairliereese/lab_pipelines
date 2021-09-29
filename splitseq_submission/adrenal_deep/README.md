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

Patch fragment size

```bash
conda activate encode_submissions
eu_register.py -m prod -p library -i ~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/adrenal_deep/adrenal_deep_sr_library_patch.tsv --patch

```

Patch submitted files due to sample swaps
```bash
tissue=adrenal
depth=deep
data_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/${tissue}/fastq/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
script=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/patch_sample_swaps.py
cp ${data_dir}${tissue}_sr_file.tsv ${bin_dir}${tissue}_${depth}_sr_file_patch.tsv
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/submission/
python ${script} \
  -d ${fastq_dir} \
  -o ${bin_dir}${tissue}_${depth} \
  --deep

conda activate encode_submissions
eu_register.py -m dev -p file -i ${bin_dir}/${tissue}_${depth}_sr_file_patch_1.tsv --patch -w
eu_register.py -m dev -p file -i ${bin_dir}/${tissue}_${depth}_sr_file_patch_2.tsv --patch -w

eu_register.py -m prod -p file -i ${bin_dir}/${tissue}_${depth}_sr_file_patch_1.tsv --patch -w
eu_register.py -m prod -p file -i ${bin_dir}/${tissue}_${depth}_sr_file_patch_2.tsv --patch -w
```
