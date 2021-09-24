```bash
sample='heart_deep'
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/heart/fastq/deep/
cd ${fastq_dir}
mkdir submission
cd submission
ln -s ../H_*fastq.gz .
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
local
```bash
tissue=heart
depth=deep
data_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/${tissue}/fastq/${depth}/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
script=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/patch_frag_size.py
scp freese@hpc3.rcic.uci.edu:${data_dir}${tissue}_${depth}_sr_library.tsv ${bin_dir}
cp ${bin_dir}${tissue}_${depth}_sr_library.tsv ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv
python ${script} ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv ${tissue}_${depth}
scp ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv freese@hpc3.rcic.uci.edu:~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
```

```bash
tissue=heart
depth=deep
data_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/${tissue}/fastq/${depth}/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
conda activate encode_submissions
eu_register.py -m prod -p library -i ${bin_dir}/${tissue}_${depth}_sr_library_patch.tsv --patch
```
