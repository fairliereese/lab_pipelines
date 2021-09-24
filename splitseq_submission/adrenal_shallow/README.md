```bash
cd /share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/
mkdir submission
cd submission
ln -s ../A_*fastq.gz .
```

File name format <tissue>_<age>_<sex>_<replicate>
* A = adrenal
* G = gastroc
* H = hippocampus
* C = cortex

```bash
# hpc
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/submission/
python ${d}make_submission_spreadsheets.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}adrenal \
  -lib_meta=${d}adrenal_metadata.tsv

# local
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
fastq_dir=/Users/fairliereese/Documents/programming/mortazavi_lab/bin/lab_pipelines/splitseq_submission/test/
python ${d}make_submission_spreadsheets.py \
  -d ${fastq_dir} \
  -o ${fastq_dir}adrenal \
  -lib_meta=${d}adrenal_metadata.tsv
```

Submit files

```bash
# first submit to test
d=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/
sbatch ${d}submit_dev.sh /share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/submission/adrenal_sr

# then submit to prod
sbatch ${d}submit_prod.sh /share/crsp/lab/seyedam/share/Heidi_Liz/adrenal/fastq/submission/adrenal_sr
```

Patch fragment size
local
```bash
tissue=adrenal
depth=shallow
data_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/${tissue}/fastq/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
script=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/patch_frag_size.py
scp freese@hpc3.rcic.uci.edu:${data_dir}${tissue}_sr_library.tsv ${bin_dir}
cp ${bin_dir}${tissue}_sr_library.tsv ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv
python ${script} ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv ${tissue}_${depth}
scp ${bin_dir}${tissue}_${depth}_sr_library_patch.tsv freese@hpc3.rcic.uci.edu:~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
```

```bash
tissue=adrenal
depth=shallow
data_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/${tissue}/fastq/${depth}/submission/
bin_dir=~/mortazavi_lab/bin/lab_pipelines/splitseq_submission/${tissue}_${depth}/
conda activate encode_submissions
eu_register.py -m prod -p library -i ${bin_dir}/${tissue}_${depth}_sr_library_patch.tsv --patch
```
