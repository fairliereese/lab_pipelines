File name format <genotype>_<tissue>_<replicate>
* B - Bl6/Cast
* X - 5xFAD/Cast
* A = adrenal
* G = gastroc
* HC = hippocampus
* H = heart
* C = cortex

```bash
fastq_dir=/share/crsp/lab/seyedam/share/Heidi_Liz/bl6_5x_cast_ctx_hipp/fastq/shallow/
lib_meta=8mo_1_shallow_metadata.tsv
lib_type=shallow
opref=8mo_1_shallow

python ../make_submission_spreadsheets_8mo.py \
  -d $fastq \
  -o $opref \
  -lib_meta $lib_meta \
  -lib_type $lib_type
```
