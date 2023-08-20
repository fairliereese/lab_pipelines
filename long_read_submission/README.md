# Submission of PacBio data

## 0. Make a file with pertinent information

Should be tab-separated and contain the following information.

One library per rep:
| Library ID | Nickname |
| ---------- | -------- |
| PB352 | PB352 |

Multiple fastqs per rep:


| Library IDs, comma-separated | Nickname |
PB260,PB267 | lrgasp_es_1 |

Save this in samples.tsv.

## 1. Prep files for submission

Fastq files for each SMRT cell in the input library IDs are concatenated and gzipped into `files/nickname/raw_fastqs/nickname.fastq.gz`.

Subreads are for each SMRT cell in the input library IDs are symlinked in `files/nickname/subreads/` with their original subreads names.

```bash
bash prep_files.sh samples.tsv
```

# gzip a bunch of fastqs
```bash
snakemake \
  -s Snakefile \
  -j 1 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
  ```
