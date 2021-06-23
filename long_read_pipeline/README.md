# Analysis of long-read RNA-seq data


## Data

### Where can I find the data on CRSP?

Long read RNA-seq libraries that have been or are pending sequencing can be seen [on this Google Spreadsheet](https://docs.google.com/spreadsheets/d/1XU0qmc9LQ-40foJ8WSXRPiUC1y9kv0U_g3nAgsuX6VM/edit?usp=sharing). Use this spreadsheet to determine which samples you want to include in your analysis. Human and mouse datasets are stored on separate tabs of the spreadsheet. You will want the `Library ID` column from these spreadsheets.

**Example:** I want to analyze the PGP1 stem cells and PGP1-derived excitatory neurons and astrocytes. I will use the following library IDs:

| Library ID | Description |
| ---------- | ----------- |
| PB103 | PGP1 |
| PB104 | PGP1 |
| PB105 | PGP1-derived excitatory neurons |
| PB106 | PGP1-derived excitatory neurons |
| PB306 | PGP1-derived astrocytes |
| PB307 | PGP1-derived astrocytes |

The data corresponding to each dataset can be found in this directory on HPC3: `/share/crsp/lab/seyedam/share/PACBIO/`. Each dataset is located in a directory corresponding to its library ID.

You will want the FASTQs that have been converted to FLNC (full-length non-chimeric) reads, which are located in either the `Classify` or the `Refine` directory. You are looking for a FASTQ file with a name of the style `*flnc.fastq`. In some cases, this will be directly in this directory. In others, this might require descending into a directory called `*01/`.

**Example:** The data for PB103 is located here: `/share/crsp/lab/seyedam/share/PACBIO/PB103/Classify/D26_isoseq_flnc.fastq`. The data for PB305 is located here: `/share/crsp/lab/seyedam/share/PACBIO/PB305/Refine/B01/flnc.fastq`.

f you're going to use this data, I recommend **symlinking** the FASTQ files that you're going to use in a new directory, which creates a reference to the file but does not copy the underlying data. You can do this with the following command

```bash
mkdir my_analysis_folder/
cd my_analysis_folder
ln -s /share/crsp/lab/seyedam/share/PACBIO/PB103/Classify/D26_isoseq_flnc.fastq pgp1_1.fastq
```

### Where can I find the data on the ENCODE portal?

The data for datasets that have been submitted to ENCODE can be found by searching for the `ENCLB#######` term in the `submitted to DCC` column from the spreadsheet.

If there is a `D##` entry in this column, look at the tab for Human/Mouse datasets and find the corresponding `D##` number, with the ENCODE library ID in the `ENCODE library accession` column. This is more common with older datasets.

The datasets that I'm interested in have the following ENCODE library IDs.

| Library ID | ENCODE Library ID |
| ---------- | ----------------- |
| PB103 | ENCLB266AQD |
| PB104 | ENCLB644KMO |
| PB105 | ENCLB389UCF |
| PB106 | ENCLB532XRH |
| PB306 | ENCLB652QMD |
| PB307 | ENCLB486CFV |

If you don't wish to use the data already on CRSP, you can download the corresponding FASTQ files from the portal using these ENCODE library IDs.

**Example:** If I want to download the FASTQ for PB103, I will search `ENCLB266AQD` on the ENCODE portal with [this result](https://www.encodeproject.org/search/?searchTerm=ENCLB266AQD). The [resulting experiment](https://www.encodeproject.org/experiments/ENCSR676IWT/) contains the following FASTQ files: `ENCFF251CBB`and `ENCFF954UFG`, where the latter corresponds to PB103/ENCLB266AQD. If you right click on the download symbol by the FASTQ name, an option to "Copy Link Address" will pop up. Click it and then use the following strategy to download the FASTQ for the dataset you're interested in.

```bash
wget https://www.encodeproject.org/files/ENCFF954UFG/@@download/ENCFF954UFG.fastq.gz .
```

### Rename the files and make a samples file

After you've downloaded or moved all the data to the spot you'd like to analyze it, rename the files with more human-readable names.

```bash
mv ENCFF954UFG.fastq > pgp1_1.fastq
mv ENCFF251CBB.fastq > pgp1_2.fastq
mv ENCFF919JFJ.fastq > excite_neuron_1.fastq
mv ENCFF982WKN.fastq > excite_neuron_2.fastq
mv ENCFF316EZQ.fastq > astro_1.fastq
mv ENCFF474GEK.fastq > astro_2.fastq
```

Make a samples file which will make it easy to iteratively process each dataset.

```bash
touch samples.txt
printf "pgp1_1\n" >> samples.txt
printf "pgp1_2\n" >> samples.txt
printf "excite_neuron_1\n" >> samples.txt
printf "excite_neuron_2\n" >> samples.txt
printf "astro_1\n" >> samples.txt
printf "astro_2\n" >> samples.txt
```

## Map reads to the genome

Long-read RNA-seq data requires specialized mappers able to deal with the long reads and the high error rate of long reads. We use minimap2, which is already installed on HPC3.

Use the following commands to run minimap2 on every sample
```bash
n=`wc -l samples.txt | cut -d' ' -f1`
sbatch --array=1-${n} map_reads.sh samples.txt <genome>
```

The code in `map_reads.sh` is also reproduced below.
```bash
module load samtools
module load minimap2/2.17

# input arguments
samples=$1
genome=$2

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

# run minimap on this sample
fastq=${p}.fastq
sam=${p}_mapped.sam
log=${p}_minimap.log
minimap2 \
    -t 16 \
    -ax splice:hq \
    -uf \
    --MD \
    --secondary=no \
    $genome \
    $fastq > \
    $sam 2> \
    $log
```

Count the number of reads that were mapped

```bash
n_total=0
while read sample
do
  f=${sample}_mapped.sam
  n_curr=`samtools view -c $f`
  n_total=$((n_total + n_curr))
done < samples.txt
echo "$n_total reads after Minimap2"
```

## TranscriptClean
Correct common long-read sequencing artifacts. Meant to be run on the cluster.

Download TranscriptClean from here: https://github.com/mortazavilab/TranscriptClean, and use the path where you downloaded it to as `tc_path` in the following block.

```bash
n=`wc -l samples.txt | cut -d' ' -f1`
sbatch --array=1-${n} tc.sh samples.txt <genome> <tc_path>
```

The code in `tc.sh` is also reproduced below.
```bash
module load samtools

# input arguments
samples=$1
genome=$2
tc_path=$3

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

# input and output filenames
sam=${p}_mapped.sam
bam=${p}_mapped.bam
sort_bam=${p}_sorted.bam
sort_sam=${p}_sorted.sam

# first sort the sam file
samtools view -Sb $sam > $bam
samtools sort $bam > $sort_bam
samtools view -h $sort_bam > $sort_sam

# run TranscriptClean
mkdir -p ${p}_tmp/
python ${tc_path}/TranscriptClean.py \
   --sam $sort_sam \
   --genome $genome \
   -t 16 \
   --canonOnly \
   --tmpDir ${p}_tmp \
   --deleteTmp \
   --outprefix $p
```
rm -r ${p}_tmp/

Count the number of reads that passed TranscriptClean

```bash
n_total=0
while read sample
do
  f=${sample}_clean.sam
  n_curr=`samtools view -c $f`
  n_total=$((n_total + n_curr))
done < samples.txt
echo "$n_total reads after TranscriptClean"
```

## TALON

TALON is designed to annotate long reads to their transcripts of origin, as well as identify novel transcripts.

Download and install TALON according to the instructions on the TALON repository: https://github.com/mortazavilab/TALON

### TALON label reads

Before running TALON, we have to determine what the nucleotide composition of the end of each read is, which will help us filter for internal priming.

```bash
n=`wc -l samples.txt | cut -d' ' -f1`
sbatch --array=1-${n} talon_label_reads.sh samples.txt <genome>
```

The code in `talon_label_reads.sh` is also reproduced below.
```bash
# input arguments
samples=$1
genome=$2

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

sam=${p}_clean.sam
talon_label_reads \
    --f $sam \
    --g $genome \
    --t 16 \
    --ar 20  \
    --deleteTmp  \
    --o $p
```

### Create a TALON config file

Create a comma-separated file that provides the sample name, sample description, platform, and location of each input sam file. There's not a great way to automate this but the code for my example is shown below. It uses the output `*_labeled.sam` files from the previous steps.

```bash
touch talon_config.csv
printf "pgp1_1,pgp1,SequelI,pgp1_1_labeled.sam\n" >> talon_config.csv
printf "pgp1_2,pgp1,SequelI,pgp1_2_labeled.sam\n" >> talon_config.csv
printf "excite_neuron_1,excitatory_neuron,SequelI,excite_neuron_1_labeled.sam\n" >> talon_config.csv
printf "excite_neuron_2,excitatory_neuron,SequelI,excite_neuron_2_labeled.sam\n" >> talon_config.csv
printf "astro_1,astrocyte,SequelII,astro_1_labeled.sam\n" >> talon_config.csv
printf "astro_2,astrocyte,SequelII,astro_2_labeled.sam\n" >> talon_config.csv
```

### Initialize TALON db and run TALON

```bash
sbatch talon.sh talon_config.csv <annotation.gtf> <annotation_name> <genome_name> <output file prefix>
```

The code in `talon.sh` is also reproduced below.
```bash
config=$1

annot=$2
annot_name=$3

genome_name=$4

opref=$5

talon_initialize_database \
    --f $annot \
    --g $genome_name \
    --a $annot_name \
    --l 0 \
    --idprefix ENCODEH \
    --5p 500 \
    --3p 300 \
    --o ${opref}

talon \
    --f $config \
    --db ${opref}.db \
    --build $genome_name \
    --t 64 \
    --o ${opref}
```

Count the number of reads that passed TALON

```bash
n_total=`wc -l <talon.db>`
echo "$n_total reads after TALON"
```

### Filter TALON transcripts

Filter novel transcripts for internal priming and for reproducibility. Look at TALON documentation for specifics on arguments.

```bash
db=pgp1.db
talon_filter_transcripts \
    --db $db \
    -a <annot_name> \
    --maxFracA=0.5 \
    --minCount=5 \
    --minDatasets=2 \
    --o pgp1_pass_list.csv
```

### Create an unfiltered and a filtered abundance file

#### Unfiltered
```bash
db=pgp1.db
talon_abundance \
    --db $db \
    -a <annot_name> \
    -b <genome_name> \
    --o pgp1
```

#### Filtered
```bash
db=pgp1.db
talon_abundance \
    --db $db \
    -a <annot_name> \
    -b <genome_name> \
    --whitelist pgp1_pass_list.csv \
    --o pgp1
```
