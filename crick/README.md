# Pilot analysis of 8 cell lines

* Data processed according to the code in this notebook can be found [here](http://crick.bio.uci.edu/freese/rnawg/talon/)
* GENCODE v29 + SIRV set 3 reference annotation can be found [here](http://crick.bio.uci.edu/freese/rnawg/ref/)

We chose the following cell lines because they have data from many other assays on the portal. 

| Cell line | Experiment ID |
| ---------- | ----------- |
| HCT116  | ENCSR388RWE |
| PC-3    | ENCSR947JHA |
| Panc1   | ENCSR181CES |
| K562    | ENCSR526TQU |
| GM12878 | ENCSR962BVU |
| HepG2   | ENCSR634AKC |
| MCF-7   | ENCSR810XLL |
| IMR90   | ENCSR113MEY |

I compiled these experiments into a [cart on the ENCODE portal](https://www.encodeproject.org/carts/8487de1a-bb20-4ff8-a465-f514ec58035f/).

## Data download 

I will be using the filtered, mapped bams from the portal because these have been run through [TranscriptClean](https://github.com/mortazavilab/TranscriptClean) with the [human reference variants](https://www.encodeproject.org/files/ENCFF911UGW/). The links to these files are in `files.txt`. Download them using the following:

```bash
xargs -L 1 curl -O -J -L < files.txt
```

From the `metadata.tsv` file, rename the bam files to have more human-readable names.
```bash
mv ENCFF400LRT.bam hct116.bam
mv ENCFF044LIA.bam pc3.bam
mv ENCFF745DHX.bam panc1.bam
mv ENCFF322UJU.bam k562.bam
mv ENCFF219UJG.bam gm12878.bam
mv ENCFF814ABW.bam hepg2.bam
mv ENCFF118JEI.bam mcf7.bam
mv ENCFF545PJV.bam imr90.bam
```

Finally, make a `samples.txt` file to control compute cluster task array.
```bash
touch samples.txt
printf "hct116\n" >> samples.txt
printf "pc3\n" >> samples.txt
printf "panc1\n" >> samples.txt
printf "k562\n" >> samples.txt
printf "gm12878\n" >> samples.txt
printf "hepg2\n" >> samples.txt
printf "mcf7\n" >> samples.txt
printf "imr90\n" >> samples.txt
```

## TALON

TALON is designed to annotate long reads to their transcripts of origin, as well as identify novel transcripts.

Download and install TALON according to the instructions on the TALON repository: https://github.com/mortazavilab/TALON

### TALON label reads

Before running TALON, we have to determine what the nucleotide composition of the end of each read is, which will help us filter for internal priming.

```bash
genome=/dfs6/pub/freese/mortazavi_lab/ref/hg38/hg38_sirv3.fasta
n=`wc -l samples.txt | cut -d' ' -f1`
sbatch --array=1-${n} talon_label_reads.sh samples.txt $genome
```

<!-- 7/8/21 - Rerunning sample 5 and 8 (GM12878 and IMR90 because they failed) -->

### Create a TALON config file

Create a comma-separated file that provides the sample name, sample description, platform, and location of each input sam file. There's not a great way to automate this but the code for my example is shown below. It uses the output `*_labeled.sam` files from the previous steps.

```bash
touch talon_config.csv
printf "hct116,hct116,SequelI,hct116_labeled.sam\n" >> talon_config.csv
printf "pc3,pc3,SequelI,pc3_labeled.sam\n" >> talon_config.csv
printf "panc1,panc1,SequelI,panc1_labeled.sam\n" >> talon_config.csv
printf "k562,k562,SequelI,k562_labeled.sam\n" >> talon_config.csv
printf "gm12878,gm12878,SequelI,gm12878_labeled.sam\n" >> talon_config.csv
printf "hepg2,hepg2,SequelI,hepg2_labeled.sam\n" >> talon_config.csv
printf "mcf7,mcf7,SequelI,mcf7_labeled.sam\n" >> talon_config.csv
printf "imr90,imr90,SequelI,imr90_labeled.sam\n" >> talon_config.csv
```

### Initialize TALON db and run TALON

```bash
mkdir talon
annot=~/mortazavi_lab/ref/gencode.v29/gencode_v29_sirv3.gtf
sbatch talon.sh talon_config.csv $annot gencode_v29 hg38 talon/pilot
```

### Filter TALON transcripts

Filter novel transcripts for internal priming and for reproducibility. Look at the TALON documentation for specifics on arguments. I'll use the default filtering parameters that were used for the TALON paper.

```bash
db=talon/pilot.db
talon_filter_transcripts \
    --db $db \
    -a gencode_v29 \
    --maxFracA=0.5 \
    --minCount=5 \
    --minDatasets=2 \
    --o talon/pilot_pass_list.csv
```

### Create filtered and unfiltered abundance files

#### Unfiltered - used for gene-level quantification

* Link: [`pilot_talon_abundance.tsv`](http://crick.bio.uci.edu/freese/rnawg/talon/pilot_talon_abundance.tsv)

```bash
db=talon/pilot.db
talon_abundance \
    --db $db \
    -a gencode_v29 \
    -b hg38 \
    --o talon/pilot
```

#### Filtered - used for transcript-level quantification

* Link: [`pilot_talon_abundance_filtered.tsv`](http://crick.bio.uci.edu/freese/rnawg/talon/pilot_talon_abundance_filtered.tsv)

```bash
db=talon/pilot.db
talon_abundance \
    --db $db \
    -a gencode_v29 \
    -b hg38 \
    --whitelist talon/pilot_pass_list.csv \
    --o talon/pilot
```

### Create a filtered and unfiltered GTFs

#### Unfiltered - contains the entire set of GENCODE and SIRV Set 3 transcripts and all novel transcripts found by TALON

* Link: [`pilot_talon.gtf`](http://crick.bio.uci.edu/freese/rnawg/talon/pilot_talon.gtf)

```bash
db=talon/pilot.db
talon_create_GTF \
    --db $db \
    -a gencode_v29 \
    -b hg38 \
    --o talon/pilot
```

#### Filtered - contains all known transcripts and all novel transcripts that pass filtering

* Link: [`pilot_talon_observedOnly.gtf`](http://crick.bio.uci.edu/freese/rnawg/talon/pilot_talon_observedOnly.gtf)


```bash
db=talon/pilot.db
talon_create_GTF \
    --db $db \
    -a gencode_v29 \
    -b hg38 \
    --whitelist talon/pilot_pass_list.csv \
    --observed \
    --o talon/pilot
```

## Swan

For the purposes of this analysis, I used the [development branch of Swan](https://github.com/mortazavilab/swan_vis/tree/development/swan_vis). The code I used to run Swan is in `swan/swan.ipynb`.

```bash
mkdir swan
```

Output from Swan includes:
* [Saved SwanGraph that includes only non-ISM filtered transcripts](http://crick.bio.uci.edu/freese/rnawg/swan/swan.p)
* [Edge counts matrix from non-ISM filtered transcripts](http://crick.bio.uci.edu/freese/rnawg/swan/pilot_edge_abundance.tsv)
* [Saved SwanGraph that includes all transcripts](http://crick.bio.uci.edu/freese/rnawg/swan/swan_unfiltered.p)
* [Edge counts matrix from all transcripts](http://crick.bio.uci.edu/freese/rnawg/swan/pilot_unfiltered_edge_abundance.tsv)
