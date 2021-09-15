# Pre-processing LR-Split-seq data

## 0. Make a file with pertinent information

Should be tab-separated and contain the following information:

| Library ID (PB###) | Download link from GHTF | md5sum from GHTF | SMRT cell ID |

Save this in samples.tsv.

## 1. Download data

This automatically checks that the subreads that were downloaded match the md5sum computed by GHTF.

```bash
bash download.sh samples.tsv
```

## 1. CCS

Run CCS to find circular consensus reads.

```bash
bash ccs.sh samples.tsv
```

## 2. Lima

Run Lima with the Split-seq adapters.

```bash
bash lima.sh samples.tsv
```

## 3. Refine

Run Refine with the Split-seq adapters and without the poly-A tail requirement.

```bash
bash refine.sh samples.tsv
```
## 4. Concatenate fastqs
```bash
opref=<output directory + prefix>
bash concat_fastqs.sh samples.tsv $opref
```

## 4. Demultiplex Split-seq barcodes
```bash
bash demultiplex.sh $opref
```
