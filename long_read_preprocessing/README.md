# Pre-processing PacBio data

## 0. Make a file with pertinent information

Should be tab-separated and contain the following information:

| Library ID (PB###) | Download link from GHTF | md5sum from GHTF |

Save this in samples.tsv.

## 1. Download data

This automatically checks that the subreads that were downloaded match the md5sum computed by GHTF.

```bash
bash download.sh samples.tsv
mv processing_tables/download.o* processing_tables/processing_output/
mv processing_tables/download.e* processing_tables/processing_output/
```

## 1. CCS

Run CCS to find circular consensus reads.

```bash
bash ccs.sh samples.tsv
```

Check how many reads passed CCS / ensure that each thing finished running.
```bash
cat processing_tables/ccs.o*
mv processing_tables/ccs.o* processing_tables/processing_output/
mv processing_tables/ccs.e* processing_tables/processing_output/
```

## 2. Lima

Run Lima to demultiplex and find full-length reads

```bash
bash lima.sh samples.tsv
```

Check how many reads passed Lima / ensure that each thing finished running.
```bash
cat processing_tables/lima.o*
mv processing_tables/lima.o* processing_tables/processing_output/
mv processing_tables/lima.e* processing_tables/processing_output/
```

## 3. Refine
Find full-length non-chimeric reads

```bash
bash refine.sh samples.tsv
```

Check how many reads passed Refine / ensure that each thing finished running.
```bash
cat processing_tables/refine.o*
mv processing_tables/refine.o* processing_tables/processing_output/
mv processing_tables/refine.e* processing_tables/processing_output/
```
