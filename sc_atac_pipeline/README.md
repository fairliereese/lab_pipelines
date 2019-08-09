# scATAC analysis pipeline

1. Run BioRad Docker pipeline on each sample
```
bash driver.sh MB1/ biorad/
bash driver.sh MB2/ biorad/
bash driver.sh MT1/ biorad/
bash driver.sh MT2/ biorad/
```

2. Merge peaks from each sample. Make a csv file called peakfiles.csv with the path to each sample's narrowPeak file output from the Docker pipeline that looks like this
```
MB1/peaks/alignments.possorted.tagged_peaks.narrowPeak,MT1/peaks/alignments.possorted.tagged_peaks.narrowPeak,MB2/lanes1_2/peaks/alignments.possorted.tagged_peaks.narrowPeak,MT2/peaks/alignments.possorted.tagged_peaks.narrowPeak
```
Then run the merge peaks script to generate a set of peaks common to each sample. 
```
python merge_peaks.py -p peakfiles.csv --o MT1_MT2_MB1_MB2
```

3. Generate a counts matrix using the combined peaks file, barcode list, and the mapped reads for each sample.
```
# MB1
python gen_count_matrix.py -p MT1_MT2_MB1_MB2_merged.bed \
  -b  MB1/cells_filtered/alignments.possorted.tagged.final.bam \
  -bc MB1/cells_filtered/aboveKneeBarcodes.csv \
  -e MB1 \
  --o MB1/ \
  --full_csv

# MB2
python gen_count_matrix.py -p MT1_MT2_MB1_MB2_merged.bed \
  -b  MB2/cells_filtered/alignments.possorted.tagged.final.bam \
  -bc MB2/cells_filtered/aboveKneeBarcodes.csv \
  -e MB2 \
  --o MB2/ \
  --full_csv

# MT1
python gen_count_matrix.py -p MT1_MT2_MB1_MB2_merged.bed \
  -b  MT1/cells_filtered/alignments.possorted.tagged.final.bam \
  -bc MT1/cells_filtered/aboveKneeBarcodes.csv \
  -e MT1 \
  --o MT1/ \
  --full_csv

# MT2
python gen_count_matrix.py -p MT1_MT2_MB1_MB2_merged.bed \
  -b  MT2/cells_filtered/alignments.possorted.tagged.final.bam \
  -bc MT2/cells_filtered/aboveKneeBarcodes.csv \
  -e MT2 \
  --o MT2/ \
  --full_csv
```

4. If you have specific genes you're interested in investigating, you can run check_peak.py with input chromosomal coordinates or gene names to see what the coverage looks like for each scATAC sample. For instance, consider myoblast and myotube marker genes (mt_mb_genes):
```
Myog,Acta1,Ckm,Mybph,Myh3,TTN,label=Myotube
Myod1,Des,Myf5,Pax7,label=Myoblast
```
Each gene or set of coordinates will be analyzed for coverage separately, and each line of gene or coordinates will be analyzed together so you can combine groups of genes that may indicate one cell type label vs another. 

Now, you can run check_peak.
```
# MB1
python check_peak.py \ 
  -gc mt_mb_genes \
  -m MB1/MB1_full_counts.csv

# MB2 
python check_peak.py \
  -gc mt_mb_genes \
  -m MB2/MB2_full_counts.csv

# MT1 
python check_peak.py \
 -gc mt_mb_genes \
 -m MT1/MT1_full_counts.csv

# MT2 
python check_peak.py \
  -gc mt_mb_genes \
  -m MT2/MT2_full_counts.csv
```

5. Coming soon! Make upset plots 

6. Coming soon! Seurat analysis


