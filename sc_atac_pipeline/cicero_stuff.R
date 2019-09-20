library(cicero)
library(ggplot2)
library(dplyr)
library(textshape)
library(Seurat)
library(Matrix)
library(stringr)
library(rtracklayer)
library(optparse)

# combined
get_atac_peaks <- function(fname) {
	peaks <- read.table(file=fname, sep=',', header=T, row.names=1)
	return(peaks)
}

reformat_peak_id <- function(peak_id){
	peak_id <- str_replace(peak_id, ':', '_')
	peak_id <- str_replace(peak_id, '-', '_')
	return(peak_id)
}

exprs <- get_atac_peaks('MB1_full_counts.csv')

exprs <- transform(
	exprs, peak_id2=reformat_peak_id(rownames(exprs)))
row.names(exprs) <- exprs$peak_id2
exprs <- subset(exprs, select=-c(peak_id2))
# unique(sapply(exprs, class))

# feature data
peak_names <- rownames(exprs)
feature_data <- data.frame(peak_names)
feature_data <- transform(
	feature_data, gene_short_name=peak_names)
feature_data <- column_to_rownames(as.matrix(feature_data), 'peak_names')

# remove peak id from columns
# exprs <- column_to_rownames()

# sample sheet 
cell_names <- colnames(exprs)
# print(length(cell_names))
pheno_data <- data.frame(cell_names)
pheno_data <- transform(
	pheno_data, cell_type=ifelse(grepl('MB', cell_names, fixed=TRUE), 'MB', 'MT'))
pheno_data <- column_to_rownames(as.matrix(pheno_data), 'cell_names')

# create the CellDataSet object
pd <- AnnotatedDataFrame(pheno_data)
fd <- AnnotatedDataFrame(feature_data)
input_cds <- newCellDataSet(as.matrix(exprs), phenoData=pd, featureData=fd)

# create a cicero CDS using aggregates of similar cells from dim. red. or monocle 
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                      reduction_method = 'tSNE', norm_method = "none")

# get tsne coords
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

# run cicero
chrom_sizes <- read.table('mm10.chrom.sizes')
sample_genome <- subset(chrom_sizes, V1 == "chr18")
sample_genome$V2 <- sapply(sample_genome$V2, as.character)
sample_genome$V2 <- sapply(sample_genome$V2, as.numeric)
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run

# gene activity

# load in your data using rtracklayer
gene_anno <- rtracklayer::readGFF("/data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf")

# rename some columns to match plotting requirements
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#### Add a column for the fData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
pos <- pos[!duplicated(pos$transcript),] # remove all but the first exons per transcript
pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg <- neg[!duplicated(neg$transcript),] # remove all but the first exons per transcript
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,
	c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

# # if you had two datasets to normalize, you would pass both:
# # num_genes should then include all cells from both sets
# unnorm_ga2 <- unnorm_ga
# cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)






