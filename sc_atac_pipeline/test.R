library(cicero)
library(ggplot2)
library(dplyr)
library(textshape)
library(Seurat)
library(Matrix)
library(stringr)

# combined
get_atac_peaks <- function(fname) {

	# peaks <- Matrix(
	# 	as.matrix(
	# 		read.table(
	# 			file=fname,
	# 			sep=',',
	# 			header=T, 
	# 			row.names=1
	# 		), 
	# 	sparse=T
	# ), sparse=T)
	peaks <- read.table(file=fname, sep=',', header=T, row.names=1)

	return(peaks)
}

reformat_peak_id <- function(peak_id){
	peak_id <- str_replace(peak_id, ':', '_')
	peak_id <- str_replace(peak_id, '-', '_')
	return(peak_id)
}

# # load all data
# mb1_atac <- get_atac_peaks('MB1_full_counts.csv')
# mb2_atac <- get_atac_peaks('MB2_full_counts.csv')
# mt1_atac <- get_atac_peaks('MT1_full_counts.csv')
# mt2_atac <- get_atac_peaks('MT2_full_counts.csv')

# # "expression" matrix
# exprs <- Merge(mb1_atac, mb2_atac, by='peak_id')
# exprs <- Merge(exprs, mt1_atac, by='peak_id')
# exprs <- Merge(exprs, mt2_atac, by='peak_id')

# testing
mb1_atac <- get_atac_peaks('MB1_full_counts.csv')
exprs <- mb1_atac

exprs <- transform(
	exprs, peak_id2=reformat_peak_id(rownames(exprs)))
exprs <- column_to_rownames(as.matrix(exprs), 'peak_id2')

# feature data
peak_names <- rownames(exprs)
feature_data <- data.frame(peak_names)
feature_data <- transform(
	feature_data, gene_short_name=peak_names)
feature_data <- column_to_rownames(as.matrix(feature_data), 'peak_names')

# remove peak id from columns
# exprs <- column_to_rownames()
# exprs <- subset(exprs, select=-c(peak_id2))

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
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run
head(conns)






