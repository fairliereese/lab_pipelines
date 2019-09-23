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

main <- function() {

	load_packages()
	args <- parse_arguments()

	# load in and reformat input counts mat to be compatible with cicero
	print('Reading in counts matrix...')
	exprs <- get_atac_peaks(args$exprs)
	exprs <- transform(
		exprs, peak_id2=reformat_peak_id(rownames(exprs)))
	row.names(exprs) <- exprs$peak_id2
	exprs <- subset(exprs, select=-c(peak_id2))

	# feature data for cds object
	peak_names <- rownames(exprs)
	feature_data <- data.frame(peak_names)
	feature_data <- transform(
		feature_data, gene_short_name=peak_names)
	feature_data <- column_to_rownames(as.matrix(feature_data), 'peak_names')

	# cell data for cds object
	cell_names <- colnames(exprs)
	pheno_data <- data.frame(cell_names)
	pheno_data <- transform(
		pheno_data, cell_type=ifelse(grepl('MB', cell_names, fixed=TRUE), 'MB', 'MT'))
	pheno_data <- column_to_rownames(as.matrix(pheno_data), 'cell_names')

	# create the cds (CellDataSet) object
	print('Creating CellDataSet object...')
	pd <- AnnotatedDataFrame(pheno_data)
	fd <- AnnotatedDataFrame(feature_data)
	input_cds <- newCellDataSet(as.matrix(exprs), phenoData=pd, featureData=fd)

	# create a cicero CDS TSNE coords
	print('Creating Cicero CellDataSet using tSNE')
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
	print('Running Cicero...')
	genome <- read.table(args$genome_sizes, header=T) 
	genome$V2 <- sapply(genome$V2, as.character)
	genome$V2 <- sapply(genome$V2, as.numeric)
	conns <- run_cicero(cicero_cds, genome) 

	# load in annotation data using rtracklayer
	gene_anno <- rtracklayer::readGFF(args$gtffile)

	# rename some columns to match plotting requirements
	gene_anno$chromosome <- gene_anno$seqid
	gene_anno$gene <- gene_anno$gene_id
	gene_anno$transcript <- gene_anno$transcript_id
	gene_anno$symbol <- gene_anno$gene_name

	# determine which peaks are at gene promoters
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


	# save files
	activity_fname <- paste0(args$odir, '/', args$sample_name, '_cicero_activity.rds')
	saveRDS(unnorm_ga, activity_fname)

	num_gene_fname <- paste0(args$odir, '/', args$sample_name, '_cicero_num_genes.rds')
	saveRDS(num_genes, num_gene_fname)
}

get_atac_peaks <- function(fname) {
	peaks <- read.table(file=fname, sep=',', header=T, row.names=1)
	return(peaks)
}

reformat_peak_id <- function(peak_id){
	peak_id <- str_replace(peak_id, ':', '_')
	peak_id <- str_replace(peak_id, '-', '_')
	return(peak_id)
}

load_packages <- function() {
	suppressPackageStartupMessages(library(cicero))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(dplyr))
	suppressPackageStartupMessages(library(textshape))
	suppressPackageStartupMessages(library(Seurat))
	suppressPackageStartupMessages(library(Matrix))
	suppressPackageStartupMessages(library(stringr))
	suppressPackageStartupMessages(library(rtracklayer))
	suppressPackageStartupMessages(library(optparse))
	return
}

parse_arguments <- function() {

	option_list <- list(
		make_option('--counts_mat', action='store', dest='exprs',
					help='counts matrix for one sample (CSV)'),
		make_option('--genome_sizes', action='store', dest='genome_sizes',
					help='tab-delimited genome size file for organism'),
		make_option('--gtf', action='store', dest='gtffile', 
					help='GENCODE gtf annotation for organism'),
		make_option('--o', action='store', dest='odir',
					help='output directory to store gene activity matrix in'),
		make_option('--sample', action='store', dest='sample_name',
					help='sample name, to be used to save output')

	)

	opt <- parse_args(OptionParser(option_list=option_list))
	return(opt)
}

main()


