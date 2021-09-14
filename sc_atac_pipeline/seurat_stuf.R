library(ggplot2)
library(dplyr)
# library(tidyverse)
library(Seurat)
library(Matrix)


peaks <- Matrix(
	as.matrix(
		read.table(
			file="MB1_full_counts.csv",
			sep=',',
			header=T, 
			row.names=1
		), 
	sparse=T
), sparse=T)

act_matrix <- 
CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "/data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation_seurat.gtf", 
    seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)
# activity.matrix <- 
# CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "/data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf", 
#     seq.levels = c(1:21, "X", "Y"), upstream = 2000, verbose = TRUE)

atac <- CreateSeuratObject(counts=peaks, assay='ATAC', project='test')
atac[["ACTIVITY"]] <- CreateAssayObject(counts = act_matrix)


# make violin plot to look at qualities
pdf("figures/qc_violin_plot.pdf", width=10, height=10)
VlnPlot(atac, features = c('nCount_ATAC', 'nFeature_ATAC'), ncol=2, pt.size=1)
VlnPlot(atac, features = c('nCount_ACTIVITY', 'nFeature_ACTIVITY'), ncol=2, pt.size=1)
dev.off()

# keep cells with above 5000 reads
atac <- subset(atac, subset = nCount_ATAC > 5000)
atac$tech <- "atac"

# make umap
DefaultAssay(atac) <- "ATAC"
VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 100))
atac <- RunLSI(atac, n = 50, scale.max = NULL)
atac <- RunUMAP(atac, reduction = "lsi", dims = 1:50)

pdf("figures/umap.pdf", width=10, height=10)
p1 <- DimPlot(atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
print(p1)
dev.off()

pdf("figures/feature_plot.pdf", width=10, height=10)
p1 <- FeaturePlot(atac, features=c("nCount_ATAC"), reduction="umap")
print(p1)
dev.off()

# myog umap feature plot
pdf("figures/myog_plot", width=10, height=10)
p1 <- FeaturePlot(atac, features=c("MYOG"), reduction="umap")
print(p1)
dev.off()


# clusters maybe
DefaultAssay(atac) <- "ACTIVITY"
atac <- RunPCA(atac, dims=1:100)
atac <- FindNeighbors(atac, reduction='pca')
atac <- FindClusters(atac, resolution=1.0)

# lil shitty fake code

# save sparse matrix
saveRDS(peaks,  file="peaks.rds")
peaks <- readRDS("peaks.rds")

# save seurat object
saveRDS(atac,  file="atac_seurat.rds")
atac <- readRDS("atac_seurat.rds")

# save rda seurat obj
save(atac, file='seurat_obj.rda')

# generate feature plot with gene
gen_feat_plot_gene <- function(obj, gene) {
	fname = paste0('figures/',gene,'_plot.pdf')
	print(fname)
	pdf(fname, width=10, height=10)
	p1 <- FeaturePlot(obj, features=c(gene), reduction="umap")
	print(p1)
	dev.off()
}

gen_feat_plot_gene(atac, 'Myog')

# combined
gen_atac_seurat_obj <- function(fname) {

peaks <- Matrix(
	as.matrix(
		read.table(
			file=fname,
			sep=',',
			header=T, 
			row.names=1
		), 
	sparse=T
), sparse=T)

act_matrix <- 
CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "/data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation_seurat.gtf", 
    seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)

atac <- CreateSeuratObject(counts=peaks, assay='ATAC', project='test')
atac[["ACTIVITY"]] <- CreateAssayObject(counts = act_matrix)

return(atac)
}

merge_seurat_things <- function(seurats) {

	merged <- merge(seurats[1], y=seurats[2:length(seurats)], project="MB_MT")
	return(merged)
}


qc_violin_plot <- function(seurat_obj) {
	# make violin plot to look at qualities
	pdf("figures/qc_violin_plot.pdf", width=10, height=10)
	p1 <- VlnPlot(seurat_obj, features = c('nCount_ATAC', 'nFeature_ATAC'), ncol=2, pt.size=1)
	p2 <- VlnPlot(seurat_obj, features = c('nCount_ACTIVITY', 'nFeature_ACTIVITY'), ncol=2, pt.size=1)
	print(p1)
	print(p2)
	dev.off()
}

mb1_atac <- gen_atac_seurat_obj('MB1_full_counts.csv')
mb2_atac <- gen_atac_seurat_obj('MB2_full_counts.csv')
mt1_atac <- gen_atac_seurat_obj('MT1_full_counts.csv')
mt2_atac <- gen_atac_seurat_obj('MT2_full_counts.csv')

# seurats <- c(mb1_atac, mb2_atac, mt1_atac, mt2_atac)
# mb_mt <- merge_seurat_things(seurats)

mb_mt <- merge(mb1_atac, y=c(mb2_atac, mt1_atac, mt2_atac), project='MB_MT')

save(mb_mt, file='combined_atac.rda')

# generate violin plot to assess cutoff range
qc_violin_plot(mb_mt)

# subset based on output qc violin plot
mb_mt <- subset(mb_mt, subset = nCount_ATAC > 1000)

# analyze raw ATAC 
# need to set default assay to ATAC as opposed to activity
DefaultAssay(mb_mt) <- 'ATAC'
mb_mt <- mb_mt %>%
	FindVariableFeatures() %>%
	NormalizeData() %>%
	ScaleData() %>%
	RunLSI(dims=1:100) %>%
	RunTSNE(reduction="lsi", dims=1:50) %>%
	RunUMAP(reduction="lsi", dims=1:50) %>% 
	FindNeighbors(reduction="lsi", dims=1:50)

# find clusters
mb_mt <- FindClusters(mb_mt, resolution=1)

# add sample IDs to seurat obj
mb_mt@meta.data$SampleID <- do.call("rbind", strsplit(colnames(mb_mt), "_"))[,1]
mb_mt@meta.data$CellType <- substr(mb_mt@meta.data$SampleID, 1,2)

pdf("figures/atac_cool_stuff.pdf", width=10, height=10)
# umap colored by clusters
DimPlot(mb_mt, reduction="umap", group.by="seurat_clusters")
# umap colored by sample id
DimPlot(mb_mt, reduction="umap", group.by="SampleID")
DimPlot(mb_mt, reduction="umap", group.by="CellType")
dev.off()

pdf('figures/atac_cool_stuff_tsne.pdf', width=10, height=10)
# tsne colored by clusters
DimPlot(mb_mt, reduction="tsne", group.by="seurat_clusters")
# tsne colored by sample id
DimPlot(mb_mt, reduction="tsne", group.by="SampleID")
DimPlot(mb_mt, reduction="tsne", group.by="CellType")

dev.off()


save(mb_mt, file='combined_atac.rda')


# analyze activity
# need to set default assay to ATAC as opposed to activity
DefaultAssay(mb_mt) <- 'ACTIVITY'
mb_mt <- mb_mt %>%
	FindVariableFeatures() %>%
	NormalizeData() %>%
	ScaleData() %>%
	RunPCA(dims=1:100) %>%
	RunTSNE(reduction="pca", dims=1:50) %>%
	RunUMAP(reduction="pca", dims=1:50) %>% 
	FindNeighbors(reduction="pca", dims=1:50)

# find clusters
mb_mt <- FindClusters(mb_mt, resolution=1)

# known marker genes of interest
known_markers <- c("Myog", "Myod1", "Id3", "Ckm", "Mybph", "Myh3", "Ttn")

# umap colored by clusters
pdf("figures/activity_cool_stuff.pdf", width=10, height=10)
DimPlot(mb_mt, reduction="umap", group.by="seurat_clusters")
DimPlot(mb_mt, reduction="umap", group.by="SampleID")
DimPlot(mb_mt, reduction="umap", group.by="CellType")
FeaturePlot(mb_mt, features='Myog', reduction='umap')
FeaturePlot(mb_mt, features='Myod1', reduction='umap')
FeaturePlot(mb_mt, features='Id3', reduction='umap')
FeaturePlot(mb_mt, features='Ckm', reduction='umap')
dev.off()

save(mb_mt, file='combined_atac.rda')

# make violin plots of open chromatin regions by cluster
pdf("figures/violin_plot.pdf")
VlnPlot(
	mb_mt, 
	features=known_markers, 
	group.by="seurat_clusters", 
	split.by="CellType",
	pt.size=0,
	ncol=2
)
dev.off()

# make ridge plots for mt,mb genes by cluster
pdf('figures/ridge_plot.pdf')
RidgePlot(
	mb_mt,
	features=known_markers,
	group.by="seurat_clusters",
	split.by="CellType",
	# ncol=2
)
dev.off()

# find marker genes for clusters for activity 
my_markers <- FindAllMarkers(mb_mt, only.pos=T)

pdf("figures/marker_violin_plot.pdf", width=10, height=10)
clusters <- unique(mb_mt$seurat_clusters)
for(i in 1:length(clusters)){
	cluster_markers <- my_markers %>% 
		subset(cluster==clusters[i]) %>%
		top_n(n=10, wt=avg_logFC)
	p <- VlnPlot(mb_mt, features=cluster_markers$gene, group.by="seurat_clusters", split.by="CellType", pt.size=0, ncol=4)+
		ggtitle(paste0('cluster',i-1))
	print(p)
}

dev.off()

#NormalizeData
#ScaleData
#RunPCA
#RunTSNE
#RunUMAP
#FindNeighbors
#FindClusters


