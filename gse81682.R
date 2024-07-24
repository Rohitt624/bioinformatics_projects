library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSE81682_HTSeq_counts.txt")


counts <- read.table("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSE81682_HTSeq_counts.txt/HTSeq_counts.txt", header=TRUE, row.names = 1)

seu <- CreateSeuratObject(counts, project = "test", assay = "RNA", min.cells = 0, min.features = 0)
seu <- NormalizeData(object = seu)
#Processing
seu <- FindVariableFeatures(object = seu)
seu <- ScaleData(object = seu)
seu <- RunPCA(object =seu)
#Cluster
seu <- FindNeighbors(object =seu, dims = 1:30)
seu <- FindClusters(object =seu)
seu <- RunUMAP(object =seu, dims = 1:30)
DimPlot(object =seu, reduction = "umap")
rm(counts)

VlnPlot(seu, features = c("ENSMUSG00000005672", "ENSMUSG00000075602"))
FeaturePlot(seu, features = c("ENSMUSG00000027611"))
mark <- FindAllMarkers(seu)
write.csv(mark, "markers.csv")

#Applying manual annotation ----------------
Idents(seu) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seu)
#Rename clusters
cluster.ids <- c("LMPP", "MPP3", "2", "3", "4", "5", "6", "CMP", "GMP", "9")
names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)
DimPlot(object =seu, reduction = "umap")
