library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data")

#Load and pre-process scRNAseq data
rnaseq.data <- Read10X(data.dir="C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")
seurat.object <- CreateSeuratObject(counts = rnaseq.data, project = "MPC", min.cells = 3, min.features = 200)
rm(rnaseq.data)
seurat.object <- NormalizeData(object = seurat.object)
#Processing
seurat.object <- FindVariableFeatures(object = seurat.object)
seurat.object <- ScaleData(object = seurat.object)
seurat.object <- RunPCA(object =seurat.object)
#Cluster
seurat.object <- FindNeighbors(object =seurat.object, dims = 1:30)
seurat.object <- FindClusters(object =seurat.object)
seurat.object <- RunUMAP(object =seurat.object, dims = 1:30)
DimPlot(object =seurat.object, reduction = "umap")

#Run clustify to annotate clusters

#Load reference data
ref <- read.csv("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Corrected.csv")
ref <- ref[, -1]
ref2 <- ref[, -1]
table(duplicated(ref2$Gene.name))
rownames(ref2) <- make.names(ref$Gene.name, unique = TRUE)
#Calculate
res2 <- clustify(
  input = seurat.object,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)

#plot results
DimPlot(res2, group.by = c("type"))

#Subset one population
group1 <- subset(res2, subset = type == "PG_3")
group1 <- NormalizeData(object = group1)
#Processing
group1 <- FindVariableFeatures(object = group1)
group1 <- ScaleData(object = group1)
group1 <- RunPCA(object =group1)
#Cluster
group1 <- FindNeighbors(object =group1, dims = 1:30)
group1 <- FindClusters(object =group1)
group1 <- RunUMAP(object =group1, dims = 1:30)
DimPlot(object =group1, reduction = "umap")
#Calculate
res2 <- clustify(
  input = group1,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)
DimPlot(group1, group.by = c("type"))


