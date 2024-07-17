library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)

#Prepare the seurat object of your new dataset
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")

#Load and pre-process scRNAseq data --------------------
rnaseq.data <- Read10X(data.dir="C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT/Data Files")
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

saveRDS(seurat.object, "adultwt_seurat.rds")

#Load in pre-annotated dataset
annotation <- readRDS("locationoffile.rds")


py_runfile("cellproject.py")

  

