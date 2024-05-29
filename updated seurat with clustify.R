library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC")

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
ref <- read.csv("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv", row.names = 1)
#Calculate
res <- clustify(
  input = seurat.object,
  ref_mat = ref,
  cluster_col = "seurat_clusters",
  query_genes = VariableFeatures(seurat.object),
  obj_out = TRUE
)
#view results
cor_to_call(res)
#plot results
plot_best_call(
  cor_mat = res,
  cluster_col = "seurat_clusters"
)

gene_names <- rownames(seurat.object[["RNA"]]@counts)
VlnPlot(seurat.object, features = "Cd27")
VariableFeatures(seurat.object)
