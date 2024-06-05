library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")

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

FeaturePlot(seurat.object, features = c("Cd34", "Fcgr3"))
VlnPlot(seurat.object, features = c("Cd34", "Cd27", "Sell", "Flt3")) 

#Run clustify to annotate clusters

#Load reference data
ref <- read.csv("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv")
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

VlnPlot(res2, features = c("Ly6a", "Kit"))


#Subset one population
group2 <- subset(res2, subset = type %in% c("CMP_2", "CMP_3"))
group2 <- NormalizeData(object = group2)
#Processing
group2 <- FindVariableFeatures(object = group2)
group2 <- ScaleData(object = group2)
group2 <- RunPCA(object =group2)
#Cluster
group2 <- FindNeighbors(object =group2, dims = 1:30)
group2 <- FindClusters(object =group2)
group2 <- RunUMAP(object =group2, dims = 1:30)
DimPlot(object =group2, reduction = "umap")
#Calculate
res3 <- clustify(
  input = group2,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)
DimPlot(res3, group.by = c("type"))

table(res2$type)
table(res3$type)

clusters.markers <- FindAllMarkers(seurat.object)
write.csv(clusters.markers, "allmarkers.csv")

#Applying annotation
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat.object)
#Rename clusters
cluster.ids <- c("HSC/MPP", "CMP", "GMP", "E-MEP", "GMP", "E-MEP", "MK-MEP", "MK-MEP", "GMP", "GMP", "E-MEP", "MPP")
names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat.object, reduction = "umap")
