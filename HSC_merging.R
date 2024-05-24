library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HSC")
# Load the HSC dataset
HSC.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HSC/HSC_P7")

# Pre-processing
HSC <- CreateSeuratObject(counts = HSC.data, project = "HSC", min.cells = 3, min.features = 200)
rm(HSC.data)
HSC_split <- SplitObject(HSC, split.by = "orig.ident")
rm(HSC)
HSC_merged <- merge(x = HSC_split$HSC.E16.5, y = list(HSC_split$HSC.P7, HSC_split$HSC.P14, HSC_split$HSC.Adult))
rm(HSC_split)
HSC_merged <- NormalizeData(object = HSC_merged)
HSC_merged <- FindVariableFeatures(object = HSC_merged)
HSC_merged <- ScaleData(object = HSC_merged)
HSC_merged <- RunPCA(object =HSC_merged)
HSC_merged <- IntegrateLayers(object =HSC_merged, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                              verbose = FALSE)
HSC_joined <- JoinLayers(HSC_merged)
rm(HSC_merged)

#Cluster
HSC_joined <- FindNeighbors(object =HSC_joined, dims = 1:30)
HSC_joined <- FindClusters(object =HSC_joined)
HSC_joined <- RunUMAP(object =HSC_joined, dims = 1:30)
DimPlot(object =HSC_joined, reduction = "umap")
DimPlot(HSC_joined, group.by = c("orig.ident"))
DimPlot(HSC_joined, split.by = c("orig.ident"))

#identify markers in all clusters
HSC.markers <- FindAllMarkers(HSC_joined)
write.csv(HSC.markers, "allHSCmarkers.csv")
DoHeatmap(HSC, features = top10$gene) + NoLegend()

#plots
VlnPlot(HSC_joined, features = c("Cd27", "Cd34", "Esam", "Sell", "Ctnnal1"))
VlnPlot(HSC_joined, features = c("Cd34", "Cd27", "Sell"))
RidgePlot(HSC_joined, features = c("Cd34", "Cd27", "Sell"))
FeaturePlot(HSC_joined, features = c("Sell"))

#Define High vs Low for a certain gene
gene_expression <- FetchData(HSC_joined, vars = "Ctnnal1")
HSC_joined[["ctnnal1_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Ctnnal1), "High", "Low")
median(gene_expression$Ctnnal1)
DimPlot(HSC_joined, group.by = c("ctnnal1_highvslow"))

gene_expression <- FetchData(HSC_joined, vars = "Cd27")
HSC_joined[["cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27), "High", "Low")
median(gene_expression$Cd27)
DimPlot(HSC_joined, group.by = c("cd27_highvslow"))

gene_expression <- FetchData(HSC_joined, vars = "Cd34")
HSC_joined[["cd34_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd34), "High", "Low")
median(gene_expression$Cd34)
DimPlot(HSC_joined, group.by = c("cd34_highvslow"))

#Subset
CD34high <- subset(HSC_joined, subset = cd34_highvslow == "high")
CD34low <- subset(HSC_joined, subset = cd34_highvslow == "Low")

cd34.markers <- FindMarkers(HSC_joined, cd34_highvslow == "High", cd34_highvslow == "Low")

