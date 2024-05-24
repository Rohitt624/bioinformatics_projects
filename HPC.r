library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC")
# Load the HPC dataset
HPC.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC/Data")

# Pre-processing
HPC <- CreateSeuratObject(counts = HPC.data, project = "HPC", min.cells = 3, min.features = 200)
rm(HPC.data)
HPC_split <- SplitObject(HPC, split.by = "orig.ident")
rm(HPC)
HPC_merged <- merge(x = HPC_split$HPC.E16.5, y = list(HPC_split$HPC.P7, HPC_split$HPC.P14, HPC_split$HPC.Adult))
rm(HPC_split)
HPC_merged <- NormalizeData(object = HPC_merged)
HPC_merged <- FindVariableFeatures(object = HPC_merged)
HPC_merged <- ScaleData(object = HPC_merged)
HPC_merged <- RunPCA(object =HPC_merged)

#Integrate Layers
HPC_merged <- IntegrateLayers(object =HPC_merged, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                              verbose = FALSE)
HPC_joined <- JoinLayers(HPC_merged)
rm(HPC_merged)

#Cluster
HPC_joined <- FindNeighbors(object =HPC_joined, dims = 1:30)
HPC_joined <- FindClusters(object =HPC_joined)
HPC_joined <- RunUMAP(object =HPC_joined, dims = 1:30)
DimPlot(object =HPC_joined, reduction = "umap")
DimPlot(HPC_joined, group.by = c("orig.ident"))
DimPlot(HPC_joined, split.by = c("orig.ident"))

#Define High vs Low for a certain gene
gene_expression <- FetchData(HPC_joined, vars = "Ctnnal1")
HPC_joined[["ctnnal1_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Ctnnal1, na.rm = TRUE), "High", "Low")
median(gene_expression$Ctnnal1)
DimPlot(HPC_joined, group.by = c("ctnnal1_highvslow"))

gene_expression <- FetchData(HPC_joined, vars = "Cd27")
HPC_joined[["cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
DimPlot(HPC_joined, group.by = c("cd27_highvslow"))

gene_expression <- FetchData(HPC_joined, vars = "Cd34")
HPC_joined[["cd34_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd34, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd34)
DimPlot(HPC_joined, group.by = c("cd34_highvslow"))


#use the FindMarkers function to calculate differential gene expression between the two layers:
HPC_joined[["identity"]] <- HPC_joined$orig.ident
de_results <- FindMarkers(HPC_joined$orig.ident, cells.1 = "HPC.P7", cells.2 = "HPC.Adult")


#identify markers in all clusters
HPC.markers <- FindAllMarkers(HPC_joined)
write.csv(HPC.markers, "allHPCmarkers.csv")
DoHeatmap(HPC_joined, features = c("Flt3", "Slamf1", "Cd48")) + NoLegend()


#plots
VlnPlot(HPC_joined, features = c("Cd34", "Cd38", "Ptprc", "Thy1", "Itga6", "Kit", "Ly6a", "Flt3", "Slamf1", "Cd48"))
VlnPlot(HPC_joined, features = c("Sell"), split.by = c("cd27_highvslow"))
RidgePlot(HPC_joined, features = c("Cd34"))
FeaturePlot(HSC_joined, features = c("Sell"))

#Subset
CD34high <- subset(HPC_joined, subset = cd34_highvslow == "high")
CD34low <- subset(HPC_joined, subset = cd34_highvslow == "Low")

cd34.markers <- FindMarkers(HPC_joined, cd34_highvslow == "High", cd34_highvslow == "Low")