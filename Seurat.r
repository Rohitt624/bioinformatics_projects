library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC")
# Load the HPC dataset
rnaseq.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC/Data")

# Pre-processing
seurat.object <- CreateSeuratObject(counts = rnaseq.data, project = "HPC", min.cells = 3, min.features = 200)
rm(rnaseq.data)
seurat.object <- NormalizeData(object = seurat.object)

#Subset
group1 <- subset(seurat.object, subset = orig.ident == "HPC.E16.5")
group2 <- subset(seurat.object, subset = orig.ident == "HPC.P7")
group3  <- subset(seurat.object, subset = orig.ident == "HPC.P14")
group4 <- subset(seurat.object, subset = orig.ident == "HPC.Adult")

#Splitting by age
split.object <- SplitObject(seurat.object, split.by = "orig.ident")
rm(seurat.object)
merged.object <- merge(x = split.object$HPC.E16.5, y = list(split.object$HPC.P7, split.object$HPC.P14, split.object$HPC.Adult))
rm(split.object)

#Differential expression by condition
gene_expression <- FetchData(group1, vars = "Cd27")
group1[["Cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
Idents(group1) <- "Cd27_highvslow"
de.markers <- FindMarkers(group1, ident.1 = "High", ident.2 = "Low")
write.csv(de.markers, "Cd27_highvslow_markers_E16.5.csv")

#Processing
group1 <- FindVariableFeatures(object = group1)
group1 <- ScaleData(object = group1)
group1 <- RunPCA(object =group1)

#Integrate Layers (if necessary)
merged.object <- IntegrateLayers(object =merged.object, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 verbose = FALSE)
joined.object <- JoinLayers(merged.object)
rm(merged.object)

#Cluster
group1 <- FindNeighbors(object =group1, dims = 1:30)
group1 <- FindClusters(object =group1)
group1 <- RunUMAP(object =group1, dims = 1:30)
DimPlot(object =group1, reduction = "umap")
DimPlot(group1, group.by = c("ctnnal1_highvslow"))
DimPlot(group1, split.by = c("ctnnal1_highvslow"))

#Define High vs Low for a certain gene
gene_expression <- FetchData(joined.object, vars = "Ctnnal1")
joined.object[["ctnnal1_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Ctnnal1, na.rm = TRUE), "High", "Low")
median(gene_expression$Ctnnal1)
DimPlot(joined.object, group.by = c("ctnnal1_highvslow"))

gene_expression <- FetchData(joined.object, vars = "Cd27")
joined.object[["cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
DimPlot(joined.object, group.by = c("cd27_highvslow"))

gene_expression <- FetchData(joined.object, vars = "Cd34")
joined.object[["cd34_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd34, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd34)
DimPlot(joined.object, group.by = c("cd34_highvslow"))

#identify markers in all clusters
clusters.markers <- FindAllMarkers(joined.object)
write.csv(clusters.markers, "allHPCmarkers.csv")
DoHeatmap(joined.object, features = c("Flt3", "Slamf1", "Cd48")) + NoLegend()
Adult.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(joined.object, features = top10$gene) + NoLegend()


#plots
VlnPlot(joined.object, features = c("Cd34", "Cd38", "Ptprc", "Thy1", "Itga6", "Kit", "Ly6a", "Flt3", "Slamf1", "Cd48"))
VlnPlot(joined.object, features = c("Sell"), split.by = c("cd27_highvslow"))
RidgePlot(joined.object, features = c("Cd34"))
FeaturePlot(HSC_joined, features = c("Sell"))
DimPlot(object =group1, reduction = "umap")
DimPlot(group1, group.by = c("ctnnal1_highvslow"))
DimPlot(group1, split.by = c("ctnnal1_highvslow"))