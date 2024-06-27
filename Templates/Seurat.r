library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2) #increases max amount of RAM the program can use (this is 2GB)

#Set working directory ---------------------------------------
setwd(".../output/")
# Load the HPC dataset
data <- Read10X(data.dir = ".../Data/") #folder with barcodes, matrix, and features files

# Pre-processing-------------------------------------
#min.cells filters for features that are present in at least that many cells
#min.features filters for cells that have at least that many features
seurat_object <- CreateSeuratObject(counts = data, project = "HPC", min.cells = 3, min.features = 200)
rm(data)
seurat_object <- NormalizeData(object = seurat_object)

#Subset (can be done before or after processing+clustering)----------------------------------
g1 <- subset(seurat_object, subset = orig.ident == "HPC.E16.5")
g2 <- subset(seurat_object, subset = orig.ident == "HPC.P7")
g3  <- subset(seurat_object, subset = orig.ident == "HPC.P14")
g4 <- subset(seurat_object, subset = orig.ident == "HPC.Adult")

#Splitting by age (can be done before or after processing+clustering)
#creates layers for the different original identities
split.object <- SplitObject(seurat_object, split.by = "orig.ident")
rm(seurat_object)
merged.object <- merge(x = split.object$HPC.E16.5, y = list(split.object$HPC.P7, split.object$HPC.P14, split.object$HPC.Adult))
rm(split.object)

#Differential expression by condition (can be done before or after processing+clustering)
#split by median
gene_expression <- FetchData(g1, vars = "Cd27")
g1[["Cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
DimPlot(object =seurat.object, reduction = "umap", group.by = "Cd27_highvslow")

#split at defined percentile
gene_expression <- FetchData(group1, vars = "Cd27")
g1[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.7), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.3), "Cd27Low", "Cd27Mid"))
DimPlot(object =group1, reduction = "umap", group.by = "Cd27_highvslow")

Idents(g1) <- "Cd27_highvslow" #set which metadata column it looks at for idents in the next line
de.markers <- FindMarkers(g1, ident.1 = "High", ident.2 = "Low")
write.csv(de.markers, "Cd27_highvslow_markers_E16.5.csv")

#Processing------------------------------------------
seurat_object <- FindVariableFeatures(object = seurat_object)
seurat_object <- ScaleData(object = seurat_object)
seurat_object <- RunPCA(object =seurat_object)

#Integrate Layers (if you merged objects into multiple layers)
merged.object <- IntegrateLayers(object =merged.object, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 verbose = FALSE)
joined.object <- JoinLayers(merged.object)
rm(merged.object)

#Cluster
seurat_object <- FindNeighbors(object =seurat_object, dims = 1:30)
seurat_object <- FindClusters(object =seurat_object)
seurat_object <- RunUMAP(object =seurat_object, dims = 1:30)
#Generate UMAP to visualize clusters
DimPlot(object =seurat_object, reduction = "umap") #regular UMAP
DimPlot(seurat_object, group.by = c("ctnnal1_highvslow")) #UMAP by predefined group instead of cluster
DimPlot(seurat_object, split.by = c("ctnnal1_highvslow")) #Splits UMAP by these predefined groups

#Identify markers in all clusters
clusters.markers <- FindAllMarkers(joined.object)
write.csv(clusters.markers, "allHPCmarkers.csv")
DoHeatmap(joined.object, features = c("Flt3", "Slamf1", "Cd48")) + NoLegend()
Adult.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(joined.object, features = top10$gene) + NoLegend()


#Applying manual annotation ---------------------------------------------
Idents(seurat_object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat_object)
#Rename clusters
cluster.ids <- c("HSC/MPP/LMPP", "CMP", "GP", "CFU-E", "GMP", "CFU-E", "MEP", "MkP", "GMP", "Mast Cell Progenitor", "CFU-E", "MPP")
names(cluster.ids) <- levels(seurat.object)
seurat_object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat_object, reduction = "umap")

#Applying annotation
Idents(seurat_object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat_object)
#Merge clusters
idents[idents == "2" | idents == "3"] <- "2_3"
Idents(seurat_object) <- idents
#Rename clusters
new_idents <- seurat_object@meta.data[, "type"]
seurat_object <- RenameIdents(seurat_object, new_idents)

#Different types of plots --------------------------------------------------
VlnPlot(joined.object, features = c("Cd34", "Cd38", "Ptprc", "Thy1", "Itga6", "Kit", "Ly6a", "Flt3", "Slamf1", "Cd48"))
VlnPlot(joined.object, features = c("Sell"), split.by = c("cd27_highvslow"))
RidgePlot(joined.object, features = c("Cd34"))
FeaturePlot(HSC_joined, features = c("Sell"))
DimPlot(object =g1, reduction = "umap")
DimPlot(g1, group.by = c("ctnnal1_highvslow"))
DimPlot(g1, split.by = c("ctnnal1_highvslow"))