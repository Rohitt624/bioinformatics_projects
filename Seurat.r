library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
options(future.globals.maxSize = 2000 * 1024^2)

#Set working directory
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC")
# Load the HPC dataset
rnaseq.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC/Data")

# Pre-processing
seurat.object <- CreateSeuratObject(counts = rnaseq.data, project = "HPC", min.cells = 3, min.features = 200)
rm(rnaseq.data)
seurat.object <- NormalizeData(object = seurat.object)

#Subset (can be done before or after processing+clustering)
group1 <- subset(seurat.object, subset = orig.ident == "HPC.E16.5")
group2 <- subset(seurat.object, subset = orig.ident == "HPC.P7")
group3  <- subset(seurat.object, subset = orig.ident == "HPC.P14")
group4 <- subset(seurat.object, subset = orig.ident == "HPC.Adult")

#Splitting by age (can be done before or after processing+clustering)
split.object <- SplitObject(seurat.object, split.by = "orig.ident")
rm(seurat.object)
merged.object <- merge(x = split.object$HPC.E16.5, y = list(split.object$HPC.P7, split.object$HPC.P14, split.object$HPC.Adult))
rm(split.object)

#Differential expression by condition (can be done before or after processing+clustering)
gene_expression <- FetchData(group1, vars = "Ctnnal1")
write.csv(gene_expression, "ctnnal1.csv")
group1[["Cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
group1[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.67), "Cd34High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34, 0.33), "Cd34Low", "Cd34Mid"))

Idents(group1) <- "Cd27_highvslow" #set which metadata column it looks at for idents in the next line
de.markers <- FindMarkers(group1, ident.1 = "High", ident.2 = "Low")
write.csv(de.markers, "Cd27_highvslow_markers_E16.5.csv")

#Processing
seurat.object <- FindVariableFeatures(object = seurat.object)
seurat.object <- ScaleData(object = seurat.object)
seurat.object <- RunPCA(object =seurat.object)

#Integrate Layers (if you merged objects into multiple layers)
merged.object <- IntegrateLayers(object =merged.object, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 verbose = FALSE)
joined.object <- JoinLayers(merged.object)
rm(merged.object)

#Cluster
seurat.object <- FindNeighbors(object =seurat.object, dims = 1:30)
seurat.object <- FindClusters(object =seurat.object)
seurat.object <- RunUMAP(object =seurat.object, dims = 1:30)
#Generate UMAP to visualize clusters
DimPlot(object =seurat.object, reduction = "umap") #regular UMAP
DimPlot(seurat.object, group.by = c("ctnnal1_highvslow")) #UMAP by predefined group instead of cluster
DimPlot(seurat.object, split.by = c("ctnnal1_highvslow")) #Splits UMAP by these predefined groups

#Define High vs Low for a certain gene (can also be done before processing)
gene_expression <- FetchData(joined.object, vars = "Ctnnal1")
joined.object[["ctnnal1_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Ctnnal1, na.rm = TRUE), "High", "Low")
median(gene_expression$Ctnnal1)
group4[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.67), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.33), "Cd27Low", "Cd27Mid"))
DimPlot(joined.object, group.by = c("ctnnal1_highvslow"))

gene_expression <- FetchData(joined.object, vars = "Cd27")
joined.object[["cd27_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd27, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd27)
DimPlot(joined.object, group.by = c("cd27_highvslow"))

gene_expression <- FetchData(joined.object, vars = "Cd34")
joined.object[["cd34_highvslow"]] <- ifelse(gene_expression > median(gene_expression$Cd34, na.rm = TRUE), "High", "Low")
median(gene_expression$Cd34)
DimPlot(joined.object, group.by = c("cd34_highvslow"))

#Run clustify to annotate clusters if you have prelabeled scRNAseq data or bulk RNA-seq data
#Load reference data (the following is for bulk seq data)
ref <- read.csv("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv")
ref <- ref[, -1] #only necessary if bulk seq data has ensembl id before gene names
ref2 <- ref[, -1]
table(duplicated(ref2$Gene.name)) #check for duplicate gene names
rownames(ref2) <- make.names(ref$Gene.name, unique = TRUE) #forces gene names to become unique
#Calculate (it will automatically add a metadata column with the label)
res2 <- clustify(
  input = seurat.object,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)
#Plot results
DimPlot(res2, group.by = c("type"))
#Subset one population (you can recluster this subset for more analysis)
group2 <- subset(res2, subset = type %in% c("CMP_2", "CMP_3"))

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

#Applying annotation
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat.object)
#Merge clusters
idents[idents == "2" | idents == "3"] <- "2_3"
Idents(seurat_object) <- idents
#Rename clusters
new_idents <- seurat_object@meta.data[, "type"]
seurat_object <- RenameIdents(seurat_object, new_idents)

#Plots
VlnPlot(joined.object, features = c("Cd34", "Cd38", "Ptprc", "Thy1", "Itga6", "Kit", "Ly6a", "Flt3", "Slamf1", "Cd48"))
VlnPlot(joined.object, features = c("Sell"), split.by = c("cd27_highvslow"))
RidgePlot(joined.object, features = c("Cd34"))
FeaturePlot(HSC_joined, features = c("Sell"))
DimPlot(object =group1, reduction = "umap")
DimPlot(group1, group.by = c("ctnnal1_highvslow"))
DimPlot(group1, split.by = c("ctnnal1_highvslow"))