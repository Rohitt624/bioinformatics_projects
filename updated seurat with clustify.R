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
VlnPlot(seurat.object, features = c("Cd27", "Cd34", "Sell")) 
RidgePlot(seurat.object, features = c("Itga2b", "Eng")) 

#Subset one population
group1 <- subset(seurat.object, subset = idents %in% c("CMP"))
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
DimPlot(object =group1, reduction = "umap", group.by = "idents")
FeaturePlot(group1, features = c("Cd34", "Fcgr3"))
VlnPlot(group1, features = c("Cd27", "Cd34", "Sell")) 
RidgePlot(group1, features = c("Cd27", "Cd34", "Sell")) 

gene_expression <- FetchData(group1, vars = "Cd27")
group1[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.7), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.3), "Cd27Low", "Cd27Mid"))
DimPlot(object =group1, reduction = "umap", split.by = "Cd27_highvslow")
gene_expression <- FetchData(group1, vars = "Cd34")
group1[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.7), "Cd34High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34, 0.3), "Cd34Low", "Cd34Mid"))
DimPlot(object =group1, reduction = "umap", split.by = "Cd34_highvslow")
gene_expression <- FetchData(group1, vars = "Sell")
group1[["Cd62L_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Sell, 0.7), "Cd62LHigh", 
                                     ifelse(gene_expression < quantile(gene_expression$Sell, 0.3), "Cd62LLow", "Cd62LMid"))
DimPlot(object =group1, reduction = "umap", split.by = "Cd62L_highvslow")

Idents(group1) <- "Cd34_highvslow"
s2.1 <- subset(group1, idents = "Cd34Low")
Idents(s2.1) <- "Cd62L_highvslow"
s2 <- subset(s2.1, idents = "Cd62LLow")
rm(s2)
Idents(group1) <- "Cd27_highvslow"
s3.1 <- subset(group1, idents = "Cd27Low")
Idents(s3.1) <- "Cd62L_highvslow"
s3 <- subset(s3.1, idents = "Cd62LLow")
rm(s3)

Idents(group1) <- "Cd62L_highvslow"
s2 <- subset(group1, idents = "Cd62LLow")
s2$PM <- ifelse(s2$Cd27_highvslow == 'Cd27Low' & s2$Cd34_highvslow == 'Cd34Low', 'BothLow', 
                ifelse(s2$Cd27_highvslow == "Cd27Low" & s2$Cd34_highvslow != "Cd34Low", "Cd27only",
                ifelse(s2$Cd27_highvslow != "Cd27Low" & s2$Cd34_highvslow == "Cd34Low", "Cd34only", 
                "Misc")))
DimPlot(s2, group.by = 'PM')
table(s2$PM)

ref <- read.csv("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv")
ref <- ref[, -1] #only necessary if bulk seq data has ensembl id before gene names
ref2 <- ref[, -1]
table(duplicated(ref2$Gene.name)) #check for duplicate gene names
rownames(ref2) <- make.names(ref$Gene.name, unique = TRUE) #forces gene names to become unique
#Calculate (it will automatically add a metadata column with the label)
res <- clustify(
  input = s2,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)

clusters.markers <- FindAllMarkers(s2)
write.csv(clusters.markers, "allmarkers.csv")

DimPlot(res, group.by = c("type"))

Idents(s2) <- "PM"
s2.markers <- FindAllMarkers(s2)
write.csv(s2.markers, "alls2markers.csv")
top_genes <- FindAllMarkers(s2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by("PM") %>%
  top_n(n = 10, wt = avg_logFC)

# Generate the heatmap
DoHeatmap(s2, features = top_genes$gene, slot = )
DimPlot(group1, group.by = c("Cd34_highvslow"))

#Applying annotation
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat.object)
#Rename clusters
cluster.ids <- c("HSC/MPP/LMPP", "CMP", "GP", "E-MEP", "GMP", "E-MEP", "MEP", "MK-MEP", "GMP", "MOP", "E-MEP", "MPP")
names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat.object, reduction = "umap")


idents <- Idents(seurat.object)
seurat.object[["idents"]] <- idents

