library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(AUCell)
library(GSEABase)
options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")

#Load and pre-process scRNAseq data --------------------
rnaseq.data <- Read10X(data.dir="C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")
seurat.object <- CreateSeuratObject(counts = rnaseq.data, project = "MPC")
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
umap <- DimPlot(object =seurat.object, reduction = "umap")

#Visualization ----------------------------------------------
FeaturePlot(seurat.object, features = c("Cd34", "Fcgr3"))
VlnPlot(seurat.object, features = c("Cd27", "Cd34", "Sell")) 
RidgePlot(seurat.object, features = c("Itga2b", "Eng")) 

#Further processing -------------------------------------------
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

#clustifyr -------------------------
ref <- read.csv("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv")
ref <- ref[, -1] #only necessary if bulk seq data has ensembl id before gene names
ref2 <- ref[, -1]
table(duplicated(ref2$Gene.name)) #check for duplicate gene names
rownames(ref2) <- make.names(ref$Gene.name, unique = TRUE) #forces gene names to become unique
#Calculate (it will automatically add a metadata column with the label)
res <- clustify(
  input = seurat.object,
  ref_mat = ref2,
  cluster_col = "seurat_clusters",
  obj_out = TRUE
)
DimPlot(res, group.by = c("type"))


# Generate the heatmap
DoHeatmap(s2, features = top_genes$gene, slot = )
DimPlot(group1, group.by = c("Cd34_highvslow"))

#AUCell -----------------------------------
gHSC <- c("Procr", "Pdzk1ip1", "Ltb", "Mllt3", "Ifitm1", "Gimap1", "Gimap6", "Limd2", "Trim47", "Neil2", "Vwf",
          "Pde1b", "Neo1", "Sqrdl", "Sult1a1", "Cd82", "Ramp2", "Ubl3", "Ly6a", "Cdkn1c", "Fgfr3", "Cldn10", "Ptpn14", 
          "Mettl7a1", "Smtnl1", "Ctsf", "Gstm1", "Sox18", "Fads3")
gCMP <- c("Mpo", "Lars2", "Gm23935", "Rpl18a", "Gpx1", "Gm42418", "Gm26917", "Eef1a1", "Rps2", "Rpsa", "Rplp0", "Ppia", 
          "Ftl1", "Gapdh", "Eif5a", "Prtn3", "Rpl13", "Rps18", "Rps3", "Elane", "Ly6e", "Rps14", "Rps19", "H2-D1", 
          "Rack1")
HSCSets <- GeneSet(gHSC, setName="HSC")
CMPSets <- GeneSet(gCMP, setName="CMP")
geneSets <- GeneSetCollection(c(HSCSets, CMPSets))

# Convert Seurat object to SingleCellExperiment
sce_complete <- as.SingleCellExperiment(seurat.object)

# Calculate cell rankings
cell_rankings <- assay(sce_complete, "counts") %>% AUCell_buildRankings()

# Calculate signature scores
signature_scores_sce_complete <- AUCell_calcAUC(geneSets, cell_rankings)

# Add signature scores to colData
colData(sce_complete) <- cbind(colData(sce_complete), t(signature_scores_sce_complete@assays@data$AUC)[colnames(sce_complete),])

# Normalize counts
sce_complete <- logNormCounts(sce_complete)

# Get top HVGs
hvg_complete <- getTopHVGs(dec_complete, n = 2000)

# Run PCA on HVGs
sce_complete <- runPCA(sce_complete, subset_row = hvg_complete)

# Run UMAP on PCA
sce_complete <- runUMAP(sce_complete, dimred="PCA", verbose = T, min_dist = 0.2, n_neighbors = 30, name= "Uncorrected_UMAP")

# Plot UMAP
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "HSC") 
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "CMP")



#Applying manual annotation ----------------
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat.object)
#Rename clusters
cluster.ids <- c("HSC/MPP/LMPP", "CMP", "GP", "E-MEP 1", "GMP 1", "E-MEP 2", "MEP", "MK-MEP", "GMP 2", "MOP", "E-MEP 3", "MPP")
names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat.object, reduction = "umap")


idents <- Idents(seurat.object)
seurat.object[["idents"]] <- idents
