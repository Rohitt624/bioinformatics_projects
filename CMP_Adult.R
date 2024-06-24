library(clustifyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(AUCell)
library(GSEABase)
library(scuttle)
library(scran)
library(scater)

options(future.globals.maxSize = 2000 * 1024^2)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")

#Load and pre-process scRNAseq data --------------------
rnaseq.data <- Read10X(data.dir="C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")
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
umap <- DimPlot(object =seurat.object, reduction = "umap")
markers <- FindAllMarkers(seurat.object)



#Add some metadata column ----------------------------------
gene_expression <- FetchData(seurat.object, vars = "Cd27")
seurat.object[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.8), "Cd27High", 
                                     ifelse(gene_expression > 0, "Cd27Mid", "Cd27Low"))
DimPlot(object =seurat.object, reduction = "umap", group.by = "Cd27_highvslow")
gene_expression <- FetchData(seurat.object, vars = "Cd34")
seurat.object[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.8), "Cd34High", 
                                     ifelse(gene_expression >0, "Cd34Mid", "Cd34Low"))
DimPlot(object =seurat.object, reduction = "umap", group.by = "Cd34_highvslow")
gene_expression <- FetchData(seurat.object, vars = "Sell")
seurat.object[["Cd62L_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Sell, 0.8), "Cd62LHigh", 
                                      ifelse(gene_expression > 0, "Cd62LMid", "Cd62LMidLow"))
DimPlot(object =seurat.object, reduction = "umap", group.by = "Cd62L_highvslow")





#Visualization ----------------------------------------------
FeaturePlot(seurat.object, features = c("H2afy", "Procr", "mt-Nd1", "Pdzph1", "Rhd", "Csf1r", "Clec4a2", "Gp1bb", "Gzmb", "Gm15915", "Cd34", "Fcgr3"))
VlnPlot(seurat.object, features = c("Cd27", "Cd34", "Sell")) 
RidgePlot(seurat.object, features = c("Itga2b", "Eng")) 

#Further processing -------------------------------------------
#Subset one population
group1 <- subset(seurat.object, idents = "CMP")
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
FeaturePlot(group1, features = c("Cd34", "Sell"))
VlnPlot(group1, features = c("Cd27", "Cd34", "Sell")) 
RidgePlot(group1, features = c("Cd27", "Cd34", "Sell")) 

gene_expression <- FetchData(group1, vars = "Cd27")
group1[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.7), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.3), "Cd27Low", "Cd27Mid"))
DimPlot(object =group1, reduction = "umap", group.by = "Cd27_highvslow")
gene_expression <- FetchData(group1, vars = "Cd34")
group1[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.7), "Cd34High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34, 0.3), "Cd34Low", "Cd34Mid"))
DimPlot(object =group1, reduction = "umap", group.by = "Cd34_highvslow")
gene_expression <- FetchData(group1, vars = "Sell")
group1[["Cd62L_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Sell, 0.7), "Cd62LHigh", 
                                     ifelse(gene_expression < quantile(gene_expression$Sell, 0.3), "Cd62LLow", "Cd62LMid"))
DimPlot(object =group1, reduction = "umap", group.by = "Cd62L_highvslow")



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
ref <- read.csv("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Corrected.csv")
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
gCMP <- c("D930028M14Rik","Ermap","Atp1b2","Sox5","Kel","Zg16","Slc25a21","Ptpn14","Tmem246","Tspan8","Prss50",
          "Sphk1","Trib2","Ccl17","Tspo2","Mogat2","Mmp14","Aqp1","Car1","Gm15915","Ces2g","Mt2","Il1rl1","Col5a1",
          "Klf1")
gGMP <- c("Cybb","Epx","Dcbld1","Orm1","Prom1","Prg3","Glt1d1","Olr1","Nlrp12","B4galt6","Rag2","Hsd11b1",
          "Rhoj","Far2","Gda","Gm5294","Mtus1","Rnase10","Kcnk6","Alcam","Celsr3","Tmem184c","Cd5","Prg2","Vcam1")
gPG <- c("Ly6c1","Chd7","Fcnb","Ly6c2","Fam181b","Mss51","Wfdc21","Cracr2b","Tmem25","Gpx3","Cebpb","S100a8",
         "Oosp1","Gm16104","Gm7308","E030030I06Rik","Thy1","Gm4853","Gm13736","Rbm3os","Il11ra2","Gm11934",
         "Gm12286","Gm13186","Sptbn1")
gPM <- c("Gpr37l1","Gzmb","Pcsk1n","Nnat","Cacnb3","Mt3","S100b","Nags","Cpe","X6330403K07Rik","Cma1","Sparc",
          "Tcrg.C4","Chrna1os","Hes1","Plekhb1","Ccp110","Tagln3","Cd7","Myo1e","Tek","Gfap","Klf1","Ces2g","Pdcd1")
gMEP <- c("Slc17a4","Gm33104","Gsdma","Ntn4","Car1","Aldh1a7","Zpbp2","Pde10a","Pf4","F2rl2","Cplx2","Ces2g","Klf1",
          "Samd14","Aldh1a1","Epor","Ppbp","Slamf1","Col1a2","C5ar2","Trbv13.2","Dcaf12l1","Epdr1","F2r","Optn")
gMkP <- c("Cd41", "Vwf", "Gata2", "Gata1", "Gp1bb")
HSCSets <- GeneSet(gHSC, setName="HSC")
CMPSets <- GeneSet(gCMP, setName="CMP")
GMPSets <- GeneSet(gGMP, setName="GMP")
PGSets <- GeneSet(gPG, setName="PG")
PMSets <- GeneSet(gPM, setName="PM")
MEPSets <- GeneSet(gMEP, setName="MEP")
MkPSets <- GeneSet(gMkP, setName="MkP")

geneSets <- GeneSetCollection(c(HSCSets, CMPSets, GMPSets, PGSets, PMSets, MEPSets, MkPSets))

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
hvg_complete <- getTopHVGs(sce_complete, n = 2000)

# Run PCA on HVGs
sce_complete <- runPCA(sce_complete, subset_row = hvg_complete)

# Run UMAP on PCA
sce_complete <- runUMAP(sce_complete, dimred="PCA", verbose = T, min_dist = 0.2, n_neighbors = 30, name= "Uncorrected_UMAP")

# Plot UMAP
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "HSC") 
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "CMP")
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "GMP")
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "PG")
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "PM")
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "MEP")
plotReducedDim(sce_complete, dimred = "UMAP", colour_by = "MkP")

#Applying annotations
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
seurat.object[["idents"]] <- Idents(seurat.object)
#Rename clusters
cluster.ids <- c("HSC", "unlabeled" , "PM/PG", "MEP", "unlabeled", "CMP", "unlabeled", "MkP", "PG/GMP", "PM", "MEP", "unlabeled")
names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat.object, reduction = "umap")

group1 <- subset(seurat.object, idents = c("CMP", "PM", "MEP"))
DimPlot(object =group1, reduction = "umap")
FeaturePlot(group1, features = c("Cd27", "Cd34"))
VlnPlot(group1, features = c("Cd27", "Cd34", "Sell"))

# Get the top 10 genes from each cluster
top10_genes <- FindAllMarkers(group1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_genes <- top10_genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Create a vector of the gene names
genes_use <- top10_genes$gene
# Generate the heatmap
DoHeatmap(group1, features = genes_use)



#Applying manual annotation ----------------
Idents(seurat.object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat.object)
#Rename clusters
cluster.ids <- c("HSC/MPP/LMPP", "CMP", "GP", "CFU-E", "GMP", "CFU-E", "MEP", "MkP", "GMP", "Mast Cell Progenitor", "CFU-E", "MPP")
names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
DimPlot(object =seurat.object, reduction = "umap")




idents <- Idents(seurat.object)
seurat.object[["idents"]] <- idents
