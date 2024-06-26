library(AUCell)
library(GSEABase)
library(scuttle)
library(scran)
library(scater)
library(Seurat)
setwd(".../output/") #replace this with the file path to whatever folder you're using

#AUCell -----------------------------------
#the HSC and MkP sets are from papers that I found, but the others are the most differentially expressed genes from the MPC bulk seq data
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
sce_complete <- as.SingleCellExperiment(seurat_object) #preferrably use a seurat object that you have already clustered if you want to use this for annotation

#from here to "Run UMAP on PCA" is copy pasted from the code used in this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11147777/
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
Idents(seurat_object) <- "seurat_clusters" #make sure the clusters are the identity metadata column
idents <- Idents(seurat_object) #moves these cluster identities into their own variable
#Rename clusters (decided based on which parts were most enriched in the DimPlots from above)
cluster.ids <- c("HSC", "unlabeled" , "PM/PG", "CMP", "unlabeled", "MEP", "unlabeled", "MkP", "PG/GMP", "PM", "MEP", "unlabeled")
names(cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, cluster.ids)
DimPlot(object =seurat_object, reduction = "umap")
