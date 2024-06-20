library(Seurat)
library(clustifyr)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC")

#Run clustifyr to annotate clusters if you have prelabeled scRNAseq data or bulk RNA-seq data
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