library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(biomaRt)
options(future.globals.maxSize = 4000 * 1024^2)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/Cite-seq/GSE145491/")

# Load the dataset and create Seurat Object with two assays ---------------------------------
data <- as.sparse(read.table(file = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/Cite-seq/GSE145491/GSE145491_UMIsMatrix.txt.gz",
                               sep = "", header = TRUE, row.names = 1))
adtfeat <- c("CD9", "CD41", "CD48", "CD55", "CD105", "CD115", "CD135", "CD150", "CXCR4", "ESAM", "CD4")
adtmat <- data[adtfeat, ]
allfeat <- rownames(data)
rnafeat <- setdiff(allfeat, adtfeat)
rnamat <- data[rnafeat, ]
rm(data, adtfeat, allfeat)

# Replace RNA Ensembl IDs with gene IDs (This section maxes out my RAM so I can't test it. Just skip)----------------
ensembl <- useEnsembl(biomart = "genes")
listDatasets(ensembl)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_ids <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  values = rnafeat,
  mart = mart
)

id_list <- split(gene_ids[["mgi_symbol"]], gene_ids[["ensembl_gene_id"]])
rnaname <- character(nrow(rnamat))
for (i in seq_len(nrow(rnamat))) {
  ensembl_id <- rownames(rnamat)[i]
  if (ensembl_id %in% names(id_list)) {
    rnaname[i] <- id_list[[ensembl_id]]
  } else{
    rnaname[i] <- ensembl_id
  }
}
rnaname[rnaname == ""] <- NA

rownames(rnamat) <- rnaname
rm(ensembl, gene_ids, id_list, mart, ensembl_id, i, rnaname)
rownames(rnamat) <- make.names(rnamat, unique = TRUE)


#Create Seurat Object -------------------------------------------------------------------
seurat_object <- CreateSeuratObject(counts = rnamat)
seurat_object[["ADT"]] <- CreateAssay5Object(counts = adtmat)
rm(adtmat, rnamat, rnafeat)
Assays(seurat_object)
seurat_object <- NormalizeData(object = seurat_object, assay = "RNA")
seurat_object <- NormalizeData(object = seurat_object, assay = "ADT", method = "CLR")



#Check features in each assay
rownames(seurat_object[["ADT"]])
rownames(seurat_object[["RNA"]])

#test ADT data by choosing any two ADT features
FeatureScatter(seurat_object, feature1 = "adt_CD41", feature2 = "adt_CD105")

# List the current default assay
DefaultAssay(seurat_object)
#Change default assay
DefaultAssay(seurat_object) <- "RNA" #change it to either "RNA" or "ADT"
DefaultAssay(seurat_object)

#Process RNA data
seurat_object <- FindVariableFeatures(object = seurat_object)
seurat_object <- ScaleData(object = seurat_object)
seurat_object <- RunPCA(object =seurat_object)
seurat_object <- FindNeighbors(object =seurat_object, dims = 1:30)
seurat_object <- FindClusters(object =seurat_object)
seurat_object <- RunUMAP(object =seurat_object, dims = 1:30)
DimPlot(seurat_object)


# Visualize levels for RNA and protein 
DefaultAssay(seurat_object) <- "ADT"
p1 <- FeaturePlot(seurat_object, "CD41", cols = c("lightgrey", "darkgreen")) + ggtitle("CD41 protein")
DefaultAssay(seurat_object) <- "RNA"
p2 <- FeaturePlot(seurat_object, "Cd41") + ggtitle("Cd41 RNA")
# place plots side-by-side
p1 | p2

#alternative to above
Key(seurat_object[["RNA"]]) #use this key so that you don't need to switch default assay every time
Key(seurat_object[["ADT"]])
p1 <- FeaturePlot(seurat_object, "adt_CD34", cols = c("lightgrey", "darkgreen")) + ggtitle("CD34 protein")
p2 <- FeaturePlot(seurat_object, "rna_Cd34") + ggtitle("CD34 RNA")

#Similarly, to identify surface markers
VlnPlot(seurat_object, "adt_CD11b")
#You can also FindMarkers for each assay
adt_markers <- FindMarkers(seurat_object, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(seurat_object, ident.1 = 6, assay = "RNA")
#And create ADT scatter plots while showing cluster
FeatureScatter(seurat_object, feature1 = "adt_CD150", feature2 = "adt_CD105")
#Can "gate" cells using HoverLocator and FeatureLocator
# view relationship between protein and RNA
FeatureScatter(seurat_object, feature1 = "adt_CD34", feature2 = "rna_CD34")

