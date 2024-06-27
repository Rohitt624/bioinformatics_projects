library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
options(future.globals.maxSize = 4000 * 1024^2)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/Cite-seq/GSE268141_RAW/")

# Load the dataset
data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/Cite-seq/GSE268141_RAW/")
rownames(x = data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqA", replacement = "",
                                                         x = rownames(x = data[["Antibody Capture"]]))
# Initialize the Seurat object with the raw (non-normalized data).
seurat_object <- CreateSeuratObject(counts = data[["Gene Expression"]], min.cells = 3, min.features = 200)
seurat_object[["ADT"]] <- CreateAssay5Object(counts = data[["Antibody Capture"]][, colnames(x = seurat_object)])
Assays(seurat_object) #check assays in object
seurat_object <- NormalizeData(object = seurat_object, assay = "RNA")
seurat_object <- NormalizeData(object = seurat_object, assay = "ADT", method = "CLR")
rm(data)

FeatureScatter(seurat_object, feature1 = "adt_CD11b", feature2 = "adt_CD127")

# Extract a list of features measured in the ADT assay
rownames(seurat_object[["ADT"]])
rownames(seurat_object[["RNA"]])

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


# Now, we will visualize CD34 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(seurat_object) <- "ADT"
p1 <- FeaturePlot(seurat_object, "CD27", cols = c("lightgrey", "darkgreen")) + ggtitle("CD27 protein")
DefaultAssay(seurat_object) <- "RNA"
p2 <- FeaturePlot(seurat_object, "Cd27") + ggtitle("Cd27 RNA")
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

