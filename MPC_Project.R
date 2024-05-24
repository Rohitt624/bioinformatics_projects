library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq")
# Load the HSC dataset
Adult.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/GSM4645164_Adult_WT")

# Pre-processing
Adult <- CreateSeuratObject(counts = Adult.data, project = "MPC", min.cells = 3, min.features = 200)
rm(Adult.data)

Adult <- NormalizeData(object = Adult)
Adult <- FindVariableFeatures(object = Adult)
Adult <- ScaleData(object = Adult)
Adult <- RunPCA(object =Adult)

#Cluster
Adult <- FindNeighbors(object =Adult, dims = 1:30)
Adult <- FindClusters(object =Adult)
Adult <- RunUMAP(object =Adult, dims = 1:30)
DimPlot(object =Adult, reduction = "umap")


#identify markers in all clusters
Adult.markers <- FindAllMarkers(Adult)
write.csv(Adult.markers, "allAdultmarkers.csv")
DoHeatmap(Adult, features = c("Ctla2a", "F2r", "Mpl", "Dntt", "Zfpm1", "Shank3", "Samd14", "Mfsd2b", "Ache", "Pf4")) + NoLegend()
Adult.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Adult, features = top10$gene) + NoLegend()

VlnPlot(Adult, features = c("Cd34", "Cd38", "Ptprc", "Thy1", "Itga6"))
VlnPlot(Adult, features = c("Kit", "Ly6a", "Flt3", "Slamf1", "Cd48"))
VlnPlot(Adult, features = c("Ctla2a", "F2r", "Mpl", "Dntt", "Zfpm1", "Shank3"))
RidgePlot(Adult, features = c("Ctla2a", "F2r", "Mpl", "Dntt", "Zfpm1", "Shank3"))
FeaturePlot(Adult, features = c("Flt3"))


#plots
VlnPlot(HSC_joined, features = c("Cd27", "Cd34", "Esam", "Sell", "Ctnnal1"))
VlnPlot(HSC_joined, features = c("Cd34"))
FeaturePlot(HSC_joined, features = c("Sell"))