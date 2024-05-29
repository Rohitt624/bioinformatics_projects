library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 2000 * 1024^2)

setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC/Output")
# Load the HPC dataset
rnaseq.data <- Read10X(data.dir = "C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/HPC/Data")

# Pre-processing
seurat.object <- CreateSeuratObject(counts = rnaseq.data, project = "HPC", min.cells = 3, min.features = 200)
rm(rnaseq.data)
seurat.object <- NormalizeData(object = seurat.object)

#Subset
group1 <- subset(seurat.object, subset = orig.ident == "HPC.E16.5")
group2 <- subset(seurat.object, subset = orig.ident == "HPC.P7")
group3  <- subset(seurat.object, subset = orig.ident == "HPC.P14")
group4 <- subset(seurat.object, subset = orig.ident == "HPC.Adult")

#Differential expression by CD27

gene_expression <- FetchData(group4, vars = "Cd27")
group4[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.8), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.2), "Cd27Low", "Cd27Mid"))
Idents(group4) <- "Cd27_highvslow"
de.markers <- FindMarkers(group4, ident.1 = "High", ident.2 = "Low")
write.csv(de.markers, "Cd27_highvslow_markers_Adult.csv")

#Differential expression by CD34

gene_expression <- FetchData(group4, vars = "Cd34")
group4[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.8), "Cd34High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34, 0.2), "Cd34Low", "Cd34Mid"))
Idents(group4) <- "Cd34_highvslow"
de.markers <- FindMarkers(group4, ident.1 = "High", ident.2 = "Low")
write.csv(de.markers, "Cd34_highvslow_markers_Adult.csv")

group4[["Cd27_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd27, 0.67), "Cd27High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd27, 0.33), "Cd27Low", "Cd27Mid"))

group4[["Cd34_highvslow"]] <- ifelse(gene_expression > quantile(gene_expression$Cd34, 0.67), "Cd34High", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34, 0.33), "Cd34Low", "Cd34Mid"))


group4[["Cd27vs34"]] <- ifelse(gene_expression < quantile(gene_expression$Cd27$Cd27, 0.2), "Cd27Low", 
                                     ifelse(gene_expression < quantile(gene_expression$Cd34$Cd34, 0.2), "Cd34Low", "misc"))
Idents(group4) <- "Cd27vs34"
de.markers <- FindMarkers(group4, ident.1 = group4$Cd27_highvslow == "Low", ident.2 = group4$Cd34_highvslow == "Low")
write.csv(de.markers, "Cd34_highvslow_markers_Adult.csv")

write.csv(gene_expression, "gene_expression.csv")

#testing

