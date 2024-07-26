library(DESeq2)
library(tximport)
library(tidyverse)
library(rnaseqGene)
library(pheatmap)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MDS/GSE63569/")

#DESeq2  analysis ------------------------
#Load data into DESeq dataset
cts2 <- read.csv("GSE63569_raw_counts_GRCh38.p13_NCBI.csv") #file path to Counts csv document
cts2 <- cts2[, -1] #only necessary if bulk seq data has ensembl id before gene names
cts <- cts2[, -1]
table(duplicated(cts$Gene.name)) #check for duplicate gene names
rownames(cts) <- make.names(cts2$Gene.name, unique = TRUE) #forces gene names to become unique
rm(cts2)
coldata <- read.csv("coldata.csv") #make coldata file with sample names in one column and then condition in second column
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

#Pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Differential Expression
dds <- DESeq(dds)
rescmp <- results(dds) 

#write results as csv
write.csv(rescmp, "completeresults.csv")

#display as heatmap
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

colnames(assay(ntd)[select,]) <- str_sub(colnames(assay(ntd)[select,]), 1, -3)
rownames(df) <- colnames(assay(ntd)[select,])
# For displaying the plot
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

# For saving the plot
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename="heatmap.pdf")
