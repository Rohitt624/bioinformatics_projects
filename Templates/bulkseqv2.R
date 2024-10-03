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

#Differential Expression
rlog <- rlog(dds, blind=FALSE)
dds <- DESeq(dds)
rescmp <- results(dds) 

# Comparison 1 of 3:  DKO-Tet2KO
res <- DESeq2::results(dds, 
                       contrast = c("groups", "DKO", "Tet2KO"),
                       independentFiltering = TRUE,
                       alpha = alpha
)

#write results as csv
write.csv(rescmp, "completeresults.csv")

#filter for p-value
res_filtered1 <- res1[which(res1$padj < 0.1 & !is.na(res1$padj)), ]

write.csv(res_filtered1, "C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/Kanak Tet2 Ripk3/dko_wt_filter.csv")
