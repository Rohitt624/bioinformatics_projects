library(DESeq2)
library(tximport)
library(tidyverse)
library(rnaseqGene)
library(pheatmap)
setwd("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/Kanak Tet2 Ripk3/")

#DESeq2  analysis ------------------------
#Load data into DESeq dataset
cts2 <- read.csv("raw_counts modified corrected triplicates.csv") #file path to Counts csv document
cts2 <- cts2[, -1] #only necessary if bulk seq data has ensembl id before gene names
cts <- cts2[, -1]
table(duplicated(cts$Gene.name)) #check for duplicate gene names
rownames(cts) <- make.names(cts2$Gene.name, unique = TRUE) #forces gene names to become unique
rm(cts2)
cd <- read.csv("coldata.csv") #make coldata file with sample names in one column and then condition in second column
dds <- DESeqDataSetFromMatrix(countData = cts, colData = cd, design = ~ condition)


#Differential Expression
dds <- DESeq(dds)
rescmp <- results(dds) 
#write results as csv
write.csv(rescmp, "completeresults.csv")

# Comparison 1 of 3:  DKO-WT
res <- DESeq2::results(dds, 
                       contrast = c("condition", "DKO", "WT"),
                       independentFiltering = TRUE,
                       alpha = 0.05
)

# Filter for p-value < 0.05
res_filtered <- res[which(res$pvalue < 0.05), ]

write.csv(res_filtered, "dko_wt_filter.csv")

write.csv(res, "dko_wt.csv")



# Comparison 2 of 3:  Tet2KO-WT
res <- DESeq2::results(dds, 
                       contrast = c("condition", "Tet2KO", "WT"),
                       independentFiltering = TRUE,
                       alpha = 0.05
)

# Filter for p-value < 0.05
res_filtered <- res[which(res$pvalue < 0.05), ]

write.csv(res_filtered, "tet2_wt_filter.csv")

write.csv(res, "tet2_wt.csv")

# Comparison 3 of 3:  DKO-Tet2KO
res <- DESeq2::results(dds, 
                       contrast = c("condition", "DKO", "Tet2KO"),
                       independentFiltering = TRUE,
                       alpha = 0.05
)

# Filter for p-value < 0.05
res_filtered <- res[which(res$pvalue < 0.05), ]

write.csv(res_filtered, "dko_tet2ko_filter.csv")

write.csv(res, "dko_tet2ko.csv")



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
