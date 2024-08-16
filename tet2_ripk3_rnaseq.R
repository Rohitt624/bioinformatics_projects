library(DESeq2)
library(tximport)
library(tidyverse)
library(rnaseqGene)
library(pheatmap)
setwd("C:/Users/rohit/OneDrive - Loyola University Chicago/Zhang Lab/Kanak Tet2 Ripk3/")

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
#------------------------------------------------------------------------------------------

# Comparison 1 of 3:  DKO-WT
res1 <- DESeq2::results(dds, 
                       contrast = c("condition", "DKO", "WT"),
                       independentFiltering = TRUE,
)

# Filter for p-value 
res_filtered1 <- res1[which(res1$padj < 0.1), ]

write.csv(res_filtered1, "dko_wt_filter.csv")

write.csv(res1, "dko_wt.csv")



# Comparison 2 of 3:  Tet2KO-WT
res2 <- DESeq2::results(dds, 
                       contrast = c("condition", "Tet2KO", "WT"),
                       independentFiltering = TRUE,
)

# Filter for p-value 
res_filtered2 <- res2[which(res2$padj < 0.1), ]

write.csv(res_filtered2, "tet2_wt_filter.csv")

write.csv(res2, "tet2_wt.csv")

# Comparison 3 of 3:  DKO-Tet2KO
res3 <- DESeq2::results(dds, 
                       contrast = c("condition", "DKO", "Tet2KO"),
                       independentFiltering = TRUE,
)

# Filter for p-value
res_filtered3 <- res3[which(res3$padj < 0.1), ]

write.csv(res_filtered3, "dko_tet2ko_filter.csv")

write.csv(res3, "dko_tet2ko.csv")
#-------------------------------------------------------------------------------------


#display top 20 as heatmap
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
#------------------------------------------------------------------------------------------

#Display selected Genes
# Assuming 'dds' is your DESeqDataSet object and 'res' is your results object from DESeq2
# Extract normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# List of specific genes
genes_of_interest <- c("Rgs1", "Sh2d1a", "Crip1", "S100a13", "Gabbr1", "Dtnb", "Rsad2", "Cxcl9", 
                       "Kenj10", "Gstm5", "C1qb", "Slc40a1", "Slc43a1", "C1qc", "Veam1", "Marchf1", 
                       "Pla2g2d", "C1qa", "Slc25a37", "Ckap4", "Rgl1", "Mafb", "Pou2af1", "Ighm", 
                       "Il18", "Slco2b1", "Marcks", "Itga9", "Tns3", "Abcc3", "Wdfy3", "Arhgap23")

# Define the desired column order
desired_order <- c("WT_1", "WT_2", "WT_3", "Tet2KO_1", "Tet2KO_2", "Tet2KO_3", "DKO_1", "DKO_2", "DKO_3")

# Reorder the columns
subset_norm_counts <- subset_norm_counts[, desired_order]

# Extract column data (sample information) and reorder it
col_data <- colData(dds)[desired_order, ]

# Define color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define annotation colors for the conditions
annotation_colors <- list(
  group = c(WT = "blue", Tet2KO = "orange", DKO = "purple")
)

# Create a heatmap
pheatmap(subset_norm_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         cluster_cols = FALSE, 
         clustering_method = "complete", 
         color = color_palette,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap of AITL Genes")

