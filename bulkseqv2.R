library(DESeq2)
library(tximport)
library(tidyverse)
library(rnaseqGene)

ref <- read.csv("C:/Users/rthalla/OneDrive - Loyola University Chicago/Zhang Lab/RNASeq/MPC Bulk Seq Data/MPC_Counts_Ref.csv")
dds <- DESeqDataSetFromMatrix(countData = ref, colData = )