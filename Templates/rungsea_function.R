run_gsea <- function(dds, condition1, condition2, gmt_path, output_gsea, output_deg) {
  res1 <- DESeq2::results(dds, 
                          contrast = c("condition", condition1, condition2),
                          independentFiltering = TRUE,
                          alpha = 0.1)
  res1 <- subset(res1, padj < 0.1)
  rnk <- setNames(res1$log2FoldChange, rownames(res1))
  rnk <- rnk[!is.na(rnk)]
  gmt <- gmtPathways(gmt_path)
  gsea <- fgsea(pathways = gmt, stats = rnk, minSize = 15, maxSize = 500)
  gseadf <- as.data.frame(gsea)
  list_columns <- sapply(gseadf, is.list)
  gseadf[list_columns] <- lapply(gseadf[list_columns], function(col) sapply(col, paste, collapse = ";"))
  write.csv(gseadf, output_gsea)
  res1_filtered <- res1[complete.cases(res1), ]
  write.csv(res1_filtered, output_deg)
}