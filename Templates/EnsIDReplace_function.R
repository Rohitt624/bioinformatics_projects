EnsIDReplace <- function(input, output){
  rna <- as.sparse(read.csv(file = input, header = TRUE, row.names = 1))
  rnanames <- rownames(rna)
  ensembl <- useEnsembl(biomart = "genes")
  listDatasets(ensembl)
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  gene_ids <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    values = rnanames,
    mart = mart
  )
  
  id_list <- split(gene_ids[["mgi_symbol"]], gene_ids[["ensembl_gene_id"]])
  rnaname <- character(nrow(rna))
  for (i in seq_len(nrow(rna))) {
    ensembl_id <- rownames(rna)[i]
    if (ensembl_id %in% names(id_list)) {
      rnaname[i] <- id_list[[ensembl_id]]
    } else{
      rnaname[i] <- ensembl_id
    }
  }
  rnaname[rnaname == ""] <- NA
  
  rownames(rna) <- rnaname
  rm(ensembl, gene_ids, id_list, mart, ensembl_id, i, rnaname)
  write.csv(rna, output)
}