#ensembl id for ref matrix
# Load the necessary library
library(biomaRt)
# Specify the dataset
ensembl <- useEnsembl(biomart = "genes")
listDatasets(ensembl)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
#Get gene names
ensembl_ids <- read.csv("MPC_Ensembl.csv", header = TRUE, stringsAsFactors = FALSE)
gene_ids <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  values = ensembl_ids,
  mart = mart
)

write.csv(gene_ids, "gene_ids.csv")
