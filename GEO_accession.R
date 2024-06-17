#download data directly from GEO

if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)
# gse <- getGEO('GSE60361') # does not work, the matrix is in a suppl file
geoFile <- "GSE60361_C1-3005-Expression.txt.gz"
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60361/suppl/GSE60361_C1-3005-Expression.txt.gz", destfile = geoFile)
