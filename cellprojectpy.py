import scanpy as sc
import cellproject as cp

#load the data
adata_ref = sc.datasets.pbmc3k_processed()
adata = sc.datasets.pbmc68k_reduced()

#make sure that the gene names are the same
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

#calculate pca
sc.pp.pca(adata_ref, n_comps=15)
sc.pp.neighbors(adata_ref, n_neighbors=10)

#umap of reference
umap_ref = cp.quick_umap(adata_ref)
sc.pl.umap(adata_ref, color='louvain')

#projection
cp.project_cells(adata, 
                 adata_ref,
                 obs_columns=['louvain'],
                 k=10,
                 scale_data=False,
                 use_vargenes=False,
                 umap_ref=umap_ref,
                 categorical_how='distribution')
                 
#plot projection
sc.pl.umap(adata, color=['ref_louvain', 'bulk_labels'])
comb = adata_ref.concatenate(adata)
sc.pl.umap(comb, color=['batch'])
