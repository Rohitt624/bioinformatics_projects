import scanpy as sc
import cellproject as cp
import scanpy.external as sce

adata_ref = sc.datasets.pbmc3k_processed()
adata = sc.datasets.pbmc68k_reduced()
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.pp.pca(adata_ref, n_comps = 15)
sc.pp.neighbors(adata_ref, n_neighbors = 10)
umap_ref = cp.quick_umap(adata_ref)
sc.pl.umap(adata_ref, color = 'louvain')

cp.project_cells(adata,
                 adata_ref,
                 obs_columns=['louvain'],
                 k=10,
                 scale_data=False,
                 use_vargenes=False,
                 umap_ref=umap_ref,
                 categorical_how='distribution')