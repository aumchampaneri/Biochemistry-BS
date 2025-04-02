import scanpy as sc
import anndata as ad

adata_sc = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA.h5ad")
adata_sn = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_snRNA.h5ad")

#%%
adata = ad.concat([adata_sn, adata_sc], merge='same')
adata