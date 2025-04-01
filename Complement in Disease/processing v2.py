import scanpy as sc
import anndata as ad

adata_sc = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA.h5ad")
adata_sn = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_snRNA.h5ad")

#%%

# Ensure both datasets have the same columns in .obs
common_obs_columns = set(adata_sc.obs.columns).intersection(set(adata_sn.obs.columns))
adata_sc.obs = adata_sc.obs[common_obs_columns]
adata_sn.obs = adata_sn.obs[common_obs_columns]

# Concatenate the datasets
adata = ad.concat(
    [adata_sc, adata_sn],
    join='outer'  # Keep all genes from both datasets
)

#%%
'''
Data is from CellxGene. Can expect it has been preprocessed.
Out of a abundance of caution, we will reprocess the data to ensure it is clean and ready for analysis.
1. Remove cells with less than 500 genes and genes expressed in less than 3 cells
2. Normalize the data to counts per 10k cells
3. Log transform the data
'''

# Remove cells with less than 500 genes and genes expressed in less than 3 cells
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize the data to counts per 10k cells
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform the data
sc.pp.log1p(adata)

# Save the processed data
sc.write("/Users/aumchampaneri/Databases/Triple/Hs_Nor-CKD-AKF_scRNA_processed.h5ad", adata)