import scanpy as sc
import pandas as pd

'''
Goal: Use the GSE184653 dataset for Mus musculus (DKD) to filter the cells based on metadata.

1. Load the existing AnnData object (.h5ad file). -> Metadata is already appended in the previous step.
2. Filter the cells based on the metadata criteria.
    - "db/db" = Diabetic Mice
    - "db/m" = Normal Mice
    - Should be contained in 'group2' column of the metadata.
3. Save the filtered AnnData object.
'''

# Load the existing AnnData object with metadata
adata_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/adata_with_metadata.h5ad'
adata = sc.read(adata_path)

# Print the AnnData object to verify its structure
print(adata)
print(adata.obs.head())

# Filter the cells based on the 'group2' column in the metadata
filtered_adata = adata[adata.obs['group2'].isin(['db/db', 'db/m'])].copy()

# Print the number of cells before and after filtering
print(f"Number of cells before filtering: {adata.n_obs}")
print(f"Number of cells after filtering: {filtered_adata.n_obs}")

# Save the filtered AnnData object
output_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/filtered_adata.h5ad'
filtered_adata.write(output_path)
print(f"Filtered AnnData object saved to {output_path}")