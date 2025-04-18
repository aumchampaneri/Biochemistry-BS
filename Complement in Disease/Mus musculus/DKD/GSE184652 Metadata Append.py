import pandas as pd
import scanpy as sc
'''
Goal: Append metadata to the GSE184652 dataset for Mus musculus (DKD).
1. Load the existing AnnData object. (.h5ad file)
2. Load the metadata from the CSV file.
3. Merge the metadata into the AnnData object based on the cell barcodes.
4. Save the updated AnnData object.
'''

# Load the existing AnnData object
adata_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/adata.h5ad'
adata = sc.read(adata_path)

# Load the metadata from the CSV file
metadata_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/GSE184652_metadata.csv'
metadata = pd.read_csv(metadata_path)

# Ensure the 'cell_barcode' column in metadata matches the index of adata
metadata['cell_barcode'] = metadata['cell_barcode'].astype(str)

# Ensure the index of adata is also in string format
adata.obs_names = adata.obs_names.astype(str)

# Merge the metadata into the AnnData object
adata.obs = adata.obs.merge(metadata, left_index=True, right_on='cell_barcode', how='left')

# Drop the 'cell_barcode' column from obs if it exists
if 'cell_barcode' in adata.obs.columns:
    adata.obs.drop(columns=['cell_barcode'], inplace=True)

# Save the updated AnnData object
output_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/adata_with_metadata.h5ad'
adata.write(output_path)
print(f"Updated AnnData object saved to {output_path}")

# Print the updated AnnData object to verify the metadata has been added
print(adata)
# Print the first few rows of the obs DataFrame to verify the metadata
print(adata.obs.head())