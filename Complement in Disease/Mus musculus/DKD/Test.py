import scanpy as sc
import pandas as pd
import anndata as ad

# List of file paths to .h5 files
h5_files = [
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594468_E3019_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594469_A3020_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594470_F3021_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594471_G3022_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594472_H3023_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594473_N3024_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594474_B3025_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594475_A3026_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594476_B3027_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594477_C3028_raw_feature_bc_matrix.h5',
    '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/GSM5594478_A3029_raw_feature_bc_matrix.h5'
]

# Read metadata
metadata = pd.read_csv('/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_metadata.csv', index_col=0)

# Create list for AnnData objects
adata_list = []

# Load each h5 file using scanpy
for i, h5_file in enumerate(h5_files):
    adata = sc.read_10x_h5(h5_file)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Optional: You can try matching metadata by file name or barcode prefixes if needed
    # For now, we just append it
    adata_list.append(adata)

# Combine all AnnData objects
combined_adata = ad.concat(adata_list, label='experiment', keys=[f'experiment_{i + 1}' for i in range(len(h5_files))])

# Save combined data
combined_adata.write('/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW/combined_experiments.h5ad')

print("✅ Combined .h5 file data saved as .h5ad!")