import scanpy as sc

# Define the path to your data directory
data_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/count_matrix_combined_after_cellbender/'

# Read the data
adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', make_unique=True)

print(adata)
print(adata.obs)
print(adata.var)
print(adata.X)
print(adata.raw)
print(adata.uns)
print(adata.obsm)
print(adata.layers)

# Define the output file path
output_file = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/adata.h5ad'

# Save the AnnData object
adata.write(output_file)
