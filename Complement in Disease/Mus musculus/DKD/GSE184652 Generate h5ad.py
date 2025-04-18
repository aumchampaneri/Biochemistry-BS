import scanpy as sc

# Define the path to your data directory (Specifically the flder containing the 10x matrix (2.tsv files and 1 .mtx file))
data_path = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/count_matrix_combined_after_cellbender/'

# Read the data
adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', make_unique=True)

# Print the AnnData object to verify its structure
print(adata)
print(adata.obs)
print(adata.var)


# Define the output file path
output_file = '/Users/aumchampaneri/Databases/Mus musculus/DKD/GSE184652/adata.h5ad'

# Save the AnnData object
adata.write(output_file)
