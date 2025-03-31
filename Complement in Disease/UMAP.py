import scanpy as sc
from BioMart_Query import gene_dict

# Load the processed data
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA_processed.h5ad")
