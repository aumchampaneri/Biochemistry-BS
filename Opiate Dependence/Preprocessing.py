import scanpy as sc

adata = sc.read_h5ad('/Users/aumchampaneri/Databases/Opiate Dependence/665dc4a2-7f38-4d78-b632-f575e3dd550e.h5ad')

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
sc.write('/Users/aumchampaneri/Databases/Opiate Dependence/GSE233279_CxG_pp.h5ad', adata)