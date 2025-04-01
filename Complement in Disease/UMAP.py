import scanpy as sc
import csv
import yaml

# Load the processed data
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA_processed.h5ad")

# Load the gene dictionary from the csv file
gene_dict = {}
with open('complement_gene_dictionary.csv', newline='') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        gene_dict[row[0]] = row[1]

# Load the tissue type dictionary from the yaml file
with open("Tissue Type Dictionary.yaml", "r") as file:
    loaded_cell_type_groups = yaml.safe_load(file)

#%%
'''
UMAP Plotting
This code will plot UMAP for the processed data and color the points by cell type and tissue type.
'''

# Plot UMAP colored by cell type
sc.pl.umap(adata,
           color='cell_type',
           save='_cell-type.pdf')

# Plot the UMAP colored by tissue type
sc.pl.umap(adata,
           color=['cell_type_group'],
           save='_cell-type-group.pdf')

#%%
'''
UMAP Plot – Comparing Gene Expression
This code will plot UMAP for the processed data and color the points by gene expression.
'''
