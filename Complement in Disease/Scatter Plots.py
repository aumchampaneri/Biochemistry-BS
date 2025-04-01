import scanpy as sc
import csv
import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import seaborn as sns

# Load the processed data -- Load only the dataset you want to plot
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Hs_Nor-CKD-AKF_scRNA_processed.h5ad") # Full dataset

# Load the gene dictionary from the csv file
gene_dict = {}
with open('complement_gene_dictionary.csv', newline='') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        gene_dict[row[0]] = row[1]

# Load the tissue type dictionary from the yaml file
with open("Tissue Type Dictionary.yaml", "r") as file:
    cell_type_group = yaml.safe_load(file)

# Map cell types to groups
adata.obs['cell_type_group'] = 'Other'
for group, cell_types in cell_type_group.items():
    adata.obs.loc[adata.obs['cell_type'].isin(cell_types), 'cell_type_group'] = group

#%%

def scatter_gene_expression(adata, disease, gene1, gene2, gene_dict, color_by='cell_type', save_path=None):
    """
    Creates a scatter plot plot of two genes' expression patterns colored by cell type or cell type group.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing single-cell RNA-seq data.
    disease : str
        The disease state to subset (e.g., 'AKI', 'CKD', 'Reference').
    gene1 : str
        Name of the first gene (x-axis).
    gene2 : str
        Name of the second gene (y-axis).
    gene_dict : dict
        Dictionary mapping gene names to Ensembl IDs for lookup in the expression matrix.
    color_by : str, optional
        Column name in adata.obs to use for coloring ('cell_type' or 'cell_type_group'). Default is 'cell_type'.
    save_path : str, optional
        Path to save the figure. If None, the figure is displayed but not saved.

    Returns
    -------
    None
    """
    # Subset data based on disease type
    adata = adata[adata.obs['diseasetype'] == disease, :].copy()

    # Fetch Ensembl ID from dictionary
    def get_ensembl_id(gene_name):
        return gene_dict.get(gene_name, gene_name)

    # Extract and flatten gene expression values
    def get_expression(ensembl_id):
        exp = adata[:, ensembl_id].X
        return exp.toarray().ravel() if scipy.sparse.issparse(exp) else exp.ravel()

    # Get gene expressions
    gene1_ensembl = get_ensembl_id(gene1)
    gene2_ensembl = get_ensembl_id(gene2)
    g1_vals = get_expression(gene1_ensembl)
    g2_vals = get_expression(gene2_ensembl)

    # Validate color_by column
    if color_by not in adata.obs.columns:
        raise ValueError(f"'{color_by}' column not found in adata.obs")

    # Get categories for coloring
    categories = adata.obs[color_by]
    unique_categories = list(categories.unique())

    # Use a discrete color palette (tab20 + additional colors if needed)
    palette = sns.color_palette("tab20", min(20, len(unique_categories)))
    if len(unique_categories) > 20:
        palette += sns.color_palette("tab20b", len(unique_categories) - 20)

    # Create scatter plot plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=g1_vals, y=g2_vals, hue=categories, palette=palette, alpha=0.7, edgecolor=None)

    # Set labels and title
    plt.xlabel(f'{gene1} Expression')
    plt.ylabel(f'{gene2} Expression')
    plt.title(f'{gene1} vs {gene2} Expression in {disease} Cells')

    # Adjust legend for clarity
    plt.legend(title=color_by.replace('_', ' ').title(), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, markerscale=1)

    # Save or show plot
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()
        


#%%

disease = 'AKI' # Choose the disease you want to plot (AKI, CKD, Reference)
color = 'cell_type_group' # Choose the color you want to plot (cell_type, cell_type_group)

# Plot scatter plot of C3 and CFH expression
scatter_gene_expression(adata, f"{disease}", "CFH", "C3", gene_dict, save_path=f"{disease}_{color}_scatter_CFHvC3.pdf")

# Plot scatter plot of CFH and C3ar1 expression
scatter_gene_expression(adata, f"{disease}", "CFH", "C3AR1", gene_dict, save_path=f"{disease}_{color}_scatter_CFHvC3ar1.pdf")

# Plot scatter plot of CFH and C5ar1 expression
scatter_gene_expression(adata, f"{disease}", "CFH", "C5AR1", gene_dict, save_path=f"{disease}_{color}_scatter_CFHvC5ar1.pdf")

# Plot scatter plot of CFH and C5ar2 expression
scatter_gene_expression(adata, f"{disease}", "CFH", "C5AR2", gene_dict, save_path=f"{disease}_{color}_scatter_CFHvC5ar2.pdf")

# Plot scatter plot of C3ar1 and C5ar1 expression
scatter_gene_expression(adata, f"{disease}", "C3AR1", "C5AR1", gene_dict, save_path=f"{disease}_{color}_scatter_C3ar1vC5ar1.pdf")

# Plot scatter plot of C3ar1 and C5ar2 expression
scatter_gene_expression(adata, f"{disease}", "C3AR1", "C5AR2", gene_dict, save_path=f"{disease}_{color}_scatter_C3ar1vC5ar2.pdf")

# Plot scatter plot of C5ar1 and C5ar2 expression
scatter_gene_expression(adata, f"{disease}", "C5AR1", "C5AR2", gene_dict, save_path=f"{disease}_{color}_scatter_C5ar1vC5ar2.pdf")