# Prepare the environment
import gget
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

gget.setup("cellxgene")

# Uncomment the following line to see the documentation
# help(gget.cellxgene)

'''
Query the cellxgene database for the genes of interest
Filters:
Genes: C1QA, C1QC, C1QB, FCN3, CFH, CD55, CR1, CD46, C1D, CFI, C7, C2, CFP, VCP, C5, FCN2, FCN1, CD59, C1S, C1R, C3AR1, C1QBP, CFD, C3, C5AR1, C5AR2, C6, C9, CR2, C4B, CFB, C4A
disease: normal
tissue_general: kidney

> Query Time: 3m 30s
'''

adata = gget.cellxgene(
    ensembl=True,
    verbose=True,
    gene=["ENSMUSG00000000581", "ENSMUSG00000036887", "ENSMUSG00000036905", "ENSMUSG00000018446", "ENSMUSG00000036896", "ENSMUSG00000024371", "ENSMUSG00000024164", "ENSMUSG00000040552", "ENSMUSG00000015451", "ENSMUSG00000073418", "ENSMUSG00000026874", "ENSMUSG00000049130", "ENSMUSG00000074361", "ENSMUSG00000022181", "ENSMUSG00000079105", "ENSMUSG00000035031", "ENSMUSG00000029656", "ENSMUSG00000015083", "ENSMUSG00000022149", "ENSMUSG00000016493", "ENSMUSG00000026399", "ENSMUSG00000032679", "ENSMUSG00000090231", "ENSMUSG00000061780", "ENSMUSG00000026365", "ENSMUSG00000058952", "ENSMUSG00000001128", "ENSMUSG00000022037", "ENSMUSG00000016481", "ENSMUSG00000026616", "ENSMUSG00000023224", "ENSMUSG00000028452", "ENSMUSG00000017344"],
    disease='normal',
    tissue_general='kidney',
    species='mus_musculus',
)
"""
Dictionary mapping gene names to their corresponding Ensembl IDs.

Keys:
    str: Gene names (e.g., 'C1D', 'C1QA', etc.)
Values:
    str: Ensembl IDs (e.g., 'ENSG00000197223', 'ENSG00000173372', etc.)
"""
all_genes = {
    "ENSMUSG00000000581": "C1d",
    "ENSMUSG00000036887": "C1qa",
    "ENSMUSG00000036905": "C1qb",
    "ENSMUSG00000018446": "C1qbp",
    "ENSMUSG00000036896": "C1qc",
    "ENSMUSG00000024371": "C2",
    "ENSMUSG00000024164": "C3",
    "ENSMUSG00000040552": "C3ar1",
    "ENSMUSG00000015451": "C4a",
    "ENSMUSG00000073418": "C4b",
    "ENSMUSG00000026874": "Hc",
    "ENSMUSG00000049130": "C5ar1",
    "ENSMUSG00000074361": "C5ar2",
    "ENSMUSG00000022181": "C6",
    "ENSMUSG00000079105": "C7",
    "ENSMUSG00000035031": "C8a",
    "ENSMUSG00000029656": "C8b",
    "ENSMUSG00000015083": "C8g",
    "ENSMUSG00000022149": "C9",
    "ENSMUSG00000016493": "Cd46",
    "ENSMUSG00000026399": "Cd55",
    "ENSMUSG00000032679": "Cd59a",
    "ENSMUSG00000090231": "Cfb",
    "ENSMUSG00000061780": "Cfd",
    "ENSMUSG00000026365": "Cfh",
    "ENSMUSG00000058952": "Cfi",
    "ENSMUSG00000001128": "Cfp",
    "ENSMUSG00000022037": "Clu",
    "ENSMUSG00000016481": "Cr1l",
    "ENSMUSG00000026616": "Cr2",
    "ENSMUSG00000023224": "Serping1",
    "ENSMUSG00000028452": "Vcp",
    "ENSMUSG00000017344": "Vtn"
}



"""
Dictionary mapping a subset of gene names to their corresponding Ensembl IDs.

Keys:
   str: Gene names (e.g., 'C3', 'C3AR1', etc.)

Values:
   str: Ensembl IDs (e.g., 'ENSG00000125730', 'ENSG00000171860', etc.)
"""
subset_genes = {
    "ENSMUSG00000024164": "C3",
    "ENSMUSG00000040552": "C3ar1",
    "ENSMUSG00000026874": "Hc",
    "ENSMUSG00000049130": "C5ar1",
    "ENSMUSG00000074361": "C5ar2",
    "ENSMUSG00000026365": "Cfh"
}


"""
Preprocess the kidney dataset by filtering genes, normalizing counts,
log-transforming the data, and scaling the data. Then, perform PCA
and plot an overview of the PCA results.

Steps:
1. Filter genes with at least 1 count.
2. Normalize the total counts to 1e6.
3. Log-transform the data.
4. Scale the data.
5. Perform PCA.
6. Plot an overview of the PCA results.
"""

sc.pp.filter_genes(adata, min_counts=1)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.scale(adata)

sc.tl.pca(adata)
sc.pl.pca_overview(adata)

"""
Check if 'feature_id' exists in `adata.var` and contains unique values.
If so, set `adata.var_names` to `adata.var['feature_id']`.
"""

if 'feature_id' in adata.var.columns and not adata.var['feature_id'].duplicated().any():
    adata.var_names = adata.var['feature_id']

"""
Filter the AnnData object to include only cells from the top 50 cell types
based on total gene expression. Verify the filtering by printing the 
original and filtered number of cells.
"""

# Check if cell_type annotation exists
if 'cell_type' not in adata.obs.columns:
    print("No cell_type annotation found in adata.obs")
    print("Available annotations:", adata.obs.columns)
else:
    # Calculate mean expression of our genes per cell type
    # First verify that the gene IDs are in the dataset
    valid_genes = [gene for gene in adata.var_names if gene in adata.var_names]
    # print(f"Found {len(valid_genes)} genes out of 28")

    # Group by cell type and calculate mean expression
    cell_type_expression = pd.DataFrame()
    for cell_type in adata.obs['cell_type'].unique():
        # Get indices for this cell type
        cells_idx = adata.obs['cell_type'] == cell_type
        # Calculate mean expression across these cells for each gene
        mean_expr = adata[cells_idx, valid_genes].X.mean(axis=0)
        # Sum the expression across all genes
        total_expr = np.sum(mean_expr)
        cell_type_expression.loc[cell_type, 'total_expression'] = total_expr

    # Sort by total expression and get top 50
    top_50_cell_types = cell_type_expression.sort_values('total_expression', ascending=False).head(50)

    # print("Top 50 cell types with highest gene expression:")
    # print(top_50_cell_types)

    # Filter the AnnData object to include only cells from these top 50 cell types
    top_50_cells = adata.obs['cell_type'].isin(top_50_cell_types.index)
    filtered_adata = adata[top_50_cells, :]

    # Verify the filtering
    # print(f"Original number of cells: {adata.n_obs}")
    # print(f"Number of cells in top 50 cell types: {filtered_adata.n_obs}")

"""
Generate a heatmap of gene expression for the top 50 cell types.

Steps:
1. Extract the expression data for the top 50 cell types based on total gene expression.
2. Filter the AnnData object to include only cells from these top 50 cell types.
3. Calculate the mean expression of each gene for the top 50 cell types.
4. Rename the columns from Ensembl IDs to gene names using the `all_genes` dictionary.
5. Create and display a heatmap with the transposed data.

The heatmap shows the gene expression levels for the top 50 cell types, with cell types on the x-axis and genes on the y-axis.
"""

# Extract the expression data for the top 50 cell types
top_50_cell_types = cell_type_expression.sort_values('total_expression', ascending=False).head(50)
top_50_cells = adata.obs['cell_type'].isin(top_50_cell_types.index)
filtered_adata = adata[top_50_cells, :]

# Calculate mean expression of each gene for the top 50 cell types
mean_expression = filtered_adata.to_df().groupby(filtered_adata.obs['cell_type']).mean()

# Rename the columns from Ensembl IDs to gene names
mean_expression.columns = [all_genes.get(col, col) for col in mean_expression.columns]

# Create a heatmap with transposed data -- ALL GENES
plt.figure(figsize=(16, 16))
sns.heatmap(mean_expression.T, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True, cbar_kws={"shrink":0.45})
plt.title('Gene Expression Heatmap for Top 50 Cell Types')
plt.xlabel('Cell Types')
plt.ylabel('Genes')
plt.tight_layout()
plt.savefig('Mouse_census_all-genes_cell-type_hm.pdf')

# Create a cluster with transposed data -- ALL GENES
plt.figure(figsize=(20, 35))
sns.clustermap(mean_expression, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True, cbar_pos=(0.01,0.8,0.0,0.3))
plt.tight_layout()
plt.savefig('Mouse_census_all-genes_cell-type_cm.pdf')

# Filter the mean expression data to only include the subset of genes
mean_expression = mean_expression[subset_genes.values()]

# Create a heatmap with transposed data -- SUBSET GENES
plt.figure(figsize=(16, 16))
sns.heatmap(mean_expression, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True, cbar_kws={"shrink":0.9})
plt.title('Gene Expression Heatmap for Top 50 Cell Types')
plt.xlabel('Cell Types')
plt.ylabel('Genes')
plt.tight_layout()
plt.savefig('Mouse_census_subset-genes_cell-type_hm.pdf')

# Create a clustermap with transposed data -- SUBSET GENES
plt.figure(figsize=(20, 35))
sns.clustermap(mean_expression, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True, cbar_pos=(0.01,0.8,0.0,0.3))
plt.tight_layout()
plt.savefig('Mouse_census_subset-genes_cell-type_cm.pdf')

"""
Calculate mean expression of our genes per tissue and generate heatmaps and clustermaps.

Steps:
1. Identify valid genes present in the dataset.
2. Initialize a DataFrame to store mean expression values for each tissue type.
3. For each tissue type:
   a. Extract cells belonging to the tissue type.
   b. Calculate mean expression of each gene across these cells.
   c. Handle sparse matrices if needed.
   d. Store the mean expression values in the DataFrame.
4. Convert all values to numeric and fill NAs with zeros.
5. Rename columns from Ensembl IDs to gene names using the `all_genes` dictionary.
6. Filter the tissues to remove the kidney data.
7. Generate and save heatmaps and clustermaps for all genes and the subset of genes.
"""

# Calculate mean expression of our genes per tissue
valid_genes = [gene for gene in all_genes.keys() if gene in adata.var_names]
print(f"Found {len(valid_genes)} out of {len(all_genes)} genes in the dataset")

# Initialize DataFrame with proper structure
tissue_types = adata.obs['tissue'].unique()
tissue_expression = pd.DataFrame(index=tissue_types, columns=valid_genes)

for tissue_type in tissue_types:
    # Get indices for this tissue
    cells_idx = adata.obs['tissue'] == tissue_type
    # Calculate mean expression across these cells for each gene
    subset = adata[cells_idx, valid_genes]

    # Handle sparse matrices if needed
    import scipy.sparse
    if scipy.sparse.issparse(subset.X):
        mean_expr = subset.X.mean(axis=0).A1  # Convert to 1D array
    else:
        mean_expr = np.array(subset.X.mean(axis=0))

    # Fill in the values explicitly
    for i, gene in enumerate(valid_genes):
        tissue_expression.loc[tissue_type, gene] = float(mean_expr[i])

# Convert all values to numeric and fill NAs
tissue_expression = tissue_expression.apply(pd.to_numeric, errors='coerce')
tissue_expression = tissue_expression.fillna(0)

# Rename columns from Ensembl IDs to gene names
tissue_expression.columns = [all_genes.get(col, col) for col in tissue_expression.columns]


# Filter the tissues to remove the kidney data
tissue_expression = tissue_expression.drop('kidney')

# Create a heatmap
plt.figure(figsize=(5, 10))
sns.heatmap(tissue_expression.T, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True)
plt.title('Gene Expression Across Tissues')
plt.ylabel('Genes')
plt.xlabel('Tissues')
plt.tight_layout()
plt.savefig('Mouse_census_all-genes_tissue_hm.pdf')

# Create a clustermap
plt.figure(figsize=(20, 20))
sns.clustermap(tissue_expression.T, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True)
plt.tight_layout()
plt.savefig('Mouse_census_all-genes_tissue_cm.pdf')

# Filter the mean expression data to only include the subset of genes
tissue_expression = tissue_expression[subset_genes.values()]

# Create a heatmap
plt.figure(figsize=(5, 5))
sns.heatmap(tissue_expression.T, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True)
plt.title('Gene Expression Across Tissues')
plt.ylabel('Genes')
plt.xlabel('Tissues')
plt.tight_layout()
plt.savefig('Mouse_census_subset-genes_tissue_hm.pdf')

# Create a clustermap
plt.figure(figsize=(20, 20))
sns.clustermap(tissue_expression.T, cmap='inferno', xticklabels=True, yticklabels=True, square=True, annot=False, robust=True)
plt.tight_layout()
plt.savefig('Mouse_census_subset-genes_tissue_cm.pdf')