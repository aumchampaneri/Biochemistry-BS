# Import packages
import scanpy as sc

# Load the kidney data
kidney = sc.read_h5ad("/Volumes/CHAMPANERI/Databases/MouseKidneyData.h5ad")

####
# Each run will REDO the PCA analysis and clustering
# Use a Jupyter Notebook for testing and debugging
####

'''
Define Dictionary of Genes for each set of plots

> 'Human_Census Gene Query.ipynb' has the gene names for the complement genes

all_genes: Dictionary of all genes to be plotted
    - Any genes not in the dataset are commented out
subset_genes: Dictionary of subset of genes to be plotted

'''

all_genes = {
    "C1d": "ENSMUSG00000000581",
    # "C1qa": "ENSMUSG00000036887",
    "C1qb": "ENSMUSG00000036905",
    "C1qbp": "ENSMUSG00000018446",
    # "C1qc": "ENSMUSG00000036896",
    # "C2_ENSMUSG00000024371": "ENSMUSG00000024371",
    "C2": "ENSMUSG00000024371",
    # "C3_ENSMUSG00000024164": "ENSMUSG00000024164",
    "C3": "ENSMUSG00000024164",
    "C3ar1": "ENSMUSG00000040552",
    "C4a": "ENSMUSG00000015451",
    # "C4b": "ENSMUSG00000073418",
    # "Hc": "ENSMUSG00000026874",
    # "C5": "ENSMUSG00000026874",
    "C5ar1": "ENSMUSG00000049130",
    "C5ar2": "ENSMUSG00000074361",
    # "C6_ENSMUSG00000022181": "ENSMUSG00000022181",
    # "C6": "ENSMUSG00000022181",
    # "C7_ENSMUSG00000079105": "ENSMUSG00000079105",
    # "C7": "ENSMUSG00000079105",
    "C8a": "ENSMUSG00000035031",
    # "C8b": "ENSMUSG00000029656",
    "C8g": "ENSMUSG00000015083",
    # "C9_ENSMUSG00000022149": "ENSMUSG00000022149",
    # "C9": "ENSMUSG00000022149",
    "Cd46": "ENSMUSG00000016493",
    "Cd55": "ENSMUSG00000026399",
    "Cd59a": "ENSMUSG00000032679",
    "Cfb": "ENSMUSG00000090231",
    # "Cfd": "ENSMUSG00000061780",
    "Cfh": "ENSMUSG00000026365",
    "Cfi": "ENSMUSG00000058952",
    "Cfp": "ENSMUSG00000001128",
    "Clu": "ENSMUSG00000022037",
    "Cr1l": "ENSMUSG00000016481",
    "Cr2": "ENSMUSG00000026616",
    "Serping1": "ENSMUSG00000023224",
    "Vcp": "ENSMUSG00000028452",
    # "Vtn": "ENSMUSG00000017344"
}

subset_genes = {
    "C3": "ENSMUSG00000024164",
    "C3ar1": "ENSMUSG00000040552",
    # "C5": "ENSMUSG00000026874",
    "C5ar1": "ENSMUSG00000049130",
    "C5ar2": "ENSMUSG00000074361",
    "Cfh": "ENSMUSG00000026365",
}

'''
PLOTTING SECTION:

Dotplots:
    - Cell Type
    - Tissue
    - Development Stage

Heatmaps:
    - Cell Type
    - Tissue

'''

# Define the dendrogram location in anndata
sc.tl.pca(kidney)
sc.pl.pca_variance_ratio(kidney, n_pcs=50, log=True)
sc.tl.dendrogram(kidney, groupby='cell_type')

# Plot the dotplots for all_genes
sc.pl.dotplot(
    kidney,
    all_genes,
    groupby='cell_type',
    title='Kidney Cell Type Complement Gene Expression',
    dendrogram=True,
    mean_only_expressed=True,
    log=True,
    colorbar_title='mean expression in group (log scale)',
    save='mouse-h5ad_all-gene_cell-type_dp.pdf'
)

sc.pl.dotplot(
    kidney,
    all_genes,
    groupby='tissue',
    title='Kidney Tissue Complement Gene Expression',
    dendrogram=True,
    mean_only_expressed=True,
    log=True,
    colorbar_title='mean expression in group (log scale)',
    save='mouse-h5ad_all-gene_tissue_dp.pdf'
)

# Plot the dotplots for subset_genes
sc.pl.dotplot(
    kidney,
    subset_genes,
    groupby='cell_type',
    title='Kidney Cell Type Complement Gene Expression',
    dendrogram=True,
    mean_only_expressed=True,
    log=True,
    colorbar_title='mean expression in group (log scale)',
    save='mouse-h5ad_subset-genes_cell-type_dp.pdf'
)

sc.pl.dotplot(
    kidney,
    subset_genes,
    groupby='tissue',
    title='Kidney Tissue Complement Gene Expression',
    dendrogram=True,
    mean_only_expressed=True,
    log=True,
    colorbar_title='mean expression in group (log scale)',
    save='mouse-h5ad_subset-genes_tissue_dp.pdf'
)


# Plot the heatmap for all_genes
sc.pl.heatmap(
    kidney,
    all_genes,
    groupby='cell_type',
    swap_axes=False,
    log=False,
    use_raw=False,
    figsize=(12,20),
    dendrogram=True,
    save='mouse-h5ad_all-genes_cell-type_hm.pdf'
)
sc.pl.heatmap(
    kidney,
    all_genes,
    groupby='tissue',
    swap_axes=False,
    log=False,
    use_raw=False,
    figsize=(12,20),
    dendrogram=True,
save='mouse-h5ad_all-genes_tissue_hm.pdf'
)

# Plot the heatmap for subset_genes
sc.pl.heatmap(
    kidney,
    subset_genes,
    groupby='cell_type',
    swap_axes=False,
    log=False,
    figsize=(3,22),
    dendrogram=True,
    use_raw=False,
    save='mouse-h5ad_subset-genes_cell-type_hm.pdf'
)
sc.pl.heatmap(
    kidney,
    subset_genes,
    groupby='tissue',
    swap_axes=False,
    log=False,
    figsize=(3,6),
    dendrogram=True,
    use_raw=False,
    save='mouse-h5ad_subset-genes_tissue_hm.pdf'
)