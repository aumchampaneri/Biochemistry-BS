# Import packages
import scanpy as sc

# Load the kidney data
kidney = sc.read_h5ad("/Volumes/CHAMPANERI/Databases/KidneyData.h5ad")

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
    'C1D': 'ENSG00000197223',
    'C1QA': 'ENSG00000173372',
    'C1QB': 'ENSG00000173369',
    'C1QBP': 'ENSG00000108561',
    'C1QC': 'ENSG00000159189',
    # 'C1R': 'ENSG00000159403',
    'C1S': 'ENSG00000182326',
    # 'C2': 'ENSG00000166278',
    'C3': 'ENSG00000125730',
    'C3AR1': 'ENSG00000171860',
    # 'C4A': 'ENSG00000244731',
    # 'C4B': 'ENSG00000224389',
    'C5': 'ENSG00000106804',
    'C5AR1': 'ENSG00000197405',
    'C5AR2': 'ENSG00000134830',
    'C6': 'ENSG00000039537',
    'C7': 'ENSG00000112936',
    'C9': 'ENSG00000113600',
    'CD46': 'ENSG00000117335',
    'CD55': 'ENSG00000196352',
    'CD59': 'ENSG00000085063',
    # 'CFB': 'ENSG00000243649',
    # 'CFD': 'ENSG00000197766',
    'CFH': 'ENSG00000000971',
    'CFI': 'ENSG00000205403',
    'CFP': 'ENSG00000126759',
    'CR1': 'ENSG00000203710',
    'CR2': 'ENSG00000117322',
    'FCN1': 'ENSG00000085265',
    'FCN2': 'ENSG00000160339',
    'FCN3': 'ENSG00000142748',
    'VCP': 'ENSG00000165280'
}

subset_genes = {
    "C3": "ENSG00000125730",
    "C3AR1": "ENSG00000171860",
    "C5": "ENSG00000106804",
    "C5AR1": "ENSG00000197405",
    "C5AR2": "ENSG00000134830",
    "CFH": "ENSG00000000971",
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
    save='human-h5ad_all-gene_cell-type_dp.pdf'
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
    save='human-h5ad_all-gene_tissue_dp.pdf'
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
    save='human-h5ad_subset-genes_cell-type_dp.pdf'
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
    save='human-h5ad_subset-genes_tissue_dp.pdf'
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
    save='human-h5ad_all-genes_cell-type_hm.pdf'
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
save='human-h5ad_all-genes_tissue_hm.pdf'
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
    save='human-h5ad_subset-genes_cell-type_hm.pdf'
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
    save='human-h5ad_subset-genes_tissue_hm.pdf'
)