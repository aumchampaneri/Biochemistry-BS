# Import packages
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import csv
import yaml
import re
import seaborn as sns
import scipy.cluster.hierarchy as sch

# Load the processed data
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Hs_Nor-CKD-AKF_scRNA_processed.h5ad")

# Load the gene dictionary from the csv file
gene_dict = {}
with open('complement_gene_dictionary.csv', newline='') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        gene_dict[row[0]] = row[1]

# Extract keys and values into separate lists
gene_dict_names = list(gene_dict.keys())
gene_dict_keys = list(gene_dict.values())

# Change the name of some entries in gene_dict_names to fix plotting errors
gene_dict_names = [re.sub(r'\bC2\b', 'C2_ENSG00000166278', name) for name in gene_dict_names]
gene_dict_names = [re.sub(r'\bC3\b', 'C3_ENSG00000125730', name) for name in gene_dict_names]
gene_dict_names = [re.sub(r'\bC6\b', 'C6_ENSG00000039537', name) for name in gene_dict_names]
gene_dict_names = [re.sub(r'\bC7\b', 'C7_ENSG00000112936', name) for name in gene_dict_names]
gene_dict_names = [re.sub(r'\bC9\b', 'C9_ENSG00000113600', name) for name in gene_dict_names]

# Load the tissue type dictionary from the yaml file
with open("Tissue Type Dictionary.yaml", "r") as file:
    cell_type_group = yaml.safe_load(file)

# Map cell types to groups
adata.obs['cell_type_group'] = 'Other'
for group, cell_types in cell_type_group.items():
    adata.obs.loc[adata.obs['cell_type'].isin(cell_types), 'cell_type_group'] = group

#%%

def differential_expression_heatmap(
        adata, gene_dict, groupby="cell_type", group1="Reference", group2="CKD",
        save_path=None, title_fontsize=16, tick_fontsize=12, annot_fontsize=10,
        colorbar_position=[0.85, 0.75, 0.03, 0.2], figsize=(15, 12), dendrogram_ratio=(0.05, 0.2),
        colormap="bwr", show_dendrogram=True  # Added parameter
):
    """
    Generates a heatmap showing log fold change (logFC) in gene expression
    between two groups across cell types.

    Parameters
    ----------
    adata : AnnData
        Single-cell RNA-seq data.
    gene_dict : dict
        Dictionary mapping Gene Names to Ensembl IDs.
    groupby : str, optional
        Column in `adata.obs` defining cell types (default: "cell_type").
    group1 : str, optional
        Reference group (default: "Reference").
    group2 : str, optional
        Experimental/disease group (default: "CKD").
    save_path : str, optional
        Path to save the heatmap figure (default: None).
    title_fontsize : int, optional
        Font size of the heatmap title (default: 16).
    tick_fontsize : int, optional
        Font size of the tick labels (default: 12).
    annot_fontsize : int, optional
        Font size of heatmap annotations (default: 10).
    colorbar_position : list, optional
        Position of the colorbar as [x, y, width, height] (default: [0.85, 0.75, 0.03, 0.2]).
    figsize : tuple, optional
        Size of the heatmap figure (default: (15, 12)).
    dendrogram_ratio : tuple, optional
        Ratio of space allocated to the dendrograms (default: (0.05, 0.2)).
    colormap : str, optional
        Colormap to use for the heatmap (default: "bwr").
    show_dendrogram : bool, optional
        Whether to show the dendrogram (True) or create a simple heatmap (False).

    Returns
    -------
    matplotlib.figure.Figure
        The heatmap figure.
    """
    # Convert Ensembl IDs to gene names
    ensembl_to_gene = {v: k for k, v in gene_dict.items()}
    valid_ensembl_ids = [gene_id for gene_id in ensembl_to_gene.keys() if gene_id in adata.var_names]

    if not valid_ensembl_ids:
        raise ValueError("None of the specified genes are found in the dataset.")

    gene_labels = [ensembl_to_gene[ensembl_id] for ensembl_id in valid_ensembl_ids]

    # Extract expression data for selected genes
    expression_data = adata[:, valid_ensembl_ids].to_df()
    expression_data = expression_data.join(adata.obs[[groupby, "diseasetype"]])

    # Normalize data to 1e4 and log transform
    expression_data[valid_ensembl_ids] = np.log1p(expression_data[valid_ensembl_ids] / 1e4)

    # Check if both groups exist in the dataset
    if group1 not in expression_data["diseasetype"].values or group2 not in expression_data["diseasetype"].values:
        raise ValueError(f"One or both groups ({group1}, {group2}) not found in the dataset.")

    # Compute mean expression per cell type for both groups
    group1_mean = expression_data.loc[expression_data["diseasetype"] == group1].groupby(groupby)[
        valid_ensembl_ids].mean()
    group2_mean = expression_data.loc[expression_data["diseasetype"] == group2].groupby(groupby)[
        valid_ensembl_ids].mean()

    # Use adaptive pseudocount
    min_nonzero = max(expression_data[valid_ensembl_ids].replace(0, np.nan).min().min(), 1e-6)
    pseudocount = min_nonzero

    # Compute log fold change with careful handling of edge cases
    logFC_data = np.log2((group2_mean + pseudocount) / (group1_mean + pseudocount))

    # Rename columns to gene names
    logFC_data.index.name = "Cell Type"
    logFC_data.columns = gene_labels

    if show_dendrogram:
        # Cluster genes using hierarchical clustering
        row_linkage = sch.linkage(logFC_data.T, method="ward")

        # Plot heatmap with clustering
        g = sns.clustermap(
            logFC_data.T,
            cmap=colormap, center=0, linewidths=0.5, annot=True, fmt=".2f",
            row_cluster=True, col_cluster=False, row_linkage=row_linkage,
            figsize=figsize, annot_kws={"size": annot_fontsize},
            xticklabels=True, yticklabels=True,
            dendrogram_ratio=dendrogram_ratio,
            cbar_pos=colorbar_position
        )

        # Set title with adjustable font size
        g.ax_heatmap.set_title(f"Log Fold Change (logFC) in Gene Expression: {group2} vs {group1}",
                               fontsize=title_fontsize)

        # Rotate x-axis labels and set font sizes
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=tick_fontsize)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=tick_fontsize)

        fig = g.fig
    else:
        # Create a simple heatmap without dendrogram
        fig, ax = plt.subplots(figsize=figsize)
        g = sns.heatmap(
            logFC_data.T,
            cmap=colormap, center=0, linewidths=0.5, annot=True, fmt=".2f",
            ax=ax, cbar=True,
            xticklabels=True, yticklabels=True,
            annot_kws={"size": annot_fontsize}
        )

        # Set title with adjustable font size
        ax.set_title(f"Log Fold Change (logFC) in Gene Expression: {group2} vs {group1}", fontsize=title_fontsize)

        # Rotate x-axis labels and set font sizes
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=tick_fontsize)
        plt.setp(ax.yaxis.get_majorticklabels(), fontsize=tick_fontsize)

    # Save or show plot
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
    else:
        plt.show()

    return fig

#%%
'''
Generate genes of interest lists for plotting
- Use the BioMart query for ideas

gene_dict_names -> all the complement genes

'''

complement_receptors = ['CR1', 'CR2', 'CR3', 'CR4', 'C3AR1', 'C5AR1', 'C5AR2']
complement_regulatory_proteins = ['CD46', 'CD55', 'CD59', 'SERPING1', 'CLU', 'VTN', 'PLG', 'CD35', 'THBD', 'VWF']
complement_alternative_pathway = ['CFB', 'CFD', 'CFH', 'CFHR1', 'CFHR2', 'CFHR3', 'CFHR4', 'CFHR5', 'CFI']
core_complement = ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C2_ENSG00000166278', 'C3_ENSG00000125730', 'C4A', 'C4B', 'C5', 'C6_ENSG00000039537', 'C7_ENSG00000112936', 'C8A', 'C8B', 'C8G', 'C9_ENSG00000113600']
lectin_pathway = ['MBL2', 'MASP1', 'MASP2', 'MASP3', 'FCN1', 'FCN2', 'FCN3']
complosome_core = ['C3_ENSG00000125730', 'C5', 'C3AR1', 'C5AR1']
intracellular_activation = ['C3AR1', 'C5AR1', 'C5AR2']
metabolism_autophagy = ['ATG5', 'ATG7', 'ATG12', 'ATG16L1', 'BECN1', 'LC3B', 'ULK1']
mitochondiral_response = ['NLRP3', 'CASP1', 'CASP4', 'CASP5', 'CASP8', 'CASP9', 'CASP12']

#%%

# group1 = "Reference" what is the control group
disease = 'CKD' # Choose the disease you want to plot (AKI, CKD, Reference)

# Change the groupby variable to the one you want to plot
subset = 'cell_type' # Choose the groupby variable you want to plot (cell_type, cell_type_group)

# change the gene list to the one you want to plot
gene_list_f = 'gene_dict_names' # Change this to match the gene list you are plotting for the filename
gene_list = gene_dict_names

# Fix the function call by passing gene_dict instead of gene_list
differential_expression_heatmap(
    adata, gene_dict, groupby=f"{subset}", group1="Reference", group2=f"{disease}",
    save_path=f"{disease}_{subset}_logFC-heatmap_{gene_list_f}.pdf", title_fontsize=14,
    tick_fontsize=8, annot_fontsize=6, colorbar_position=[0.98, 0.30, 0.01, 0.50],
    figsize=(12, 17), dendrogram_ratio=(0.05, 0.2), colormap="RdYlBu", show_dendrogram=True  # Use a diverging colormap
)