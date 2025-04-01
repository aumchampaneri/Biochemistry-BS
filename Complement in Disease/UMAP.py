import scanpy as sc
import csv
import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse

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

The function `plot_gene_expression` takes two gene keys and a dictionary mapping Ensembl gene IDs to gene names.
It generates a UMAP plot with a custom colormap based on the expression levels of the two genes.
The function also creates a legend showing the expression levels of both genes.
The UMAP plot is saved to a specified path if provided.

Get the gene key from 'complement_gene_dictionary.csv' file
'''


def plot_gene_expression(adata, gene1, gene2, gene_dict, save_path=None):
    """
    Plots UMAP visualization of two gene expressions with a custom colormap.

    Parameters:
        adata (AnnData): The annotated data matrix.
        gene1 (str): First gene key (e.g., "ENSG00000125730").
        gene2 (str): Second gene key (e.g., "ENSG00000000971").
        gene_dict (dict): A dictionary mapping Ensembl gene IDs to gene names.
        save_path (str): Path to save the figure (optional). If None, the plot won't be saved.
    """

    def get_expression(key):
        exp = adata[:, key].X
        return exp.toarray().flatten() if scipy.sparse.issparse(exp) else exp.flatten()

    # Extract expression data
    g1_vals = get_expression(gene1)
    g2_vals = get_expression(gene2)

    # Determine the data range for each gene
    g1_min, g1_max = g1_vals.min(), g1_vals.max()
    g2_min, g2_max = g2_vals.min(), g2_vals.max()

    # Generate five evenly spaced tick values
    g1_ticks = np.linspace(g1_min, g1_max, 5)
    g2_ticks = np.linspace(g2_min, g2_max, 5)

    # Rescale the raw data to [0,1] for color mapping
    X_scatter = (g1_vals - g1_min) / (g1_max - g1_min)
    Y_scatter = (g2_vals - g2_min) / (g2_max - g2_min)

    # Define the white basis as light grey
    white_val = 0.9

    # Compute weight factors
    w = (1 - X_scatter) * (1 - Y_scatter)
    r = X_scatter * (1 - Y_scatter)
    b = (1 - X_scatter) * Y_scatter
    p = X_scatter * Y_scatter

    # Compute scatter plot colors
    scatter_R = white_val * w + 1 * r + 0 * b + 1 * p
    scatter_G = white_val * w + 0 * r + 0 * b + 0 * p
    scatter_B = white_val * w + 0 * r + 1 * b + 1 * p

    scatter_colors = np.column_stack([scatter_R, scatter_G, scatter_B])
    scatter_colors = np.clip(scatter_colors, 0, 1)

    # Create figure and axes
    fig = plt.figure(figsize=(10, 10), facecolor='white')
    ax_umap = fig.add_subplot(111)
    ax_legend = fig.add_axes([0.85, 0.83, 0.10, 0.10])
    ax_umap.set_facecolor('white')

    # Plot UMAP scatter
    ax_umap.scatter(
        adata.obsm['X_umap'][:, 0],
        adata.obsm['X_umap'][:, 1],
        c=scatter_colors,
        s=0.75,
        alpha=1.0
    )

    # Get gene names from the gene_dict
    gene1_name = gene_dict.get(gene1, gene1)  # Use gene ID if name is not found
    gene2_name = gene_dict.get(gene2, gene2)  # Use gene ID if name is not found

    # Update plot title and labels
    ax_umap.set_title(f'{gene1_name} and {gene2_name} Expression in Cells')
    ax_umap.set_xlabel('UMAP 1')
    ax_umap.set_ylabel('UMAP 2')
    ax_umap.grid(False)
    ax_umap.tick_params(axis='both', which='both', length=0)
    ax_umap.set_xticklabels([])
    ax_umap.set_yticklabels([])

    for spine in ax_umap.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.2)

    # ---- Build the legend colormap ----
    n_bins = 100
    x_lin = np.linspace(g1_min, g1_max, n_bins)
    y_lin = np.linspace(g2_min, g2_max, n_bins)
    X_leg, Y_leg = np.meshgrid(x_lin, y_lin)

    # Rescale legend meshgrid values to [0,1]
    X_leg_scaled = (X_leg - g1_min) / (g1_max - g1_min)
    Y_leg_scaled = (Y_leg - g2_min) / (g2_max - g2_min)

    # Compute weight factors for the legend
    w_leg = (1 - X_leg_scaled) * (1 - Y_leg_scaled)
    r_leg = X_leg_scaled * (1 - Y_leg_scaled)
    b_leg = (1 - X_leg_scaled) * Y_leg_scaled
    p_leg = X_leg_scaled * Y_leg_scaled

    # Compute legend color mapping
    legend_R = white_val * w_leg + 1 * r_leg + 0 * b_leg + 1 * p_leg
    legend_G = white_val * w_leg + 0 * r_leg + 0 * b_leg + 0 * p_leg
    legend_B = white_val * w_leg + 0 * r_leg + 1 * b_leg + 1 * p_leg

    legend_color = np.stack([legend_R, legend_G, legend_B], axis=2)
    legend_color = np.clip(legend_color, 0, 1)

    ax_legend.imshow(legend_color, origin='lower', extent=[g1_min, g1_max, g2_min, g2_max])
    ax_legend.set_xlabel(f'{gene1_name} Expression', fontsize=8)
    ax_legend.set_ylabel(f'{gene2_name} Expression', fontsize=8)
    ax_legend.set_title('Expression Level Legend', fontsize=8)
    ax_legend.grid(False)

    # Apply the new tick values
    ax_legend.set_xticks(g1_ticks)
    ax_legend.set_yticks(g2_ticks)
    ax_legend.set_xticklabels([f"{val:.2f}" for val in g1_ticks], fontsize=6)
    ax_legend.set_yticklabels([f"{val:.2f}" for val in g2_ticks], fontsize=6)

    for spine in ax_legend.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1)

    # Save the figure if save_path is provided
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')

    # Show the plot
    plt.tight_layout()
    plt.show()

# Plot UMAP of C3 and CFH expression
plot_gene_expression(adata, "ENSG00000125730", "ENSG00000000971", gene_dict, save_path="umap_C3vCFH.pdf")

# Plot UMAP of CFH and C3ar1 expression
plot_gene_expression(adata, "ENSG00000000971", "ENSG00000171860", gene_dict, save_path="umap_CFHvC3ar1.pdf")

# Plot UMAP of CFH and C5ar1 expression
plot_gene_expression(adata, "ENSG00000000971", "ENSG00000197405", gene_dict, save_path="umap_CFHvC5ar1.pdf")

# Plot UMAP of CFH and C5ar2 expression
plot_gene_expression(adata, "ENSG00000000971", "ENSG00000134830", gene_dict, save_path="umap_CFHvC5ar2.pdf")

# Plot UMAP of C3ar1 and C5ar1 expression
plot_gene_expression(adata, "ENSG00000171860", "ENSG00000197405", gene_dict, save_path="umap_C3ar1vC5ar1.pdf")

# Plot UMAP of C3ar1 and C5ar2 expression
plot_gene_expression(adata, "ENSG00000171860", "ENSG00000134830", gene_dict, save_path="umap_C3ar1vC5ar2.pdf")

# Plot UMAP of C5ar1 and C5ar2 expression
plot_gene_expression(adata, "ENSG00000197405", "ENSG00000134830", gene_dict, save_path="umap_C5ar1vC5ar2.pdf")