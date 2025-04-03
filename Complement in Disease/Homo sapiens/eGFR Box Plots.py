# Import packages
import scanpy as sc
import csv
import yaml
import re
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

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
def plot_egfr_distribution(adata, disease_types, show_stats=False, save_path=None, figsize=(8, 6), stats_position=(0.48, 0.13)):
    """
    Plot the distribution of eGFR across different disease types.

    Parameters:
    - adata: AnnData object containing the data.
    - disease_types: List of disease types to include in the plot.
    - show_stats: Boolean indicating whether to show Chi-Square test results on the plot.
    - save_path: Path to save the plot. If None, the plot will be shown.
    - figsize: Tuple indicating the size of the plot (width, height).
    - stats_position: Tuple indicating the position of the stats results on the plot (x, y).

    Returns:
    - None
    """
    # Extract the relevant data and remove rows with unknown eGFR category and 'Reference' diseasetype
    eGFR_data = adata.obs[['diseasetype', 'eGFR']].dropna()
    eGFR_data = eGFR_data[(eGFR_data['eGFR'] != 'unknown') & (eGFR_data['diseasetype'] != 'Reference')]

    # Remove 'Reference' category from the 'diseasetype' column
    eGFR_data['diseasetype'] = eGFR_data['diseasetype'].cat.remove_categories(['Reference'])

    # Filter data by disease types
    if disease_types != 'all':
        eGFR_data = eGFR_data[eGFR_data['diseasetype'].isin(disease_types)]

    # Create a boxplot
    plt.figure(figsize=figsize)
    sns.boxplot(x='diseasetype', y='eGFR', data=eGFR_data)
    plt.title('eGFR Distribution Across Different Disease Types')
    plt.xlabel('Disease Type')
    plt.ylabel('eGFR')
    plt.xticks(rotation=45)

    # Reverse the y-axis
    plt.gca().invert_yaxis()

    if show_stats:
        # Create a contingency table
        contingency_table = pd.crosstab(eGFR_data['diseasetype'], eGFR_data['eGFR'])

        # Perform Chi-Square test
        chi2, p, dof, expected = stats.chi2_contingency(contingency_table)

        # Add Chi-Square test results to the plot
        plt.text(stats_position[0], stats_position[1], f'Chi-Square test statistic: {chi2:.2f}\n'
                                                       f'p-value: {p:.2e}\n'
                                                       f'Degrees of freedom: {dof}',
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=plt.gca().transAxes,
                 bbox=dict(facecolor='white', alpha=0.5))

    # Save or show plot
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()

#%%
# Plotting eGFR distribution for AKI and CKD
plot_egfr_distribution(adata, disease_types=['AKI', 'CKD'], show_stats=True,
                       save_path='egfr_boxplot_AKI-CKD_stats.pdf',
                       figsize=(8, 8), stats_position=(0.43, 0.12))
# plot_egfr_distribution(adata, disease_types=['AKI', 'CKD'], show_stats=False,
#                        save_path='egfr_boxplot_AKI-CKD.pdf',
#                        figsize=(8, 8), stats_position=(0.5, 0.15))