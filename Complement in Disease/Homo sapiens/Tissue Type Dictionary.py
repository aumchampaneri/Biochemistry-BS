import scanpy as sc
import yaml

# Load the processed data
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA_processed.h5ad")

# Define the dictionary for cell types
cell_type_groups = {
    'Kidney Cells': ['kidney distal convoluted tubule epithelial cell',
                     'epithelial cell of proximal tubule',
                     'kidney connecting tubule epithelial cell',
                     'kidney loop of Henle thick ascending limb epithelial cell',
                     'kidney collecting duct intercalated cell', 'kidney collecting duct principal cell',
                     'kidney loop of Henle thin descending limb epithelial cell',
                     'kidney loop of Henle thin ascending limb epithelial cell',
                     'kidney interstitial cell', 'kidney interstitial alternatively activated macrophage',
                     'podocyte', 'parietal epithelial cell'],

    'Myeloid Cells': ['monocyte',
                      'non-classical monocyte',
                      'mononuclear phagocyte',
                      'conventional dendritic cell',
                      'plasmacytoid dendritic cell, human'],

    'Lymphoid Cells': ['T cell',
                       'cytotoxic T cell',
                       'natural killer cell',
                       'mature NK T cell',
                       'B cell',
                       'plasma cell'],

    'Mast Cells': ['mast cell'],

    'Endothelial Cells': ['endothelial cell'],

    'Epithelial Cells': ['podocyte',
                         'parietal epithelial cell']
}

# Save the dictionary to a YAML file
with open("Tissue Type Dictionary.yaml", "w") as file:
    yaml.dump(cell_type_groups, file, default_flow_style=False, sort_keys=False)