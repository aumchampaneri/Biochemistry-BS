import scanpy as sc
import yaml

# Load the processed data
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Mm_MKA-pp.h5ad")

# Define the dictionary for cell types
cell_type_groups = {
    'Kidney Cells': ['kidney distal convoluted tubule epithelial cell',
                     'epithelial cell of proximal tubule',
                     'kidney connecting tubule epithelial cell',
                     'kidney loop of Henle cortical thick ascending limb epithelial cell',
                     'kidney loop of Henle medullary thick ascending limb epithelial cell',
                     'kidney loop of Henle thin descending limb epithelial cell',
                     'kidney loop of Henle thin ascending limb epithelial cell',
                     'kidney collecting duct principal cell',
                     'renal alpha-intercalated cell',
                     'renal beta-intercalated cell',
                     'kidney collecting duct epithelial cell',
                     'macula densa epithelial cell',
                     'kidney cortex tubule cell',
                     'kidney glomerular epithelial cell',
                     'mesangial cell',
                     'podocyte'],

    'Myeloid Cells': ['macrophage',
                      'neutrophil',
                      'dendritic cell'],

    'Lymphoid Cells': ['T cell',
                       'B cell',
                       'natural killer cell'],

    # 'Mast Cells': [],

    'Endothelial Cells': ['endothelial cell',
                          'glomerular endothelial cell',
                          'kidney afferent arteriole endothelial cell',
                          'kidney efferent arteriole endothelial cell'],

    'Perivascular Cells': ['pericyte',
                           'vasa recta descending limb cell'],

    'Fibroblasts': ['fibroblast'],

    'Unknown/Unclassified': ['unknown']
}


# Save the dictionary to a YAML file
with open("Tissue Type Dictionary.yaml", "w") as file:
    yaml.dump(cell_type_groups, file, default_flow_style=False, sort_keys=False)