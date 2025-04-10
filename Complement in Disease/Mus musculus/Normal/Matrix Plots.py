import scanpy as sc
import csv
import yaml
import re

# Load the processed data -- Load only the dataset you want to plot
adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Mm_MKA-pp.h5ad") # Normal dataset

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
gene_dict_names = [re.sub(r'\bC2\b', 'C2_ENSMUSG00000024371', name) for name in gene_dict_names]
gene_dict_names = [re.sub(r'\bC3\b', 'C3_ENSMUSG00000024164', name) for name in gene_dict_names]

# Load the tissue type dictionary from the yaml file
with open("Tissue Type Dictionary.yaml", "r") as file:
    cell_type_group = yaml.safe_load(file)

# Map cell types to groups
adata.obs['cell_type_group'] = 'Other'
for group, cell_types in cell_type_group.items():
    adata.obs.loc[adata.obs['cell_type'].isin(cell_types), 'cell_type_group'] = group

#%%
'''
Generate genes of interest lists for plotting
- Use the BioMart query for ideas

gene_dict_names -> all the complement genes

!!! Currently lists below are human genes adjust for mouse gene names
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
# Denedrogram == .add_dendrogram(show=True, size=0.9)

# Subset data based on disease type
disease = 'PATO:0000461' # Only option for this dataset
adata = adata[adata.obs['disease_ontology_term_id'] == disease, :].copy()

# Change the groupby variable to the one you want to plot
subset = 'cell_type_group' # Choose the groupby variable you want to plot (cell_type, cell_type_group)

# change the gene list to the one you want to plot
gene_list_f = 'gene_dict_names' # Change this to match the gene list you are plotting for the filename
gene_list = gene_dict_names


sc.pl.MatrixPlot(adata,
                 gene_list,
                 groupby=f'{subset}',
                 gene_symbols='feature_name',
                 use_raw=False,
                 log=False,
                 # title='Complement Genes Expression in Kidney Cells',
                 ).style(cmap='plasma', edge_color='none').savefig(f'Mm-Normal_{subset}_matrix-nd_{gene_list_f}.pdf')

# Change the groupby variable to the one you want to plot
subset = 'cell_type' # Choose the groupby variable you want to plot (cell_type, cell_type_group)

sc.pl.MatrixPlot(adata,
                 gene_list,
                 groupby=f'{subset}',
                 gene_symbols='feature_name',
                 use_raw=False,
                 log=False,
                 # title='Complement Genes Expression in Kidney Cells',
                 ).style(cmap='plasma', edge_color='none').savefig(f'Mm-Normal_{subset}_matrix-nd_{gene_list_f}.pdf')