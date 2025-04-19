# Prepare the environment
import gget
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from numpy.lib.npyio import savez

gget.setup("cellxgene")

# Uncomment the following line to see the documentation
# help(gget.cellxgene)

'''
Query the cellxgene database for the genes of interest
Filters:
Genes: CFH, C3AR1, C3, C5AR1, C5AR2, C5
disease: normal
tissue_general: kidney

> Query Time: 3m 30s
'''

adata = gget.cellxgene(
    ensembl=True,
    verbose=True,
    # gene=['ENSG00000000971', 'ENSG00000171860', 'ENSG00000125730', 'ENSG00000197405', 'ENSG00000134830', 'ENSG00000106804'],
    disease='normal',
    tissue_general='kidney',
)

all_genes = {
"""
Dictionary mapping a subset of gene names to their corresponding Ensembl IDs.

Keys:
   str: Gene names (e.g., 'C3', 'C3AR1', etc.)

Values:
   str: Ensembl IDs (e.g., 'ENSG00000125730', 'ENSG00000171860', etc.)
"""
    "C3": "ENSG00000125730",
    "C3AR1": "ENSG00000171860",
    "C5": "ENSG00000106804",
    "C5AR1": "ENSG00000197405",
    "C5AR2": "ENSG00000134830",
    "CFH": "ENSG00000000971",
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

adata.write_h5ad(
    filename="/Users/aumchampaneri/Databases/Kidney gget/Homo sapiens/Human_Census_KidneyData.h5ad",
)

sc.pp.filter_genes(adata, min_counts=1)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.scale(adata)

sc.tl.pca(adata)

# Algorithically integrate multiple experiments
sc.external.pp.harmony_integrate(adata, key='dataset_id')

# Calculate UMAP
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

adata.write_h5ad(
    filename="/Users/aumchampaneri/Databases/Kidney gget/Homo sapiens/Human_Census_KidneyData_PP.h5ad",
)