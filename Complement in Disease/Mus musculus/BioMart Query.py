import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Mm_MKA-pp.h5ad")

# Returns a DataFrame with all Ensembl gene IDs and external gene names in the Mus musculus genome
biomart_query = sc.queries.biomart_annotations("mmusculus", # Query for Mus musculus
                                              ["ensembl_gene_id", "external_gene_name"],
                                               host='www.ensembl.org')

#%%
# List of Complement Gene Names to find the ensembl IDs for
# Comprehensive list of genes associated with the complement system and the complosome
complement_genes = {
    # Core complement components (classical, alternative, and lectin pathways)
    "core_complement": [
        "C1qa", "C1qb", "C1qc",      # C1q subunits
        "C1r", "C1s", "C1rl"              # C1 proteases
        "C2", "C3", "C4", "C5", "Hc", #C5 is Hc
        "C6", "C7", "C8a", "C8b", "C8g", "C9"
    ],

    # Alternative pathway components
    "alternative_pathway": [
        "Cfb",         # Factor B
        "Cfd",         # Factor D
        "Cfp",         # Properdin (Factor P)
        "Cfh",         # Factor H
        # None,          # CFHR1 – no direct mouse ortholog
        # None,          # CFHR2 – no direct mouse ortholog
        # None,          # CFHR3 – no direct mouse ortholog
        # None,          # CFHR4 – no direct mouse ortholog
        # None,          # CFHR5 – no direct mouse ortholog
        "Cfi"          # Factor I
    ],

    # Lectin pathway components
    "lectin_pathway": [
        "Mbl2",        # Mannose-binding lectin (mouse Mbl2 is the main functional gene)
        "Masp1", "Masp2", "Masp3",  # MASP proteases
        "Fcna", "Fcnb",         # Ficolins: FCN1 -> Fcna, FCN2 -> Fcnb; FCN3 has no clear mouse ortholog
    ],

    # Complement receptors
    "complement_receptors": [
        "Cr1",   # CR1: mouse lacks a true CR1; the Cr2 gene produces both CR1/CR2 isoforms
        "Cr2",   # CR2
        "Itgam", # CR3 (CD11b, alpha chain)
        "Itgax", # CR4 (CD11c, alpha chain)
        "C3ar1", # C3a receptor
        "C5ar1", # C5a receptor 1
        "C5ar2"  # C5a receptor 2
    ],

    # Complement regulatory proteins
    "regulatory_proteins": [
        "Cd46",    # CD46 (note: mouse CD46 expression is tissue‐restricted)
        "Cd55",    # CD55
        "Cd59a",   # CD59 (mouse mainly expresses Cd59a)
        "Serping1",# SERPING1 (C1 inhibitor)
        "Clu",     # CLU (clusterin)
        "Vtn",     # VTN (vitronectin)
        "Plg",     # PLG (plasminogen)
        "Cr2",     # CD35 (human CR1/CD35; mouse uses Cr2 as the functional homolog)
        "Thbd",    # THBD (thrombomodulin)
        "Vwf",     # VWF (von Willebrand Factor)
        "C4bp"     # C4b-binding protein
    ],

    # Complosome-specific genes (intracellular complement system)
    "complosome_core": [
        "C3", "C5", "C3ar1", "C5ar1"
    ],

    # Intracellular complement activation
    "intracellular_activation": [
        "Ctsb", "Ctsd", "Ctss"  # Cathepsins B, D, S
    ],

    # Metabolic and autophagy-linked components
    "metabolism_autophagy": [
        "Mtor", "Ulk1", "Atg5", "Atg7", "Lamp1", "Tfeb",
        "Sqstm1", "Becn1", "Gabarap"
    ],

    # Mitochondrial and cellular stress response
    "mitochondrial_response": [
        "Vdac1", "Nlrp3", "Casp1", "Casp8"
    ],

    # Additional regulators of intracellular complement
    "intracellular_regulators": [
        "Cfh", "Cd46", "Cd55", "Cd59a"
    ],

    # Autoimmune and inflammatory disease-related genes
    "autoimmune_disease": [
        "Trex1", "Irf3", "Irf7", "Stat3", "Tlr2", "Tlr4"
    ],

    # Neuroinflammation and synaptic pruning genes
    "neuroinflammation": [
        "Apoe", "Hsp90aa1", "C1qtnf6", "Gpnmb"
    ],

    # DNA repair & non-canonical complement functions
    "dna_repair": [
        "Ddb1"
    ]
}


# Create a sorted one-dimensional list from the dictionary
all_complement_genes = sorted(set(gene for sublist in complement_genes.values() for gene in sublist))

# Find the complement genes in the biomart_query DataFrame
complement_genes_set = set(all_complement_genes)
biomart_genes_set = set(biomart_query["external_gene_name"].values)

# Find the intersection
intersection = complement_genes_set.intersection(biomart_genes_set)

# Filter the biomart_query DataFrame to include only the intersecting genes
filtered_biomart_query = biomart_query[biomart_query["external_gene_name"].isin(intersection)]

# Alphabetize the DataFrame by external_gene_name
filtered_biomart_query = filtered_biomart_query.sort_values(by="external_gene_name")

# Reset the index
filtered_biomart_query = filtered_biomart_query.reset_index(drop=True)

# Create a dictionary with gene names on the left and Ensembl codes on the right
gene_dict = dict(zip(filtered_biomart_query["external_gene_name"], filtered_biomart_query["ensembl_gene_id"]))

# Find the genes in the adata object
adata_genes = set(adata.var_names.values)

# Find the intersection with the gene_dict ensembl IDs
intersection_adata = set(gene_dict.values()).intersection(adata_genes)

# Filter the gene_dict to include only the intersecting genes
filtered_gene_dict = {k: v for k, v in gene_dict.items() if v in intersection_adata}

# Alphabetize the dictionary by gene names
filtered_gene_dict = dict(sorted(filtered_gene_dict.items()))

# Reset the index
filtered_gene_dict = {k: v for i, (k, v) in enumerate(filtered_gene_dict.items())}

# Print the dictionary
print(filtered_gene_dict)

# Save the dictionary to a CSV file
df = pd.DataFrame(list(filtered_gene_dict.items()), columns=["Gene Name", "Ensembl ID"])
df.to_csv("complement_gene_dictionary.csv", index=False)