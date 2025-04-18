import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("/Users/aumchampaneri/Databases/Triple/Human_Nor-CKD-AKF_scRNA_processed.h5ad")

# Returns a DataFrame with all Ensembl gene IDs and external gene names in the Homo sapiens genome
biomart_query = sc.queries.biomart_annotations("hsapiens", # Query for Homo sapiens
                                              ["ensembl_gene_id", "external_gene_name"],
                                               host='www.ensembl.org')

# List of Complement Gene Names to find the ensembl IDs for
# Comprehensive list of genes associated with the complement system and the complosome
complement_genes = {
    # Core complement components (classical, alternative, and lectin pathways)
    "core_complement": [
        "C1QA", "C1QB", "C1QC",  # C1q subunits (classical pathway initiation)
        "C1R", "C1S",  # C1 proteases
        "C2", "C3", "C4A", "C4B", "C5",  # Central complement cascade
        "C6", "C7", "C8A", "C8B", "C8G", "C9"  # Membrane attack complex (MAC)
    ],

    # Alternative pathway components
    "alternative_pathway": [
        "CFB",  # Factor B (forms alternative C3 convertase)
        "CFD",  # Factor D (activates Factor B)
        "CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5",  # Factor H and related proteins (regulators)
        "CFI"  # Factor I (C3b degradation)
    ],

    # Lectin pathway components
    "lectin_pathway": [
        "MBL2",  # Mannose-binding lectin
        "MASP1", "MASP2", "MASP3",  # MASP proteases
        "FCN1", "FCN2", "FCN3"  # Ficolins
    ],

    # Complement receptors
    "complement_receptors": [
        "CR1", "CR2", "CR3", "CR4",  # Receptors for complement components
        "C3AR1",  # C3a receptor
        "C5AR1", "C5AR2"  # C5a receptors
    ],

    # Complement regulatory proteins
    "regulatory_proteins": [
        "CD46", "CD55", "CD59",  # Membrane-bound regulators
        "SERPING1",  # C1 inhibitor
        "CLU", "VTN",  # Other regulators
        "PLG",  # Plasminogen (complement interaction)
        "CD35",  # Complement receptor 1 (CR1), important in immune complex clearance
        "THBD",  # Thrombomodulin, links complement and coagulation
        "VWF"  # von Willebrand Factor, affects complement-coagulation interactions
    ],

    # Complosome-specific genes (intracellular complement system)
    "complosome_core": [
        "C3", "C5",  # Intracellular complement pools
        "C3AR1", "C5AR1"  # Intracellular signaling receptors
    ],

    # Intracellular complement activation
    "intracellular_activation": [
        "CTSB", "CTSD", "CTSS"  # Cathepsins (intracellular C3 cleavage)
    ],

    # Metabolic and autophagy-linked components
    "metabolism_autophagy": [
        "MTOR", "ULK1", "ATG5", "ATG7", "LAMP1", "TFEB",  # Autophagy and metabolism regulators
        "SQSTM1", "BECN1", "GABARAP"  # Additional autophagy-related genes
    ],

    # Mitochondrial and cellular stress response
    "mitochondrial_response": [
        "VDAC1",  # Mitochondrial apoptosis regulator
        "NLRP3",  # Inflammasome activation
        "CASP1", "CASP8"  # Caspases (complement-apoptosis link)
    ],

    # Additional regulators of intracellular complement
    "intracellular_regulators": [
        "CFH", "CD46", "CD55", "CD59"  # Intracellular pools of regulators
    ],

    # Autoimmune and inflammatory disease-related genes
    "autoimmune_disease": [
        "TREX1",  # Links DNA damage and complement activation in lupus
        "IRF3", "IRF7",  # Interferon regulators (complement-immune system cross-talk)
        "STAT3",  # Key immune signaling protein (complement & inflammation)
        "TLR2", "TLR4"  # Toll-like receptors (synergy with complement)
    ],

    # Neuroinflammation and synaptic pruning genes
    "neuroinflammation": [
        "APOE",  # Apolipoprotein E, modulates complement in Alzheimer’s
        "HSP90AA1",  # Heat Shock Protein 90, involved in complement folding
        "C1QTNF6",  # Complement-related protein in neuronal pruning
        "GPNMB"  # Neurodegenerative complement-related protein
    ],

    # DNA repair & non-canonical complement functions
    "dna_repair": [
        "DDB1"  # Damage-Specific DNA Binding Protein 1 (links complement to DNA repair)
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