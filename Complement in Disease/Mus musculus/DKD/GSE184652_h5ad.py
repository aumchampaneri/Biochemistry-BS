# load_and_merge_gse184652.py

import scanpy as sc
import os

# === CONFIG ===
DATA_PATH = '/Users/aumchampaneri/Databases/Mm DKD Two/GSE184652_RAW'  # <-- UPDATE this path
OUTPUT_FILE = "GSE184652_subset_raw.h5ad"

# Sample IDs to include
samples = [
    "GSM5594468", "GSM5594469", "GSM5594470", "GSM5594471", "GSM5594472",  # Grp1
    "GSM5594473", "GSM5594474", "GSM5594475", "GSM5594476", "GSM5594477", "GSM5594478"  # Grp2
]

# === LOAD + CLEAN ===
adatas = []

for sample in samples:
    matched = [f for f in os.listdir(DATA_PATH) if sample in f and f.endswith(".h5")]
    if not matched:
        raise FileNotFoundError(f"No .h5 file found for {sample} in {DATA_PATH}")
    file_path = os.path.join(DATA_PATH, matched[0])
    ad = sc.read_10x_h5(file_path)
    ad.var_names_make_unique()
    ad.obs["sample"] = sample
    adatas.append(ad.copy())  # Make sure to isolate the object

# Ensure all gene sets are aligned across datasets
for i, ad in enumerate(adatas):
    ad.var_names_make_unique()
    ad.var.index = ad.var.index.astype(str)  # Force str dtype for safe merge

# === CONCATENATE ===
adata = adatas[0].concatenate(adatas[1:], batch_key="sample", batch_categories=samples)

# Check how many gene names are duplicated
import pandas as pd

gene_counts = pd.Series(adata.var_names).value_counts()
duplicates = gene_counts[gene_counts > 1]

print(f"🧠 Found {len(duplicates)} duplicated gene names.")
print("Top duplicates:")
print(duplicates.head(10))

# Just in case, apply one more time
adata.var_names_make_unique()

# === SAVE ===
adata.write(OUTPUT_FILE)
print(f"✅ Saved merged AnnData object to: {OUTPUT_FILE}")

