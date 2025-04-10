import pandas as pd
import scanpy as sc

# ==== SET YOUR PATHS HERE ====
input_file = '/Users/aumchampaneri/Databases/Mm DKD Dataset/GSE181382_Expr.txt'           # Replace with your file path
output_file = '/Users/aumchampaneri/Databases/Mm DKD Dataset/GSE181382_Expr.h5ad'    # Desired name for output

# ==== LOAD EXPRESSION DATA ====
print(f"Reading expression matrix from {input_file}...")
expr = pd.read_csv(input_file, sep="\t", index_col=0)  # rows = genes, columns = cells

# ==== CONVERT TO AnnData ====
adata = sc.AnnData(expr.T)  # Transpose so cells = rows, genes = columns
adata.var_names_make_unique()

# ==== SAVE AS .h5ad ====
print(f"Saving AnnData object to {output_file}...")
adata.write(output_file)

print("✅ Done! Saved:", output_file)

print("Shape (cells x genes):", adata.shape)
print("First few cell names:", adata.obs_names[:5].tolist())
print("First few gene names:", adata.var_names[:5].tolist())