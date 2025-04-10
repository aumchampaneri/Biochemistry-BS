library(Seurat)
library(SeuratDisk)

# Load the .rds file
seurat_obj <- readRDS('/Users/aumchampaneri/Databases/Mm DKD Dataset/GSE181382.Expr.rds')

# Save as .h5Seurat
SaveH5Seurat(seurat_obj, filename = '/Users/aumchampaneri/Databases/Mm DKD Dataset/GSE181382.Expr.h5Seurat')

# Convert to .h5ad
Convert('/Users/aumchampaneri/Databases/Mm DKD Dataset/GSE181382.Expr.h5Seurat', dest = "h5ad")
