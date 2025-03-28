# Load necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Load the .RData file
load("/Users/aumchampaneri/Databases/facs_Kidney_seurat_tiss.Robj")

# Check the structure of the loaded data
print(ls())