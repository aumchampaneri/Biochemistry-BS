import cellxgene_census

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

Run each of the download queries separately to download the data

Mouse Kidney Atlas -> Mm_MKA.h5ad
- Normal single cell dataset -> harmonized snRNA-seq and scRNA-seq data @ 59 public datasets
- 141,401 cells in total

'''
# Define the file name for saving the dataset
file_name = "Mm_MKA.h5ad"  # Name of the file to save

# Download the Mouse Kidney Atlas dataset
cellxgene_census.download_source_h5ad(
    "42bb7f78-cef8-4b0d-9bba-50037d64d8c1",
    # to_path=f"/Volumes/CHAMPANERI/Databases/{file_name}" # USB Location
    to_path=f"/Users/aumchampaneri/Databases/{file_name}" # Disk Location
)