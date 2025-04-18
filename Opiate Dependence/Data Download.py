import cellxgene_census
from pandas.io.parquet import to_parquet

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

Run each of the download queries separately to download the data

Human Kidney Data
- Harmonized data from 4 datasets
- Normal Single Cell = ~200k cells

> dataset_h5ad_path: f95d8919-1f2a-405f-8776-bfecc0ab0f3f
> https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc

'''

# cellxgene_census.download_source_h5ad(
#     "f95d8919-1f2a-405f-8776-bfecc0ab0f3f",
#     # to_path="/Volumes/CHAMPANERI/Databases/Human_Nor_Kidney_CellHint_scRNA.h5ad" # USB Location
#     to_path="/Users/aumchampaneri/Databases/Human_Nor_Kidney_CellHint_scRNA.h5ad" # Disk Location
# )
