import cellxgene_census
from pandas.io.parquet import to_parquet

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

Run each of the download queries separately to download the data

Human Dorsal Striatum in Opiate Use Disorder
- 10x 3' v3 = ~100k cells
- Caudate Nucleus and Putamen Tissue
- Normal and Opiate Dependence Disease States

> dataset_h5ad_path: c893ddc3-f25b-45e2-8c9e-155918b4261c.h5ad
> https://cellxgene.cziscience.com/collections/cec4ef8e-1e70-49a2-ae43-1e6bf1fd5978

'''

cellxgene_census.download_source_h5ad(
    "c893ddc3-f25b-45e2-8c9e-155918b4261c",
    # to_path="/Users/aumchampaneri/Databases/Opiate Dependence/GSE233279_CxG.h5ad" # USB Location
    to_path="/Users/aumchampaneri/Databases/Opiate Dependence/GSE233279_CxG.h5ad" # Disk Location
)
