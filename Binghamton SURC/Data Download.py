import cellxgene_census

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

Human Kidney Data
- Harmonized data from 4 datasets
- Normal Single Cell = ~200k cells

> dataset_h5ad_path: f95d8919-1f2a-405f-8776-bfecc0ab0f3f
> https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc

Chronic Kidney Disease Data
- Includes Normal, Chronic Kidney Disease, and Acute Kidney Failure
- Normal Single Cell = ~44k cells
- CKD Single Cell = ~130k cells
- AKF Single Cell = ~51k cells
- Normal Single Nucleus = ~120k cells
- CKD Single Nucleus = ~130k cells
- AKF Single Nucleus = ~55k cells

> dataset_h5ad_path (scRNA): dea717d4-7bc0-4e46-950f-fd7e1cc8df7d
> dataset_h5ad_path (snRNA): a12ccb9b-4fbe-457d-8590-ac78053259ef
> https://cellxgene.cziscience.com/collections/0f528c8a-a25c-4840-8fa3-d156fa11086f

'''

# cellxgene_census.download_source_h5ad(
#     "f95d8919-1f2a-405f-8776-bfecc0ab0f3f",
#     to_path="/Volumes/CHAMPANERI/Databases/Human_Nor_Kidney_CellHint_scRNA.h5ad"
# )

# cellxgene_census.download_source_h5ad(
#     "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d",
#     to_path="/Volumes/CHAMPANERI/Databases/Human_Nor-CKD-AKF_scRNA.h5ad"
# )

cellxgene_census.download_source_h5ad(
    "a12ccb9b-4fbe-457d-8590-ac78053259ef",
    to_path="/Volumes/CHAMPANERI/Databases/Human_Nor-CKD-AKF_snRNA.h5ad"
)