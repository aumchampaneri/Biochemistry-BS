import cellxgene_census

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

Human Kidney Data
> https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc

Mouse Kidney Data
> https://cellxgene.cziscience.com/collections/92fde064-2fb4-41f8-b85c-c6904000b859

'''

# cellxgene_census.download_source_h5ad(
#     "f95d8919-1f2a-405f-8776-bfecc0ab0f3f",
#     to_path="/Volumes/CHAMPANERI/Databases/KidneyData.h5ad",
# )

cellxgene_census.download_source_h5ad(
    "42bb7f78-cef8-4b0d-9bba-50037d64d8c1",
    to_path="/Volumes/CHAMPANERI/Databases/MouseKidneyData.h5ad",
)