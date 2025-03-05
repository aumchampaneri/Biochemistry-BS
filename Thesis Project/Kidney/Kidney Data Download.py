import cellxgene_census

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file

> https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc
'''

cellxgene_census.download_source_h5ad(
    "f95d8919-1f2a-405f-8776-bfecc0ab0f3f",
    to_path="KidneyData.h5ad",
)