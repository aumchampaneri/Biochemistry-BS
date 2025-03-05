import cellxgene_census

'''
Use this script to download the data from the cellxgene server
The data will be downloaded to the data folder in the root directory of the project
The data will be downloaded as a .h5ad file



'''

cellxgene_census.download_census_data(
    "xxx",
    to_path=""
)