# Generate a CSV file with the datasets information from the census

import cellxgene_census

# Open the census
census = cellxgene_census.open_soma()

# Read and concatenate the datasets
census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

# Indexing by soma_joinid for convenience
census_datasets = census_datasets.set_index("soma_joinid")

# Save the DataFrame to a CSV file
census_datasets.to_csv('CELLxGENE Datasets.csv')

print("DataFrame saved to 'CELLxGENE Datasets.csv'")