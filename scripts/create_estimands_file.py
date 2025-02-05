import yaml
import pandas as pd
import numpy as np

# TARGENE specifications
outcome_extra_covariates = ["Age-Assessment","Genetic-Sex"]
extra_treatments = []
estimands_orders = [2,3]
estimands_type = "AIE"
estimands_configuration_type = "groups"

# Final SNPs
bqtls = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_bQTLs.csv")
transactors = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_transactors.csv")

## Remove TFs in bqtls that are not contained in transactors
bqtls = bqtls[bqtls['TF'].isin(transactors.TF.unique())]

bqtl_dictionary = {key: bqtls.loc[bqtls['TF'] == key, 'SNP'].tolist() for key in bqtls['TF']}
transactors_dictionary = {key: transactors.loc[transactors['TF'] == key, 'SNP'].tolist() for key in transactors['TF']}

joined = {}

for key in bqtl_dictionary:
    joined[key] = {'bQTLs': bqtl_dictionary[key], 'transActors': transactors_dictionary[key]}

# Generate complete dictionary
if len(extra_treatments) == 0:
    data = {
        'type': estimands_configuration_type,
        'estimands': [
            {'type': estimands_type, 'orders': estimands_orders}
        ],
        'outcome_extra_covariates': outcome_extra_covariates,
        'variants': joined
    }
else:
    data = {
        'type': estimands_configuration_type,
        'estimands': [
            {'type': estimands_type, 'orders': estimands_orders}
        ],
        'outcome_extra_covariates': outcome_extra_covariates,
        'extra_treatments' : extra_treatments,
        'variants': joined
    } 

# Write to YAML
with open('estimands.yaml', 'w') as file:
    yaml.dump(data, file, default_flow_style=False, sort_keys=False)