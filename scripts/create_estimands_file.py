import yaml
import pandas as pd
import numpy as np

# Custom Dumper class to enforce double quotes on strings in output YAML
class QuotedString(str):
    pass

def quoted_presenter(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')

yaml.add_representer(QuotedString, quoted_presenter)

# Function to convert all strings in a dictionary to a QuotedString
def quote_strings(obj):
    """Recursively convert strings in a dictionary or list to QuotedString."""
    if isinstance(obj, str):
        return QuotedString(obj)
    elif isinstance(obj, list):
        return [quote_strings(item) for item in obj]
    elif isinstance(obj, dict):
        return {key: quote_strings(value) for key, value in obj.items()}
    return obj

# Custom Dumper to ensure correct indentation
class IndentedDumper(yaml.Dumper):
    """Add indentation where necessary in output YAML."""
    def increase_indent(self, flow=False, indentless=False):
        return super(IndentedDumper, self).increase_indent(flow, False)

# TARGENE specifications
outcome_extra_covariates = ["Age-Assessment","Genetic-Sex"]
extra_treatments = []
estimands_orders = [2,3]
estimands_type = ["AIE"]
estimands_configuration_type = "groups"

# Final SNPs
bqtls = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_bQTLs.csv")
transactors = pd.read_csv("/gpfs/igmmfs01/eddie/ponting-lab/breeshey/projects/TARGENE/nhr/nr3c1/work/83/436338698806fc4886943abbf9d061/final_transactors.csv")

## Remove TFs in bqtls that are not contained in transactors
bqtls = bqtls[bqtls['TF'].isin(transactors.TF.unique())]

bqtl_dictionary = {key: bqtls.loc[bqtls['TF'] == key, 'SNP'].tolist() for key in bqtls['TF']}
transactors_dictionary = {key: transactors.loc[transactors['TF'] == key, 'SNP'].tolist() for key in transactors['TF']}

joined = {key: {'bQTLs': bqtl_dictionary.get(key, []), 'transActors': transactors_dictionary.get(key, [])} for key in bqtl_dictionary}

# Construct the YAML data
estimands_entry = {'orders': estimands_orders, 'type': estimands_type}

data = {
    'type': estimands_configuration_type,
    'estimands': [estimands_entry],
    'outcome_extra_covariates': outcome_extra_covariates,
    'variants': joined
}

if extra_treatments:
    data['extra_treatments'] = extra_treatments

# Ensure quoted strings
quoted_data = quote_strings(data)

# Write to YAML with correct indentation
with open('estimands.yaml', 'w') as file:
    yaml.dump(quoted_data, file, default_flow_style=False, sort_keys=False, width=100, indent=2, Dumper=IndentedDumper)
