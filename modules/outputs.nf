process create_bqtl_lists {
    label 'bgen_python_image'

    input:
        tuple val(tf), val(bqtls), path(transactors)

    output:
        tuple val(tf), path("${tf}_bQTLs.csv"), path(transactors)

    script:
        """
        #!/usr/bin/env python

        import pandas as pd

        snps = "${bqtls}".strip("[]").replace(" ", "")
        snp_list = snps.split(',')
        df = pd.DataFrame(snp_list, columns=['SNP'])
        df["TF"] = "${tf}"
        df.to_csv("${tf}_bQTLs.csv", index = False)
        """
}

process merge_QTLs {
    label 'bgen_python_image'

    input:
        tuple val(group), val(tfs), path(transactors), path(bqtls)

    output:
        path "transactors.csv", emit: transactors
        path "bQTLs.csv", emit: bQTLs

    script:
        """
        #!/usr/bin/env python 

        import os
        import pandas as pd

        # List to hold the DataFrames
        transactor_list = []
        bqtl_list = []

        # Loop through files in the directory
        for filename in os.listdir("."):
            if "transactors" in filename: 
                df = pd.read_csv(filename)  
                transactor_list.append(df) 
            if "bQTLs" in filename:
                df = pd.read_csv(filename)  
                bqtl_list.append(df)  

        # Concatenate all DataFrames in the list
        transactor_df = pd.concat(transactor_list, ignore_index=True)
        bqtl_df = pd.concat(bqtl_list, ignore_index=True)

        # Save the combined DataFrame to a CSV file (optional)
        transactor_df.to_csv('transactors.csv', index=False)
        bqtl_df.to_csv('bQTLs.csv', index=False) 
        """
}


process generate_estimands {
    label 'bgen_python_image'
    publishDir "$params.OUTDIR/estimands", mode: 'copy'

    input:
        path bQTLs
        path transactors

    output:
        path "estimands.yaml"

    script:

        """
        #!/usr/bin/env python

        import yaml
        import pandas as pd
        import numpy as np
        import json

        # Custom Dumper class to enforce double quotes on strings in output YAML
        class QuotedString(str):
            pass

        def quoted_presenter(dumper, data):
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')

        yaml.add_representer(QuotedString, quoted_presenter)

        # Function to convert all strings in a dictionary to a QuotedString
        def quote_strings(obj):
            if isinstance(obj, str):
                return QuotedString(obj)
            elif isinstance(obj, list):
                return [quote_strings(item) for item in obj]
            elif isinstance(obj, dict):
                return {key: quote_strings(value) for key, value in obj.items()}
            return obj

        # Custom Dumper to ensure correct indentation
        class IndentedDumper(yaml.Dumper):
            def increase_indent(self, flow=False, indentless=False):
                return super(IndentedDumper, self).increase_indent(flow, False)

        # TARGENE specifications
        outcome_extra_covariates = ${params.OUTCOME_EXTRA_COVARIATES.collect{"'$it'"}}
        extra_treatments = ${params.EXTRA_TREATMENTS.collect{"'$it'"}}
        estimands_orders = ${params.ESTIMANDS_ORDERS}
        estimands_type = "${params.ESTIMANDS_TYPE}"
        estimands_configuration_type = "${params.ESTIMANDS_CONFIGURATION_TYPE}"

        # Final SNPs
        bqtls = pd.read_csv("${bQTLs}")
        transactors = pd.read_csv("${transactors}")

        ## Remove TFs in bqtls that are not contained in transactors
        bqtls = bqtls[bqtls['TF'].isin(transactors.TF.unique())]

        ## Create dictionaries for bQTLs/transActors
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

        if len(extra_treatments) > 0:
            data['extra_treatments'] = extra_treatments

        # Ensure quoted strings
        quoted_data = quote_strings(data)

        # Write to YAML with correct indentation
        with open('estimands.yaml', 'w') as file:
            yaml.dump(quoted_data, file, default_flow_style=False, sort_keys=False, width=100, indent=2, Dumper=IndentedDumper)
        """
}

process create_ld_block_input_file {
    label 'bgen_python_image'
    publishDir "$params.OUTDIR/snps", mode: 'copy'

    input:
        path final_bqtls
        path final_transactors
        path input_transactors
        path input_bqtls

    output:
        path "snps.csv"
        path "final_bQTLs.csv"
        path "final_transactors.csv"

    script:
        """
        #!/usr/bin/env python 
        import pandas as pd
        import os

        # Read in CSV files
        final_bqtls = pd.read_csv("${final_bqtls}")
        final_transactors = pd.read_csv("${final_transactors}")
        input_bqtls = pd.read_csv("${input_bqtls}")

        # Transactors are split by TF so read in each file and concatenate together
        transactor_list = []
        for filename in os.listdir("."):
            if "transactor_QTLs.csv" in filename: 
                df = pd.read_csv(filename)  
                transactor_list.append(df) 
        input_transactors = pd.concat(transactor_list, ignore_index=True)

        annot_bqtls = final_bqtls.merge(input_bqtls, 
                                        left_on = ['SNP','TF'], 
                                        right_on = ['ID','tf'], 
                                        how = 'inner')
        
        # Write annotated final bQTLs to a file
        annot_bqtls.drop(columns=["ID","tf"])
        annot_bqtls.to_csv("final_bQTLs.csv", index = False)

        # Remove chr prefix from chromosome column
        annot_bqtls["CHR_num"] = annot_bqtls['CHROM'].str.replace(r'^CHR', '', case=False, regex=True)
        annot_bqtls = annot_bqtls[["SNP","CHR_num","POS"]]
        annot_bqtls.columns = ["RSID","CHR","POS"]

        annot_transactors = final_transactors.merge(input_transactors,
                                                    left_on = ['SNP','TF'],
                                                    right_on = ['SNP','GeneSymbol'],
                                                    how = 'inner')

        ## Write annotated final transactors to a file
        annot_transactors.drop(columns=["TF"])
        annot_transactors.to_csv("final_transactors.csv", index = False)

        annot_transactors = annot_transactors[["SNP","SNPChr","SNPPos"]]
        annot_transactors.columns = ["RSID","CHR","POS"] 

        final = pd.concat([annot_bqtls, annot_transactors], ignore_index=True)
	final = final.drop_duplicates(subset = ["RSID"])
        final.to_csv("snps.csv", index = False)
        """
}
