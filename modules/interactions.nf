process check_ld {
    label 'plink_image'

    input:
        tuple val(tf), val(bqtls), path(transactors), path(bed), path(bim), path(fam)

    output:
        path "${tf}_ld_results.ld", optional: true

    script:
        """
        echo "${bqtls}" | sed 's/[][]//g' | tr ', ' '\n' > snps.txt
        
        plink --bfile merged --extract snps.txt \
              --r2 inter-chr --ld-window-r2 0 \
              --out ${tf}_ld_results
        """
}

process compile_ld_results {
    label 'bgen_python_image'
    publishDir"$params.OUTDIR/snps", mode: 'copy'

    input:
    path ld_files

    output:
    path "bQTLs_in_LD.csv", optional: true    

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import os

    r2_threshold = float(${params.R2_THRESHOLD})
    results_dir = "."
    files = [f for f in os.listdir(results_dir) if f.endswith('_ld_results.ld')]

    compiled_dataframes = []

    # Iterate through the list of files
    for file in files:
        df = pd.read_csv(os.path.join(results_dir,file), delim_whitespace=True, header=0)
        if df.empty:
            print(f"Skipping empty file: {file}. No variants in LD with one another.")
            continue
        # Add a new column for TF name
        df['TF'] = file.replace('_ld_results.ld', '')
        ld = df[df['R2'] >= r2_threshold].copy()
        compiled_dataframes.append(ld)

    if compiled_dataframes:
        compiled_df = pd.concat(compiled_dataframes, ignore_index=True)
        print("Number of bQTL pairs in LD across all TFs inspected:", compiled_df.shape[0])
        print("Saving bQTLs in LD to file: bQTLs_in_LD.csv")
        compiled_df.to_csv("bQTLs_in_LD.csv", index = False)
    else:
        compiled_df = pd.DataFrame()
        print("No bQTLs found to be in LD with one another. Exiting...")
    """

}

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
    publishDir "$params.OUTDIR/snps", mode: 'copy'

    input:
        tuple val(group), val(tfs), path(transactors), path(bqtls)

    output:
        path "final_transactors.csv", emit: transactors
        path "final_bQTLs.csv", emit: bQTLs

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
        transactor_df.to_csv('final_transactors.csv', index=False)
        bqtl_df.to_csv('final_bQTLs.csv', index=False) 
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

        bqtls = pd.read_csv("${bQTLs}")
        transactors = pd.read_csv("${transactors}")

        ## Remove TFs in bqtls that are not contained in eqtls
        bqtls = bqtls[bqtls['TF'].isin(transactors.TF.unique())]

        bqtl_dictionary = {key: bqtls.loc[bqtls['TF'] == key, 'SNP'].tolist() for key in bqtls['TF']}
        transactors_dictionary = {key: transactors.loc[transactors['TF'] == key, 'SNP'].tolist() for key in transactors['TF']}

        joined = {}

        for key in bqtl_dictionary:
            joined[key] = {'bQTLs': bqtl_dictionary[key], 'transactors': transactors_dictionary[key]}

        # Generate complete dictionary
        data = {
            'type': 'groups',
            'estimands': [
                {'type': 'IATE', 'orders': 2}
            ],
            'outcome_extra_covariates': ['Age','Sex'],
            'variants': joined
        }

        # Write to YAML
        with open('estimands.yaml', 'w') as file:
            yaml.dump(data, file, default_flow_style=False, sort_keys=False)

        """
}
