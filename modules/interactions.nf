// This process filters each BED file by bQTLs of interest
// If a bQTL is not contained in a BED file, this means it did not meet INFO_THRESHOLD, and will be excluded
process filter_beds {
    label 'plink_image'
    label 'moremem'

    input:
        tuple val(tf), val(chr), val(snps), val(prefix), path(files)

    output:
        tuple val(tf), val(chr), val(snps), val(prefix), path("${tf}_chr${chr}_filtered.{bed,bim,fam}"), optional: true

    script:
        """
        echo "${snps}" | sed 's/[][]//g' | tr ', ' '\n' > snps.txt

        # Run this even if it encounters an error
        plink --bfile ${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr} --extract snps.txt --make-bed --out ${tf}_chr${chr}_filtered 2> plink_error.log || true

        # Check specific error here, if variants do not pass thresholds, still exit with exit code 0
        if grep -q "Error: No variants remaining after --extract" plink_error.log; then
            echo "PLINK error: No variants remaining after filtering. Removing these SNPs from further analysis."
            exit 0
        fi
        """

}


process merge_beds {
    label 'moremem'
    label 'plink_image'
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        tuple path("merged.bed"), path ("merged.bim"), path("merged.fam")

    script:
        """
        # Step 1: Identify all BED files and extract their prefixes
        bed_files=(\$(ls *.bed))
        prefixes=()

        for bed_file in "\${bed_files[@]}"; do
            prefix="\${bed_file%.bed}"
            prefixes+=("\$prefix")
        done

        # Check if there are enough files to merge
        if [ \${#prefixes[@]} -lt 2 ]; then
            echo "Not enough BED files to merge. At least two BED files are required."
            exit 1
        fi

        # Step 2: Create a merge list file (excluding the first prefix)
        merge_list="mergelist.txt"
        > \$merge_list

        for prefix in "\${prefixes[@]:1}"; do
            echo "\$prefix" >> \$merge_list
        done

        # Step 3: Run PLINK to merge the datasets
        first_prefix=\${prefixes[0]}
        output_prefix="merged_dataset"

        plink --bfile "\$first_prefix" --merge-list "\$merge_list" --make-bed --out merged

        """
}

process check_ld {
    label 'plink_image'

    input:
        tuple val(tf), val(snps), path(bed), path(bim), path(fam)

    output:
        path "ld_results.ld", optional: true

    script:
        """
        echo "${snps}" | sed 's/[][]//g' | tr ', ' '\n' > snps.txt
        
        plink --bfile merged --extract snps.txt --make-bed --out extracted_data
        plink --bfile extracted_data --r2 inter-chr --out ld_results
        """

}

process create_bqtl_lists {
    label 'bgen_python_image'

    input:
        tuple val(tf), val(snps)

    output:
        tuple val(tf), path("${tf}_bQTLs.csv")

    script:
        """
        #!/usr/bin/env python

        import pandas as pd

        snps = "${snps}".strip("[]").replace(" ", "")
        snp_list = snps.split(',')
        df = pd.DataFrame(snp_list, columns=['SNP'])
        df["TF"] = "${tf}"
        df.to_csv("${tf}_bQTLs.csv", index = False)
        """
}

process merge_QTLs {
    label 'bgen_python_image'
    publishDir "$params.OUTDIR/snps"

    input:
        tuple val(group), val(tfs), path(eqtls), path(bqtls)

    output:
        path "final_eQTLs.csv", emit: eQTLs
        path "final_bQTLs.csv", emit: bQTLs

    script:
        """
        #!/usr/bin/env python 

        import os
        import pandas as pd

        # List to hold the DataFrames
        eqtl_list = []
        bqtl_list = []

        # Loop through files in the directory
        for filename in os.listdir("."):
            if "eQTLs" in filename: 
                df = pd.read_csv(filename)  
                eqtl_list.append(df) 
            if "bQTLs" in filename:
                df = pd.read_csv(filename)  
                bqtl_list.append(df)  

        # Concatenate all DataFrames in the list
        eqtl_df = pd.concat(eqtl_list, ignore_index=True)
        bqtl_df = pd.concat(bqtl_list, ignore_index=True)

        # Save the combined DataFrame to a CSV file (optional)
        eqtl_df.to_csv('final_eQTLs.csv', index=False)
        bqtl_df.to_csv('final_bQTLs.csv', index=False) 
        """

}


process generate_estimands {
    label 'bgen_python_image'
    publishDir "$params.OUTDIR/estimands"

    input:
        path bQTLs
        path eQTLs

    output:
        path "estimands.yaml"

    script:
        """
        #!/usr/bin/env python

        import yaml
        import pandas as pd
        import numpy as np

        bqtls = pd.read_csv("${bQTLs}")
        eqtls = pd.read_csv("${eQTLs}")

        ## Remove TFs in bqtls that are not contained in eqtls
        bqtls = bqtls[bqtls['TF'].isin(eqtls.TF.unique())]

        bqtl_dictionary = {key: bqtls.loc[bqtls['TF'] == key, 'SNP'].tolist() for key in bqtls['TF']}
        eqtl_dictionary = {key: eqtls.loc[eqtls['TF'] == key, 'SNP'].tolist() for key in eqtls['TF']}

        joined = {}

        for key in bqtl_dictionary:
            joined[key] = {'bQTLs': bqtl_dictionary[key], 'eQTLs': eqtl_dictionary[key]}

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
