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