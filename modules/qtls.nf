process pull_eQTLs {
    container 'roskamsh/commandlinetools:latest'

    when:
    params.eQTLGEN_DATA == ""

    output:
    path "*txt"

    script:
    """
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
    gunzip -d ./*gz
    """
}

process read_and_filter_bQTLs {
    container 'roskamsh/bgen_env:0.2.0'

    input:
    path DIR

    output:
    path "all_bQTLs.csv", emit: bqtls
    path "TF_list.csv", emit: tfs 
    path "TF_CHR_bQTL_list.csv", emit: tf_chr_bqtls

    script:
    """
    #!/usr/bin/env python
    import os
    import pandas as pd

    def remove_chr_prefix(s):
        if s.startswith('chr'):
            return s[3:]
        else:
            return s

    dir = "${DIR}"
    concordant_value = bool("${params.Concordant}")

    # Read in all CSV files, and filter by Concordant value
    files = os.listdir(dir)

    # Initialize an empty list to store DataFrames
    dataframes = []

    # Loop through the CSV files and read them into DataFrames
    for file in files:
        file_path = os.path.join(dir, file)
        df = pd.read_csv(file_path)
        if "${params.Concordant}" != "":
            df = df[df.Concordant == concordant_value] 
        dataframes.append(df)

    results = pd.concat(dataframes, ignore_index=True)

    # Define TF list from results
    tfs = pd.DataFrame(list(results.tf.unique()))
    # Define bQTL-CHR list
    filt = results[["tf","CHROM","ID"]].copy()
    filt["CHROM"] = filt["CHROM"].apply(remove_chr_prefix)
    filt = filt.drop_duplicates()

    results.to_csv("all_bQTLs.csv", index = False)
    tfs.to_csv("TF_list.csv", index = False, header=False)
    filt.to_csv("TF_CHR_bQTL_list.csv", index = False)
    """
}

process identify_eqtl_chrs {
    container 'roskamsh/bgen_env:0.2.0'

    input:
    tuple val(tf), path(eqtlgen_data)
    path script

    output:
    tuple val(tf), env(chr), path(eqtlgen_data)

    script:
    """
    python ${script} --data ${eqtlgen_data} --tf ${tf} > chr.txt
    if grep -q "Not found in eQTLGen" "chr.txt"; then
        echo "The TF, ${tf}, is not found in eQTLGen"
    else
        chr=\$(cat chr.txt)
    fi
    """ 

}
