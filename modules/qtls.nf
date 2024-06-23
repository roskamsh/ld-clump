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
    tuple val(CHR), val(TF), path(DIR)

    output:
    path "${TF}.filtered.csv"

    script:
    """
    #!/usr/bin/env python
    import os
    import pandas as pd

    dir = "${DIR}"
    tf = "${TF}"
    concordant_value = bool("${params.Concordant}")

    file_list = os.listdir(dir)
    # Find the filename in the cwd that contains our TF of interest
    tf_filename = next(
        (f for f in file_list if os.path.isfile(os.path.join(dir, f)) and tf in f),
        None
    )

    results = pd.read_csv(os.path.join(dir, tf_filename))

    # Filter by Concordant value if nextflow parameter is not empty
    if ${params.Concordant} != "":
        results = results[results.Concordant == concordant_value]

    results.to_csv(f"{tf}.filtered.csv", index = False)
    """
}

