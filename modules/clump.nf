process create_assoc_file {
    label 'bgen_python_image'
    publishDir("$params.OUTDIR/assoc_files", pattern: "*.assoc")

    input:
        tuple val(tf), path(snps), path(bed), path(bim), path(fam)
        path script

    output:
        tuple val(tf), path(bed), path(bim), path(fam), path("${tf}_transactors_hg38.assoc"), optional: true

    script:
        """
        python ${script} --snps ${snps} --tf ${tf} --prefix merged
        """
}

process ld_clump {
    label 'plink_image'
    label 'moremem'
    publishDir "$params.OUTDIR/clumps"

    input:
        tuple val(tf), path(bed), path(bim), path(fam), path(assoc)

    output:
        tuple val(tf), path("${tf}_clumped_r${params.R2_THRESHOLD}.clumped")

    script:
        """
        plink --bfile merged --clump ${assoc} --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 ${params.R2_THRESHOLD} --out ${tf}_clumped_r${params.R2_THRESHOLD}
        """
}

process create_transactor_list {
    label 'bgen_python_image'
    
    input:
        tuple val(tf), path(clumps)

    output:
        tuple val(tf), path("${tf}_independent_transactors.csv")

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        clumps = pd.read_csv("${clumps}", delim_whitespace=True)
        final_transactors = clumps[["SNP"]]
        tf = "${tf}"
        final_transactors[["TF"]] = tf
        final_transactors.to_csv(f"{tf}_independent_transactors.csv", index = False)
        """

}