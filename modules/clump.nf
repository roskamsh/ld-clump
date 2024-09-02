process create_assoc_file {
    label 'bgen_python_image'
    publishDir("$params.OUTDIR/assoc_files", pattern: "*.assoc")

    input:
        tuple val(tf), val(chr), path(snps), val(prefix), path(bed_files)
        path script

    output:
        tuple val(tf), val(chr), val(prefix), path(bed_files), path("${tf}_ciseQTLs_hg38.assoc"), optional: true

    script:
        """
        python ${script} --snps ${snps} --tf ${tf} --prefix ${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr}
        """
}

process ld_clump {
    label 'plink_image'
    publishDir "$params.OUTDIR/clumps"

    input:
        tuple val(tf), val(chr), val(prefix), path(bed_files), path(assoc)

    output:
        tuple val(tf), val(chr), val(prefix), path(bed_files), path("${tf}_clumped_r${params.R2_THRESHOLD}.clumped")

    script:
        """
        plink --bfile ${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr} --clump ${assoc} --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 ${params.R2_THRESHOLD} --out ${tf}_clumped_r${params.R2_THRESHOLD}
        """
}

process create_eqtl_list {
    label 'bgen_python_image'
    
    input:
        tuple val(tf), val(chr), val(prefix), path(bed_files), path(clumps)

    output:
        tuple val(tf), val(chr), val(prefix), path(bed_files), path("${tf}_independent_eQTLs.csv")

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        clumps = pd.read_csv("${clumps}", delim_whitespace=True)
        final_eqtls = clumps[["SNP"]]
        tf = "${tf}"
        final_eqtls[["TF"]] = tf
        final_eqtls.to_csv(f"{tf}_independent_eQTLs.csv", index = False)
        """

}