process create_assoc_file {
    container 'roskamsh/bgen_env:0.2.0'
    publishDir("${launchDir}/output/assoc_files", pattern: "*.assoc")

    input:
        tuple val(chr), val(tf), path(snps), val(prefix), path(bed_files)
        path script

    output:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path("${tf}_ciseQTLs_hg38.assoc")

    script:
        """
        python ${script} --snps ${snps} --tf ${tf} --prefix ${prefix}_info_score_0.6_chr${chr}
        """
}

process ld_clump {
    container 'roskamsh/plink1.9:0.1.1'
    publishDir "${launchDir}/output/clumps"

    input:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path(assoc)

    output:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path("${tf}_clumped_r${params.R2_THRESHOLD}.clumped")

    script:
        """
        plink --bfile ${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr} --clump ${assoc} --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 ${params.R2_THRESHOLD} --out ${tf}_clumped_r${params.R2_THRESHOLD}
        """
}

process create_eqtl_list {
    container 'roskamsh/bgen_env:0.2.0'
    
    input:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path(clumps)

    output:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path("${tf}_independent_eQTLs.csv")

    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        clumps = pd.read_csv("${clumps}", delim_whitespace=True)
        final_eqtls = clumps[["SNP"]]
        tf = "${tf}"
        final_eqtls.to_csv(f"{tf}_independent_eQTLs.csv", index = False)
        """

}