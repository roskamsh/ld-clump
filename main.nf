#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow create_input_channels {
    main:
    Channel
        .fromPath(params.TFs, checkIfExists: true)
        .splitCsv(header:true)
        .multiMap {
            row -> 
            chr_tf_ch: tuple(row.CHR, row.TF)
            chr_ch: row.CHR
        }
        .set { tfs_ch }

    Channel
        .fromFilePairs(params.BGEN_FILES, size: -1, checkIfExists: true)
        .multiMap {
            prefix, files ->
            prefix: prefix
            files: [files]
        }
        .set { bgen_files_ch }

    Channel
        .fromPath(params.SNPs, checkIfExists: true)
        .set { snps_ch }

    // Join SNPs file with each TF/CHR
    tfs_ch.chr_tf_ch
        .combine(snps_ch)
        .set { tfs_snps_ch }
    
    // Collapse on duplicate chromosomes and join with BGEN files
    tfs_ch.chr_ch
        .unique()
        .combine(bgen_files_ch.prefix)
        .combine(bgen_files_ch.files)
        .set { chr_bgen_ch }

    emit:
        tfs_bgen = tfs_snps_ch
        chr_bgen = chr_bgen_ch
    
}

workflow reconfigure_channels {
    take:
        bed_ch
        tf_ch
        chr_ch

    main:
        bed_ch
            .map { chr, bed, bim, fam -> [bed, bim, fam] }
            .collect()
            .map { files -> [files] }
            .set { bed_files_ch }

        chr_ch
            .map { chr, prefix, files -> prefix }
            .unique()
            .set { prefix_ch }
        
        tf_ch
            .combine(prefix_ch)
            .combine(bed_files_ch)
            .set { tf_bed_ch } 

    emit:
        tf_bed_ch
}

process generate_info_score {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'

    input:
        tuple val(chr), val(prefix), path(files)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}.snpstats")

    script:
        """
        qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -snp-stats -osnp chr${chr}.snpstats
        """
}

process find_exclusion_snps {
    container 'roskamsh/bgen_env:0.2.0'
    input:
        tuple val(chr), val(prefix), path(files), path(snpstats)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}_rsids2exclude_info_score0.6.txt")

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=9)
        df = df[df["info"] < 0.6]

        snps = df[["alternate_ids","rsid","chromosome","position","alleleA","alleleB"]]
        snps.columns = ["SNPID","rsid","chromosome","position","alleleA","alleleB"] 
        rsids = " ".join(snps.rsid.values)

        # Write rsids
        with open("chr${chr}_rsids2exclude_info_score0.6.txt", "w") as text_file:
            text_file.write(rsids)
        """
}

process bgen_to_bed {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'

    input:
        tuple val(chr), val(prefix), path(files), path(exclusion_snps)

    output:
        tuple val(chr), path("${prefix}_info_score_0.6_chr${chr}.bed"), path("${prefix}_info_score_0.6_chr${chr}.bim"), path("${prefix}_info_score_0.6_chr${chr}.fam")

    script:
        """
        qctool -g ${prefix}${chr}.bgen \
               -s ${prefix}${chr}.sample \
               -excl-rsids ${exclusion_snps} \
               -og "${prefix}_info_score_0.6_chr${chr}.bed"
        """
}

process create_assoc_file {
    container 'roskamsh/bgen_env:0.2.0'

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
    container 'biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
    publishDir "${launchDir}/output/clumps"

    input:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path(assoc)

    output:
        path "${tf}_clumped_r0.2*"

    script:
        """
        plink --bfile ${prefix}_info_score_0.6_chr${chr} --clump ${assoc} --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 0.2 --out ${tf}_clumped_r0.2
        """
}

workflow {
    create_input_channels()

    generate_info_score(create_input_channels.out.chr_bgen)

    find_exclusion_snps(generate_info_score.out)

    bgen_to_bed(find_exclusion_snps.out)

    reconfigure_channels(bgen_to_bed.out, create_input_channels.out.tfs_bgen, create_input_channels.out.chr_bgen)
    
    create_assoc_file(reconfigure_channels.out, file("${projectDir}/py/create_assoc_file.py"))

    ld_clump(create_assoc_file.out)
}
