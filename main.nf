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
        if [[ "${params.ASSEMBLY}" = "grch37" || "${params.ASSEMBLY}" = "hg19" ]]; then
            ranges="5:44000000-51500000 6:25000000-33500000 8:8000000-12000000 11:45000000-57000000"
        elif [[ "${params.ASSEMBLY}" = "grch38" || "${params.ASSEMBLY}" = "hg38" ]]; then
            ranges="5:43999898-52204166 6:24999772-33532223 8:8142478-12142491 11:44978449-57232526"
        else
            echo "Assembly provided does not match grch37/hg19 or grch38/hg38"
            exit 1
        fi

        qctool -g ${prefix}${chr}.bgen \
               -s ${prefix}${chr}.sample \
               -excl-rsids ${exclusion_snps} \
               -excl-range \$ranges \
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
    container 'roskamsh/plink1.9:0.1.1'
    publishDir "${launchDir}/output/clumps"

    input:
        tuple val(chr), val(tf), val(prefix), path(bed_files), path(assoc)

    output:
        path "${tf}_clumped_r${params.R2_THRESHOLD}*"

    script:
        """
        plink --bfile ${prefix}_info_score_0.6_chr${chr} --clump ${assoc} --clump-p1 5e-6 --clump-p2 0.05 --clump-r2 ${params.R2_THRESHOLD} --out ${tf}_clumped_r${params.R2_THRESHOLD}
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
