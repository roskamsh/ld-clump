#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow create_input_channels {
    main:
    Channel
        .fromPath(params.TFs, checkIfExists: true)
        .splitCsv(header:true)
        .map {
            row -> tuple(
                row.TF,
                row.CHR
                )
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

    tfs_ch
        .combine(bgen_files_ch.prefix)
        .combine(bgen_files_ch.files)
        .set { tfs_bgen_ch }

    emit:
        tfs_bgen_ch
    
}

process generate_info_score {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'

    input:
        tuple val(tf), val(chr), val(prefix), path(files)

    output:
        tuple val(tf), val(chr), path("chr${chr}.snpstats"), emit: snpstats

    script:
    """
    qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -snp-stats -osnp chr${chr}.snpstats
    """
}

workflow {
    create_input_channels()

    generate_info_score(create_input_channels.out)
}
