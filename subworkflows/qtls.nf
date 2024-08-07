include { pull_eQTLs; read_and_filter_bQTLs; identify_eqtl_chrs } from '../modules/qtls.nf'

workflow preprocess_qtl_data {
    main:
        // Create channel for output directory from baal-nf, and filter baal-nf results
        bqtl_dir_ch = Channel.fromPath(params.BAALNF_OUTDIR, type: 'dir')
        read_and_filter_bQTLs(bqtl_dir_ch)

        // Define TFs channel based on bQTLs which meet thresholds
        read_and_filter_bQTLs.out.tfs 
            .splitCsv()
            .flatten()
            .set { tfs_ch }

        // Pull eQTLGen data if not available locally and append to TFs channel
        if (params.eQTLGEN_DATA == "") {
            eqtls_ch = pull_eQTLs()
        } else {
            eqtls_ch = Channel.fromPath(params.eQTLGEN_DATA, checkIfExists: true)
        }
        tfs_eqtls_ch = tfs_ch.combine(eqtls_ch)

        // Define CHRs that need to be looked at based on CHRs which eQTLs and bQTLs belong to
        identify_eqtl_chrs(tfs_eqtls_ch, file("${projectDir}/bin/get_chr.py"))

        // Remove any TFs that don't exist in the eQTLGen database
        identify_eqtl_chrs.out
            .filter { it[1] != "" }
            .set { tf_chr_eqtls_ch }

        // Identify CHRs either eQTLs or bQTLs belong to
        read_and_filter_bQTLs.out.tf_chr_bqtls
            .splitCsv(header: true)
            .multiMap { 
                row -> 
                tf_chr_bqtl: [row.tf, row.CHROM, row.ID]
                chr: [row.CHROM]
            }
            .set { bqtls_ch }

        tf_chr_eqtls_ch
            .map { tf, chr, data -> chr }
            .mix(bqtls_ch.chr.flatten())
            .unique()
            .set { chr_ch }

        bqtls_ch.tf_chr_bqtl
            .groupTuple(by: [0,1])
            .set { tf_chr_bqtls_ch }

    emit:
    chrs = chr_ch
    tf_chr_eqtls = tf_chr_eqtls_ch 
    tf_chr_bqtls = tf_chr_bqtls_ch
}