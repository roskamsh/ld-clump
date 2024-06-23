include { pull_eQTLs; read_and_filter_bQTLs } from '../modules/qtls.nf'

workflow preprocess_qtl_data {
    take:
        tfs
    main:
        // Pull eQTLGen data if not available locally and append to TFs channel
        if (params.eQTLGEN_DATA == "") {
            eqtls_ch = pull_eQTLs()
        } else {
            eqtls_ch = Channel.fromPath(params.eQTLGEN_DATA, checkIfExists: true)
        }
        tfs_eqtls_ch = tfs.combine(eqtls_ch)

        // Create channel for output directory from baal-nf, append to TFs, then filter baal-nf results
        bqtl_dir_ch = Channel.fromPath(params.BAALNF_OUTDIR, type: 'dir')
        tfs_bqtls_ch = tfs.combine(bqtl_dir_ch)
        read_and_filter_bQTLs(tfs_bqtls_ch)

    emit:
    tfs_eqtls = tfs_eqtls_ch
    tfs_bqtls = tfs_bqtls_ch
}