include { pull_eQTLs; read_and_filter_bQTLs; identify_transactors } from '../modules/qtls.nf'

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
            eqtls_file_ch = pull_eQTLs()
        } else {
            eqtls_file_ch = Channel.fromPath(params.eQTLGEN_DATA, checkIfExists: true)
        }

        // Identify transactors 
        transactors_file_ch = Channel.fromPath(params.ADDITIONAL_TRANSACTORS)
        identify_transactors(tfs_ch.combine(eqtls_file_ch).combine(transactors_file_ch), file("${projectDir}/bin/get_transactors.py"))

        // Create channel for all transactors with CHR & SNP information
        identify_transactors.out.map { it -> it[1] }
            .splitCsv(header: true)
            .map { row -> [row.SNPChr, row.SNP] }
            .set { transactor_chr_snps_ch }

        // Create channel for all bqtls with CHR & SNP information
        read_and_filter_bQTLs.out.tf_chr_bqtls
            .splitCsv(header: true)
            .map { row -> [row.CHROM, row.ID] }
            .set { bqtl_chr_snps_ch } 

        // Define all SNPs investigated per-chromosome for filtering & merging of BED files
        snps_by_chr = transactor_chr_snps_ch.mix(bqtl_chr_snps_ch).groupTuple()
        
        // Define all bQTLs per-TF for later LD checks
        read_and_filter_bQTLs.out.tf_chr_bqtls
            .splitCsv(header: true)
            .map { row -> [row.tf, row.ID] } 
            .groupTuple()
            .set { bqtls_by_tf }
        
    emit:
    bqtls_csv_file = read_and_filter_bQTLs.out.bqtls
    bqtls_by_tf = bqtls_by_tf
    transactors_csv_by_tf = identify_transactors.out
    snps_by_chr = snps_by_chr
    
}