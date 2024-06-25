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
    
    // Collapse on duplicate chromosomes and combine with BGEN files
    tfs_ch.chr_ch
        .unique()
        .combine(bgen_files_ch.prefix)
        .combine(bgen_files_ch.files)
        .set { chr_bgen_ch }

    emit:
        tfs = tfs_ch.chr_tf_ch
        chr_bgen = chr_bgen_ch
}

workflow create_tf_bed_channel {
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


// @Olivier - best way to group bQTL-eQTL for looking at pairwise LD score?
// Currently, thinking about retaining CHR/TF information for each, having the bQTL be the index SNP
// Then have a list of eQTLs for the same TF, and if CHR for bQTL does not match CHR for eQTLs, break
// If CHR for bQTL does match CHR for eQTLs, then run LD pairwise calculation using plink2 (code in scripts/7_check_LD_pairs.sh)
/*
workflow create_bqtl_eqtl_channel {
    take:
        bqtl_ch
        eqtl_ch
    
    main:
        bqtl_ch.
    emit:

}
*/