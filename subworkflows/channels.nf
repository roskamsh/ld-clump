workflow create_tf_bed_channel {
    take:
        bed_ch
        tf_ch

    main:
        bed_ch
            .map { chr, prefix, bed, bim, fam -> [bed, bim, fam] }
            .collect()
            .map { files -> [files] }
            .set { bed_files_ch }

        bed_ch
            .map { chr, prefix, bed, bim, fam -> prefix }
            .unique()
            .set { prefix_ch }
        
        tf_ch
            .combine(prefix_ch)
            .combine(bed_files_ch)
            .set { tf_bed_ch } 

        prefix_ch 
            .combine(bed_files_ch)
            .set { bed_files_ch }

    emit:
        tf_bed = tf_bed_ch
        bed_files = bed_files_ch
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