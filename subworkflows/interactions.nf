include { filter_beds; merge_beds; check_ld; compile_ld_results; create_bqtl_lists; merge_QTLs; generate_estimands } from '../modules/interactions.nf'

workflow check_interactions {
    take:
        tf_chr_bqtls
        bed_files
        tf_eqtls

    main:
        filter_beds(tf_chr_bqtls.combine(bed_files))

        filter_beds.out
            .map { tf, chr, snps, prefix, files -> [files]}
            .flatten()
            .collect()
            .set { filtered_bed_files }
        
        filter_beds.out 
            .map { tf, chr, snps, prefix, files -> [tf, snps] }
            .groupTuple()
            .map { tf, snps -> [tf, snps.flatten()]}
            .set { bqtls_per_tf }

        merge_beds(filtered_bed_files)

        check_ld(bqtls_per_tf.combine(merge_beds.out))

        compile_ld_results(check_ld.out.collect())
        
        create_bqtl_lists(bqtls_per_tf)

        tf_eqtls 
            .join(create_bqtl_lists.out) // Only TFs with eQTLs and bQTLs will remain
            .map { tf, eqtls, bqtls -> ["group", tf, eqtls, bqtls]}
            .groupTuple()
            .set { final_qtls }

        merge_QTLs(final_qtls)

        generate_estimands(merge_QTLs.out.bQTLs, merge_QTLs.out.eQTLs)

}
