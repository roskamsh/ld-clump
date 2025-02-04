include { check_ld; compile_ld_results; create_bqtl_lists; merge_QTLs; generate_estimands } from '../modules/interactions.nf'

workflow check_interactions {
    take:
        bqtls_by_tf
        transactors_by_tf
        bed_files

    main:
        // Check whether bQTLs are in LD & save this information
        check_ld(bqtls_by_tf.join(transactors_by_tf).combine(bed_files))

        compile_ld_results(check_ld.out.collect())
        
        // by joining together, we ensure only TFs with transactors are passed on
        create_bqtl_lists(bqtls_by_tf.join(transactors_by_tf)) 

        create_bqtl_lists.out
            .map { tf, bqtls, transactors -> ["group", tf, bqtls, transactors]}
            .groupTuple()
            .set { final_qtls }
        
        merge_QTLs(final_qtls)

        generate_estimands(merge_QTLs.out.bQTLs, merge_QTLs.out.transactors)

}
