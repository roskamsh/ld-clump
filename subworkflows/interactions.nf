include { check_ld; compile_ld_results } from '../modules/interactions.nf'

workflow check_interactions {
    take:
        bqtls_by_tf
        transactors_by_tf
        bed_files

    main:
        // Check whether bQTLs are in LD & save this information
        check_ld(bqtls_by_tf.join(transactors_by_tf).combine(bed_files))

        compile_ld_results(check_ld.out.collect())

    emit:
        bqtls_by_tf = bqtls_by_tf
        transactors_by_tf = transactors_by_tf

}
