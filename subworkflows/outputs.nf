include { create_bqtl_lists; merge_QTLs; generate_estimands; create_ld_block_input_file } from '../modules/outputs.nf'

workflow create_outputs {
    take:
        bqtls_by_tf
        transactors_by_tf
        bqtls_csv
        transactors_csv_by_tf

    main:
        // by joining together, we ensure only TFs with transactors are passed on
        create_bqtl_lists(bqtls_by_tf.join(transactors_by_tf)) 

        create_bqtl_lists.out
            .map { tf, bqtls, transactors -> ["group", tf, bqtls, transactors]}
            .groupTuple()
            .set { final_qtls }
        
        merge_QTLs(final_qtls)

        generate_estimands(merge_QTLs.out.bQTLs, merge_QTLs.out.transactors)

        transactors_csv_by_tf
            .map {_tf, csv -> csv }
            .collect()
            .set{ transactors_csv_files }

        create_ld_block_input_file(merge_QTLs.out.bQTLs, merge_QTLs.out.transactors, transactors_csv_files, bqtls_csv)
}