include { create_assoc_file; ld_clump; create_transactor_list } from '../modules/clump.nf'

workflow generate_independent_snps {
    take:
        tf_transactors_ch
        bed_ch

    main:
        create_assoc_file(tf_transactors_ch.combine(bed_ch), file("${projectDir}/bin/create_assoc_file.py"))

        ld_clump(create_assoc_file.out)

        create_transactor_list(ld_clump.out)

    emit:
        create_transactor_list.out
}