include { create_assoc_file; ld_clump; create_eqtl_list } from '../modules/clump.nf'

workflow generate_independent_snps {
    take:
        tf_bed_ch

    main:
        create_assoc_file(tf_bed_ch, file("${projectDir}/bin/create_assoc_file.py"))

        ld_clump(create_assoc_file.out)

        create_eqtl_list(ld_clump.out)

    emit:
        create_eqtl_list.out
}