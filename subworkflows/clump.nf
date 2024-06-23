include { create_assoc_file; ld_clump } from '../modules/clump.nf'

workflow generate_independent_snps {
    take:
        tf_bed_ch

    main:
        create_assoc_file(tf_bed_ch, file("${projectDir}/py/create_assoc_file.py"))

        ld_clump(create_assoc_file.out)
}