include { generate_info_score; find_exclusion_snps; bgen_to_bed } from '../modules/snps.nf'

workflow preprocess_snps {
    take:
        chr_bgen_ch

    main:
        generate_info_score(chr_bgen_ch)

        find_exclusion_snps(generate_info_score.out)

        bgen_to_bed(find_exclusion_snps.out)

    emit:
        bgen_to_bed.out
}