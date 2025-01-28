include { generate_info_score; find_inclusion_snps; bgen_to_bed; merge_beds } from '../modules/bgen.nf'

workflow preprocess_genetic_data {
    take:
        snps_by_chr_ch

    main:
        // Create BGEN channel for relevant chromosomes
        Channel
            .fromFilePairs(params.BGEN_FILES, size: -1, checkIfExists: true)
            .multiMap {
                prefix, files ->
                prefix: prefix
                files: [files]
            }
            .set { bgen_files_ch }

        // Generate info score for SNPs-of-interest
        generate_info_score(snps_by_chr_ch.combine(bgen_files_ch.prefix).combine(bgen_files_ch.files))

        // Find any SNPs below INFO score threshold 
        find_inclusion_snps(generate_info_score.out.combine(bgen_files_ch.prefix).combine(bgen_files_ch.files))

        // Convert to BED file for LD clumping
        bgen_to_bed(find_inclusion_snps.out)

        // Merge BED files
        merge_beds(bgen_to_bed.out.collect())

    emit:
        merge_beds.out
}