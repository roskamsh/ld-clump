include { generate_info_score; find_exclusion_snps; bgen_to_bed } from '../modules/bgen.nf'

workflow preprocess_genetic_data {
    take:
        chr_ch

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
        
        // Combine with BGEN files
        chr_ch
            .combine(bgen_files_ch.prefix)
            .combine(bgen_files_ch.files)
            .set { chr_bgen_ch }
            
        generate_info_score(chr_bgen_ch)

        find_exclusion_snps(generate_info_score.out)

        bgen_to_bed(find_exclusion_snps.out)

    emit:
        bgen_to_bed.out
}