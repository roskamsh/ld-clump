include { filter_beds; merge_beds } from '../modules/interactions.nf'

workflow check_interactions {
    take:
        tf_chr_bqtls
        bed_files
    main:
        filter_beds(tf_chr_bqtls.combine(bed_files))

        filter_beds.out
            .map { tf, chr, prefix, files -> [files]}
            .flatten()
            .collect()
            .set { filtered_bed_files }
        
        merge_beds(filtered_bed_files)

}