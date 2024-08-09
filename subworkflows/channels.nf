workflow create_tf_bed_channel {
    take:
        bed_ch
        tf_ch

    main:
        bed_ch
            .map { chr, prefix, bed, bim, fam -> [bed, bim, fam] }
            .collect()
            .map { files -> [files] }
            .set { bed_files_ch }

        bed_ch
            .map { chr, prefix, bed, bim, fam -> prefix }
            .unique()
            .set { prefix_ch }
        
        tf_ch
            .combine(prefix_ch)
            .combine(bed_files_ch)
            .set { tf_bed_ch } 

        prefix_ch 
            .combine(bed_files_ch)
            .set { bed_files_ch }

    emit:
        tf_bed = tf_bed_ch
        bed_files = bed_files_ch
}
