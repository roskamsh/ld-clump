include { merge_beds } from '../modules/interactions.nf'

workflow check_interactions {
    take:
        bed_files
    main:
        merge_beds(bed_files)

}