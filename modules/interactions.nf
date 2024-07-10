process merge_beds {
    label 'bigmem'
    container "olivierlabayle/tl-core:0.8"
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        tuple val(prefix), path(files)
    
    output:
        path "merged*"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl \
        --input ${prefix}_info_score_${params.INFO_THRESHOLD}_chr \
        --output merged merge
        """


}