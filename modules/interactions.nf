// This process filters each BED file by bQTLs of interest
// If a bQTL is not contained in a BED file, this means it did not meet INFO_THRESHOLD, and will be excluded
process filter_beds {
    container 'roskamsh/plink1.9:0.1.1'

    input:
        tuple val(tf), val(chr), val(snps), val(prefix), path(files)

    output:
        tuple val(tf), val(chr), val(prefix), path("${tf}_chr${chr}_filtered.{bed,bim,fam}"), optional: true

    script:
        """
        echo "${snps}" | sed 's/[][]//g' | tr ', ' '\n' > snps.txt

        # Run this even if it encounters an error
        plink --bfile ${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr} --extract snps.txt --make-bed --out ${tf}_chr${chr}_filtered 2> plink_error.log || true

        # Check specific error here, if variants do not pass thresholds, still exit with exit code 0
        if grep -q "Error: No variants remaining after --extract" plink_error.log; then
            echo "PLINK error: No variants remaining after filtering. Removing these SNPs from further analysis."
            exit 0
        fi
        """

}


process merge_beds {
    label 'bigmem'
    container 'roskamsh/plink1.9:0.1.1'
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        path "merged*"

    script:
        """
        # Step 1: Identify all BED files and extract their prefixes
        bed_files=(\$(ls *.bed))
        prefixes=()

        for bed_file in "\${bed_files[@]}"; do
            prefix="\${bed_file%.bed}"
            prefixes+=("\$prefix")
        done

        # Check if there are enough files to merge
        if [ \${#prefixes[@]} -lt 2 ]; then
            echo "Not enough BED files to merge. At least two BED files are required."
            exit 1
        fi

        # Step 2: Create a merge list file (excluding the first prefix)
        merge_list="mergelist.txt"
        > \$merge_list

        for prefix in "\${prefixes[@]:1}"; do
            echo "\$prefix" >> \$merge_list
        done

        # Step 3: Run PLINK to merge the datasets
        first_prefix=\${prefixes[0]}
        output_prefix="merged_dataset"

        plink --bfile "\$first_prefix" --merge-list "\$merge_list" --make-bed --out merged_dataset

        """
}
