process generate_info_score {
    label 'moremem'
    label 'moretime'
    label 'qctool_image'
    storeDir("$params.SNPSTATS_CACHE")

    input:
        tuple val(chr), val(snps), val(prefix), path(files)

    output:
        tuple val(chr), path("chr${chr}.snpstats")

    script:
        """
        echo "${snps}" | sed 's/[][]//g' | tr ', ' '\n' > snps.txt
        if [ -e "${prefix}${chr}.bgen" ]; then
            qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -incl-rsids snps.txt -snp-stats -osnp chr${chr}.snpstats
        else
            echo "No BGEN file provided for chromosome ${chr}, Creating empty file."
            touch chr${chr}.snpstats
        fi
        """
}

process find_inclusion_snps {
    label 'bgen_python_image'

    input:
        tuple val(chr), path(snpstats), val(prefix), path(files)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}_rsids2include_info_score${params.INFO_THRESHOLD}_maf${params.MAF_THRESHOLD}.txt"), optional: true

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        # Check if file is empty (ie. BGEN file not provided)
        with open("${snpstats}", 'r') as file:
            if file.read() == '':
                print("${snpstats} is empty, exiting process...")
            else:
                # Generate snps2exclude file
                info_score_threshold = float("${params.INFO_THRESHOLD}")
                maf_threshold = float("${params.MAF_THRESHOLD}")
                df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=10)

                # SNPs that fail QC measurements 
                rsidsPassQC = list(df[(df["info"] >= info_score_threshold) | (df["minor_allele_frequency"] >= maf_threshold)].rsid.values) 

                # Remove any multiallelic SNPs
                biallelic_df = df.drop_duplicates(subset = "rsid", keep = False)

                # Only keep SNPs passing QC thresholds
                rsidsPassQC = list(biallelic_df[
                                                    (biallelic_df["info"] >= info_score_threshold) | 
                                                    (biallelic_df["minor_allele_frequency"] >= maf_threshold)
                                                ]
                                            .rsid.values
                                    )

                rsids = " ".join(rsidsPassQC)

                # Write rsids
                with open("chr${chr}_rsids2include_info_score${params.INFO_THRESHOLD}_maf${params.MAF_THRESHOLD}.txt", "w") as text_file:
                    text_file.write(rsids)
        """
}

process bgen_to_bed {
    publishDir("$params.OUTDIR/filtered_bed_files")
    label 'qctool_image'

    input:
        tuple val(chr), val(prefix), path(files), path(inclusion_snps)

    output:
        path "info${params.INFO_THRESHOLD}.maf${params.MAF_THRESHOLD}.${prefix}${chr}*"

    script:
        """
        if [[ "${params.ASSEMBLY}" = "grch37" || "${params.ASSEMBLY}" = "hg19" ]]; then
            ranges="5:44000000-51500000 6:25000000-33500000 8:8000000-12000000 11:45000000-57000000"
        elif [[ "${params.ASSEMBLY}" = "grch38" || "${params.ASSEMBLY}" = "hg38" ]]; then
            ranges="5:43999898-52204166 6:24999772-33532223 8:8142478-12142491 11:44978449-57232526"
        else
            echo "Assembly provided does not match grch37/hg19 or grch38/hg38"
            exit 1
        fi

        qctool -g ${prefix}${chr}.bgen \
               -s ${prefix}${chr}.sample \
               -incl-rsids ${inclusion_snps} \
               -excl-range \$ranges \
               -og "info${params.INFO_THRESHOLD}.maf${params.MAF_THRESHOLD}.${prefix}${chr}.bed"
        """
}

process merge_beds {
    label 'moremem'
    label 'plink_image'
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        tuple path("merged.bed"), path ("merged.bim"), path("merged.fam")

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

        # Merge all BED files using plink
        plink --bfile "\$first_prefix" --merge-list "\$merge_list" --make-bed --out merged 
        """
}