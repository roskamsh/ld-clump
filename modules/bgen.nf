process generate_info_score {
    label 'moremem'
    label 'moretime'
    label 'qctool_image'
    storeDir("$params.SNPSTATS_CACHE")

    input:
        tuple val(chr), val(prefix), path(files)

    output:
        path "chr${chr}.snpstats"

    script:
        """
        if [ -e "${prefix}${chr}.bgen" ]; then
            qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -snp-stats -osnp chr${chr}.snpstats
        else
            echo "No BGEN file provided for chromosome ${chr}, Creating empty file."
            touch chr${chr}.snpstats
        fi
        """
}

process find_exclusion_snps {
    label 'bgen_python_image'

    input:
        tuple val(chr), path(snpstats), val(prefix), path(files)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}_rsids2exclude_info_score${params.INFO_THRESHOLD}_maf${params.MAF_THRESHOLD}.txt"), optional: true

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
                df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=9)

                # SNPs that fail QC measurements 
                rsidsfailQC = list(df[(df["info"] < info_score_threshold) | (df["minor_allele_frequency"] < maf_threshold)].rsid.values)

                # Identify SNPs to include after imposing these filter
                df2include = df[~df['rsid'].isin(rsidsfailQC)]

                # Identify any remaining multiallelc SNPs
                multiallelic_mask = df2include["rsid"].duplicated(keep=False)  # keep=False marks all duplicates as True
                # Filter the DataFrame to show only the duplicated rows
                multiallelic_df = df2include[multiallelic_mask]

                # Add multiallelic SNPs remaining to rsids2exclude
                multiallelic_snps = list(multiallelic_df.rsid.unique())
                rsids2exclude = rsidsfailQC + multiallelic_snps

                rsids = " ".join(rsids2exclude)

                # Write rsids
                with open("chr${chr}_rsids2exclude_info_score${params.INFO_THRESHOLD}_maf${params.MAF_THRESHOLD}.txt", "w") as text_file:
                    text_file.write(rsids) 
        """
}

process bgen_to_bed {
    publishDir("$params.OUTDIR/filtered_bed_files")
    label 'moremem'
    label 'qctool_image'

    input:
        tuple val(chr), val(prefix), path(files), path(exclusion_snps)

    output:
        tuple val(chr), val(prefix), path("filt.${prefix}${chr}.bed"), path("filt.${prefix}${chr}.bim"), path("filt.${prefix}${chr}.fam")

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
               -excl-rsids ${exclusion_snps} \
               -excl-range \$ranges \
               -og "filt.${prefix}${chr}.bed"
        """
}