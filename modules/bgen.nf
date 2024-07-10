process generate_info_score {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'
    publishDir("$params.OUTDIR/info_scores", pattern: "*.snpstats") 

    input:
        tuple val(chr), val(prefix), path(files)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}.snpstats")

    script:
        """
        qctool -g ${prefix}${chr}.bgen -s ${prefix}${chr}.sample -snp-stats -osnp chr${chr}.snpstats
        """
}

process find_exclusion_snps {
    container 'roskamsh/bgen_env:0.2.0'

    input:
        tuple val(chr), val(prefix), path(files), path(snpstats)

    output:
        tuple val(chr), val(prefix), path(files), path("chr${chr}_rsids2exclude_info_score${params.INFO_THRESHOLD}.txt")

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        info_score_threshold = float("${params.INFO_THRESHOLD}")
        df = pd.read_csv("${snpstats}", delimiter="\t", skiprows=9)
        df = df[df["info"] < info_score_threshold]

        snps = df[["alternate_ids","rsid","chromosome","position","alleleA","alleleB"]]
        snps.columns = ["SNPID","rsid","chromosome","position","alleleA","alleleB"] 
        rsids = " ".join(snps.rsid.values)

        # Write rsids
        with open("chr${chr}_rsids2exclude_info_score${params.INFO_THRESHOLD}.txt", "w") as text_file:
            text_file.write(rsids)
        """
}

process bgen_to_bed {
    label 'moremem'
    container 'roskamsh/qctools:0.1.1'

    input:
        tuple val(chr), val(prefix), path(files), path(exclusion_snps)

    output:
        tuple val(chr), val(prefix), path("${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr}.bed"), path("${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr}.bim"), path("${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr}.fam")

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
               -og "${prefix}_info_score_${params.INFO_THRESHOLD}_chr${chr}.bed"
        """
}