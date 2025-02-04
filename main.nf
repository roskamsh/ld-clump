#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.R2_THRESHOLD = 0.8
params.INFO_THRESHOLD = 0.9
params.MAF_THRESHOLD = 0.01
params.ASSEMBLY = "hg19"
params.eQTLGEN_DATA = ""
params.ADDITIONAL_TRANSACTORS = ""
params.ASB_quality = "High"
params.OUTDIR = "${launchDir}/output"
params.SNPSTATS_CACHE = "$params.OUTDIR/info_scores"

include { preprocess_qtl_data } from './subworkflows/qtls.nf'
include { preprocess_genetic_data } from './subworkflows/bgen.nf'
include { generate_independent_snps } from './subworkflows/clump.nf'
include { check_interactions } from './subworkflows/interactions.nf'

workflow {
    preprocess_qtl_data()

    preprocess_genetic_data(preprocess_qtl_data.out.snps_by_chr)

    generate_independent_snps(preprocess_qtl_data.out.transactors_by_tf, preprocess_genetic_data.out.bed_ch)

    check_interactions(preprocess_qtl_data.out.bqtls_by_tf,  generate_independent_snps.out, preprocess_genetic_data.out.bed_ch)
}
