#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// LD thresholds
params.R2_THRESHOLD = 0.8
params.INFO_THRESHOLD = 0.9
params.MAF_THRESHOLD = 0.01
params.ASSEMBLY = "hg19"

// QTL input data and thresholds
params.eQTLGEN_DATA = ""
params.ADDITIONAL_TRANSACTORS = "${projectDir}/assets/NO_ADDITIONAL_TRANSACTORS"
params.ASB_quality = "High"

// TARGENE specifications
params.OUTCOME_EXTRA_COVARIATES = ["Age-Assessment","Genetic-Sex"]
params.ESTIMANDS_ORDERS = [2]
params.ESTIMANDS_TYPE = "AIE"
params.EXTRA_TREATMENTS = []
params.ESTIMANDS_CONFIGURATION_TYPE = "groups"

// Output specifications 
params.OUTDIR = "${launchDir}/output"

include { preprocess_qtl_data } from './subworkflows/qtls.nf'
include { preprocess_genetic_data } from './subworkflows/bgen.nf'
include { generate_independent_snps } from './subworkflows/clump.nf'
include { check_interactions } from './subworkflows/interactions.nf'
include { create_outputs } from './subworkflows/outputs.nf'

workflow {
    preprocess_qtl_data()

    preprocess_genetic_data(preprocess_qtl_data.out.snps_by_chr)

    generate_independent_snps(preprocess_qtl_data.out.transactors_csv_by_tf, preprocess_genetic_data.out.bed_ch)

    check_interactions(preprocess_qtl_data.out.bqtls_by_tf,  generate_independent_snps.out, preprocess_genetic_data.out.bed_ch)

    create_outputs(check_interactions.out.bqtls_by_tf, check_interactions.out.transactors_by_tf, 
                    preprocess_qtl_data.out.bqtls_csv_file, preprocess_qtl_data.out.transactors_csv_by_tf)
}
