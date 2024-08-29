#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.R2_THRESHOLD = 0.8
params.INFO_THRESHOLD = 0.9
params.ASSEMBLY = "hg19"
params.eQTLGEN_DATA = ""
params.Concordant = "True"
params.OUTDIR = "${launchDir}/output"
params.SNPSTATS_CACHE = "$params.OUTDIR/info_scores"

include { create_tf_bed_channel } from './subworkflows/channels.nf' 
include { preprocess_qtl_data } from './subworkflows/qtls.nf'
include { preprocess_genetic_data } from './subworkflows/bgen.nf'
include { generate_independent_snps } from './subworkflows/clump.nf'
include { check_interactions } from './subworkflows/interactions.nf'

workflow {
    preprocess_qtl_data()

    preprocess_genetic_data(preprocess_qtl_data.out.chrs)

    create_tf_bed_channel(preprocess_genetic_data.out,  preprocess_qtl_data.out.tf_chr_eqtls)

    generate_independent_snps(create_tf_bed_channel.out.tf_bed)

    check_interactions(preprocess_qtl_data.out.tf_chr_bqtls, create_tf_bed_channel.out.bed_files, generate_independent_snps.out)
}
