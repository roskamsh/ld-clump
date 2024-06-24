#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.R2_THRESHOLD = 0.8
params.INFO_THRESHOLD = 0.9
params.ASSEMBLY = "hg19"
params.eQTLGEN_DATA = ""
params.Concordant = "True"

include { create_input_channels; create_tf_bed_channel } from './subworkflows/channels.nf' 
include { preprocess_qtl_data } from './subworkflows/qtls.nf'
include { preprocess_genetic_data } from './subworkflows/bgen.nf'
include { generate_independent_snps } from './subworkflows/clump.nf'

workflow {
    create_input_channels()

    preprocess_qtl_data(create_input_channels.out.tfs)

    preprocess_genetic_data(create_input_channels.out.chr_bgen)

    create_tf_bed_channel(preprocess_genetic_data.out,  preprocess_qtl_data.out.tfs_eqtls, create_input_channels.out.chr_bgen)

    generate_independent_snps(create_tf_bed_channel.out)
}
