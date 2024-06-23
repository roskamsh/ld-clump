#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.R2_THRESHOLD = 0.8
params.INFO_THRESHOLD = 0.9

include { create_input_channels; create_tf_bed_channel } from './subworkflows/channels.nf' 
include { preprocess_snps } from './subworkflows/snps.nf'
include { generate_independent_snps } from './subworkflows/clump.nf'

workflow {
    create_input_channels()

    preprocess_snps(create_input_channels.out.chr_bgen)

    create_tf_bed_channel(preprocess_snps.out, create_input_channels.out.tfs_bgen, create_input_channels.out.chr_bgen)

    generate_independent_snps(create_tf_bed_channel.out)
}
