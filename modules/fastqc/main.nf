#!/usr/bin/env nextflow

process FASTQC {
    container "ghcr.io/bf528/fastqc:latest"
    label "process_single"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(name), path(fastqc)
    
    output:
    tuple val(name), path("*.html"), emit: html
    tuple val(name), path("*.zip"), emit: zip

    shell:
    """
    fastqc $fastqc
    """
}