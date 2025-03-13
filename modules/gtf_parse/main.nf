#!/usr/bin/env nextflow

process GTF_PARSE {
    label 'process_single'
    container 'ghcr.io/bf528/biopython:latest'
    publishDir params.outdir

    input:
    path gtf

    output:
    path('id2name.txt'), emit: id2name

    shell:
    """
    gtf_parser.py -i $gtf -o id2name.txt
    """
}