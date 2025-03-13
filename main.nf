#!/usr/bin/env nextflow

include { FASTQC } from './modules/fastqc'
include { STAR_INDEX } from './modules/star_index'
include { GTF_PARSE } from './modules/gtf_parse'
include { STAR_ALIGNMENT } from './modules/star_alignment'
include { MULTIQC } from './modules/multiqc'
include { VERSE } from './modules/verse'
include { CONCAT } from './modules/concat'

workflow {

    Channel.fromFilePairs(params.reads)
    | set { align_ch }

    Channel.fromFilePairs(params.reads).transpose()
    | set { fastqc_ch }

    GTF_PARSE(params.gtf)
    FASTQC(fastqc_ch)
    STAR_INDEX(params.genome, params.gtf)
    STAR_ALIGNMENT(align_ch, STAR_INDEX.out.index)

    FASTQC.out.zip.map { it[1] }.collect()
    | set { fastqc_out }

    STAR_ALIGNMENT.out.log.map { it[1] }.collect()
    | set { star_log } // might be different name

    fastqc_out.mix(star_log).flatten().collect()
    | set { multiqc_ch }

    MULTIQC(multiqc_ch)

    VERSE(STAR_ALIGNMENT.out.bam, params.gtf)
    
    VERSE.out.counts.map{ it[1] }.collect().set{ concat_ch }
    CONCAT(concat_ch)


}
