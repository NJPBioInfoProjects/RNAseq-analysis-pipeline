process VERSE{
    label 'process_high'
    container'ghcr.io/bf528/verse:latest'
    publishDir params.outdir
    
    input:
    tuple val(name), path(bam)
    path(gtf)

    output:
    tuple val(name), path("*exon.txt"), emit: counts

    shell:
    """
    verse -S -a $gtf -o $name $bam
    """
}