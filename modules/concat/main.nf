process CONCAT {
    label 'process_low'
    container 'ghcr.io/bf528/pandas:latest'
    publishDir params.outdir

    input:
    path df_list

    output:
    path 'verse_concat.csv'
    
    shell:
    """
    concat_df.py -i ${df_list.join(' ')} -o verse_concat.csv
    """
}