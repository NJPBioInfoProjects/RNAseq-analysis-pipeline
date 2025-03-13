#!/usr/bin/env nextflow

process STAR_ALIGNMENT {
    container 'ghcr.io/bf528/star:latest'
    label 'process_high'

    input:
    tuple val(name), path(fastq) //fastq is the reads
    path index

    output:
    tuple val(name), path ("${name}.Aligned.out.bam"), emit: bam
    tuple val(name), path ("${name}.Log.final.out"), emit: log

    

    shell:

    def reads_str = fastq.join(" ")

    """
    STAR --runThreadN $task.cpus --genomeDir $index --readFilesIn ${reads_str} --readFilesCommand gunzip -c --outFileNamePrefix ${name}. --outSAMtype BAM Unsorted
    """
}


// --runThreadN NumberOfThreads  
// --genomeDir /path/to/genomeDir  
// --readFilesIn /path/to/read1 [/path/to/read2]  

// --genomeDir specifies the path to the genome directory where genome 
// indices were generated (see Section 2. Generating genome indexes).

// --readFilesIn: name(s) (with path) of the files containing the sequences to be mapped (e.g. RNA-seq FASTQ files).  
// If using Illumina paired-end reads, the read1 and read2 files have to be supplied.  
// STAR can process both FASTA and FASTQ files.  
// Multi-line (i.e. sequence split in multiple lines) FASTA (but not FASTQ) files are supported.  

// --readFilesCommand: If the read files are compressed, use the --readFilesCommand UncompressionCommand option,  
// where UncompressionCommand is the un-compression command that takes the file name as an input parameter,  
// and sends the uncompressed output to stdout.  
// For example, for gzipped files (*.gz) use --readFilesCommand zcat OR --readFilesCommand gunzip -c.  
// For bzip2-compressed files, use --readFilesCommand bunzip2 -c.  

// --outFileNamePrefix: STAR produces multiple output files. All files have a standard name,  
// however, you can change the file prefixes using --outFileNamePrefix /path/to/output/dir/prefix.  
// By default, this parameter is ./, i.e., all output files are written in the current directory.  

// 5.3 Unsorted and sorted-by-coordinate BAM.  
// STAR can output alignments directly in binary BAM format,  
// thus saving time on converting SAM files to BAM.  
// It can also sort BAM files by coordinates, which is required by many downstream applications.  

// --outSAMtype BAM Unsorted  
// Output unsorted Aligned.out.bam file.  
// The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well.  
// This "unsorted" file can be directly used with downstream software such as HTseq,  
// without the need for name sorting.  
// The order of the reads will match that of the input FASTQ(A) files only if one thread is used  
// --runThread 1, and --outFilterType --BySJout is not used.  

// --outSAMtype BAM SortedByCoordinate  
// Output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.  
// If this option causes problems, it is recommended to reduce --outBAMsortingThreadN  
// from the default 6 to lower values (as low as 1).  

// --outSAMtype BAM Unsorted SortedByCoordinate  
// Output both unsorted and sorted files.  