# BF528 Project 1: RNAseq Analysis Nextflow Pipeline

## Overview
In this repo, a modularized Nextflow pipeline that performs read quality control, genome alignment, and gene quantification. It integrates tools such as FastQC, STAR, VERSE, and MultiQC for end-to-end processing—from raw FASTQ files to gene count matrices. Modular process design ensures reproducibility and easy customization across different reference genomes and GTF annotations. The final output includes per-sample QC reports, alignment logs, and a combined gene expression matrix ready for downstream analysis. 
