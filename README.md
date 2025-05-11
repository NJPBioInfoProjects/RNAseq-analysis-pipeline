# BF528 Project 1: RNAseq Analysis Nextflow Pipeline

## Overview
In this repo, a modularized Nextflow pipeline that performs read quality control, genome alignment, and gene quantification. It integrates tools such as FastQC, STAR, VERSE, and MultiQC for end-to-end processing—from raw FASTQ files to gene count matrices. Modular process design ensures reproducibility and easy customization across different reference genomes and GTF annotations. The final output includes per-sample QC reports, alignment logs, and a combined gene expression matrix ready for downstream analysis. 

Downstream analysis is conducted in analysis.Rmd, an R Markdown notebook that performs differential expression analysis using DESeq2, followed by functional enrichment with DAVID and fgsea. The notebook includes quality assessment, PCA, heatmaps, and enrichment plots to help interpret transcriptional changes across conditions. It also compares results with published findings and highlights key biological processes impacted by gene expression changes.

## Usage
To run the nextflow pipeline:
1. Create the conda environment
  - Ensure miniconda module is loaded and run:
   ```bash
   conda env create -f envs/base_env.yml
   ```
  - When complete run:
   ```bash
   conda activate nextflow_base
   ```
2. Run the pipeline:
  - Run:
   ```bash
   nextflow run main.nf -profile singularity,cluster
   ```
