# BF528 Project 1: RNAseq Analysis Nextflow Pipeline

## Overview
In this repo, a modularized Nextflow pipeline that performs read quality control, genome alignment, and gene quantification. It integrates tools such as FastQC, STAR, VERSE, and MultiQC for end-to-end processing—from raw FASTQ files to gene count matrices. Modular process design ensures reproducibility and easy customization across different reference genomes and GTF annotations. The final output includes per-sample QC reports, alignment logs, and a combined gene expression matrix ready for downstream analysis. 

Downstream analysis is conducted in analysis.Rmd, an R Markdown notebook that performs differential expression analysis using DESeq2, followed by functional enrichment with DAVID and fgsea. The notebook includes quality assessment, PCA, heatmaps, and enrichment plots to help interpret transcriptional changes across conditions. It also compares results with published findings and highlights key biological processes impacted by gene expression changes.

The data used in this particular analysis came from Chandra et al., 2022. 

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
   nextflow run main.nf -profile singularity,local
   ```

## Repo Structure:
- bin/: Contains scripts used for custom Nextflow modules
- envs/: Contains yml file for creating the nextflow_base conda environment that is required to the run the pipeline
- modules/: Contains Nextflow modules for each pipeline step, including FastQC, STAR indexing and alignment, GTF parsing, read quantification with VERSE, MultiQC aggregation, and count matrix concatenation
- analysis.rmd: R Markdown notebook used for downstream analysis, including differential expression with DESeq2, quality assessment, PCA, heatmaps, and functional enrichment using DAVID and fgsea
- main.nf: Main Nextflow pipeline that runs all modules contained in /modules
- nextflow.config: Configuration file that defines input paths, output directories, and reference files used in the pipeline. It also specifies execution profiles for running the workflow locally, with Conda/Singularity, or on an SGE cluster (e.g., SCC), including custom resource labels for different process intensities. The config ensures reproducibility and flexibility across computing environments
- *.csv / *.txt: Output files from downstream analysis, including normalized counts, a list of significant genes, and the top 10 DE genes with gene names
- *.png: Plots exported from the analysis notebook, such as FastQC and GSEA results, used for quick visual summaries
- c2.cp.v2024.1.Hs.symbols.gmt: Gene set file used for pathway enrichment analysis with fgsea

## References
- Chandra, V., Ibrahim, H., Halliez, C. et al. The type 1 diabetes gene TYK2 regulates β-cell development and its responses to interferon-α. Nat Commun 13, 6363 (2022). https://doi.org/10.1038/s41467-022-34069-z
