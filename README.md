# immune-lncrnas-sheep
This repository contains the pipeline and custom scripts used in the following manuscript:

"Comprehensive analysis of ovine transcriptomic data reveals novel long non-coding RNAs related to the immune response"

Each folder contains specific code for every analysis:

## RNA-seq analysis pipeline ([pipeline_rnaseq](pipeline_rnaseq))

Snakemake pipeline configured for a HPC server with Slurm. Conda environments are used for reproducible package version. This workflow follows the following steps: Download of RNA-seq data from NCBI SRA, sample preprocessing, mapping to the genome, trancriptome assembly and expression quantification. The workflow is divided into two Snakemake files because of storage memory limit but can be merged if there is enougth storage for all downloaded data and intermeiate files.

## LncRNA identification pipeline ([pipeline_lncrnas](pipeline_lncrnas))

Snakemake pipeline for the detection of unannotated lncRNAs from the transcriptome assembly. It uses custom scripts.

## CAGE-seq data analysis ([cage](cage))

Script for obtaining CAGE peaks from downloaded data.

## CHIP-seq data analysis ([chip](chip))

Script for obtaining histone CHIP-seq peaks from downloaded data.

## Differential expression and co-expression analyses ([R_analyses](R_analyses))

R scripts for the different functional analyses: Differential expression analysis, co-expression analysis and differential co-expresion network analysis.
