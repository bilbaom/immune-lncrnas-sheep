# immune-lncrnas-sheep
This repository contains the pipeline and custom scripts used in the following manuscript:

"Comprehensive analysis of ovine transcriptomic data reveals novel long non-coding RNAs related to the immune response"

Each folder contains specific code for every analysis:

## RNA-seq analysis pipeline ([pipeline_rnaseq](pipeline_rnaseq))

Snakemake pipeline configured for a HPC server with Slurm. Conda environments are used for reproducible package version. This workflow follows the following steps: Download of RNA-seq data from NCBI SRA, sample preprocessing, mapping to the genome, trancriptome assembly and expression quantification. The workflow is divided into two Snakemake files because of storage memory limit but can be merged if there is enougth storage for all downloaded data and intermediate files.

## LncRNA identification pipeline ([pipeline_lncrnas](pipeline_lncrnas))

Snakemake pipeline for the detection of unannotated lncRNAs from the transcriptome assembly. Packages needed: gffcompare, gffread, HMMER, CPAT and FEELnc.

It uses the following custom scripts:

* [filter2.py](/pipeline_lncrnas/scripts/filter2.py): Executable script to filter the gffcompare output for novel potential lncRNAs.

`filter2.py [gffcompare.gtf] [outdir]`

* [translate3frames.py](/pipeline_lncrnas/scripts/translate3frames.py): Executable script to translate sequences in the 3 possible frames.

`translate3frames.py [dna_sequences.fa] [output.fa]`

* [classification_lncrnas.py](/pipeline_lncrnas/scripts/classification_lncrnas.py): Executable script to classify lncRNAs by location to nearest gene. The [pybedtools](https://github.com/daler/pybedtools) suite is necessary.

`classification_lncrnas.py [lncrna_list.txt] [gffcompare.gtf] [ensembl.gtf] [outdir]`

## CAGE-seq data analysis ([cage](cage))

Scripts for obtaining CAGE peaks from downloaded preprocessed data (modified from Salavati et al., 2020)

## CHIP-seq data analysis ([chip](chip))

Script for downloading raw data and obtaining histone CHIP-seq peaks.

## Differential expression and co-expression analyses

* [diff_expr.R](/diff_expr.R): Differential expression analysis
* [coex.R](/coex.R): Co-expression analysis and differential co-expresion network analysis
