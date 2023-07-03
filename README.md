# README

This repository contains code and methods for the publication:

***Comprehensive saturation genome editing of RAD51C provides precise functional classification within a spectrum of variant effects***

<span style="color: red;">**Please note: These data are unpublished and not yet subject to peer review. We provide them as a service to the research community but they are embargoed until publication of our paper. They should not be used as the sole basis for clinical decision making.**</span>

***

## rad51c_analysis_code.R

This is code written in [R](https://cran.r-project.org/) to convert next generation sequencing counts into functional scores for RAD51C SGE. 

Broad functions for `rad51c_analysis_code.RR` are:

* import SGE counts, which are held as separate files for each sample output from the [QUANTS](https://github.com/cancerit/QUANTS) pipeline, into count matrices for each target region
* annotate variants with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) outputs and [VaLiAnT](https://github.com/cancerit/VaLiAnT) meta-data outputs
* filter counts and create normalisation matrices (synonymous and intronic variants)
* run [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) with adjusted size factors to calculate raw log fold changes (LFCs) and statistics
* produce QC plots for replicate consistency, dispersion, sample correlation and editing efficiency
* median scale and combine overlapping target region LFCs to generate a single LFC
* calculate adjusted z-scores, standard error and p-value
* functional classification of fast and slow depletion based on Gaussian Mixture Modelling
* FDR correction of p-value
* data frame cleanup and annotation refinement

*Note: PATHS are currently (March 2023) hard coded, with specific steps for RAD51C SGE annotation only.*

## sge_rad51c.csv

A data frame containing 9,188 unique nucleotide-levels within RAD51C with functional scores, functional classifications and annotation fields.

## dataframe_columns.xlsm
A table containing descriptions of the each of the fields in [sge_rad51c.csv](sge_rad51c.csv).
