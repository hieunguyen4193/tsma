# Tumor specific methylation atlas

Code repository for "Tissue of origin detection for cancer tumor using low-depth cfDNA samples through combination of tumor-specific methylation atlas and genome-wide methylation density in graph convolutional neural networks" (DOI https://doi.org/10.1186/s12967-024-05416-z).  

## Table of Contents

- [Tumor specific methylation atlas](#tumor-specific-methylation-atlas)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
  - [Usage](#usage)

## Requirements
Conda environment to run code in this repository: `conda env create -f env.yml`

## Usage
To replicate the analysis in the paper, proceed as follows

- Run `construct_CpG_regions.py` to generate a `BED` file containing more than 1.1m regions. Each region is a cluste of at least 5 CpGs in a 100bp window across the whole genome (hg19).

- Main functions are stored in `helper_functions.py`, 
  - `fetch_reads`: Fetch all reads overlapping a genomic region (`region`: `chr1:1-100`) from a `BAM` file.
  - `fetch_reads_multiple_bams`: Fetch reads using `fetch_reads` that overlap a genomic region in a list of `BAM` files. Output is stored in a `pandas` dataframe. See `./outputs/1:1103264-1103363.all_samples.csv` for an example. 
  - Once reads are collected from all input `BAM` files, run the function `generate_betadf` to calculate methylation regional values at the given region. See `./outputs/1_1103264_1103363_beta_values.all_samples.csv` for an example. 
- Combine the beta value dataframe with the metadata `metadata_GW_highdepth.csv`, assign samples to their corresponding labels and conduct the statistical t-test of your choice. 

- UPDATE ME: 
  - add script to convert TCGA methylation array to match the sequencing data format. 
  - add deconvolution functions.
  - add summary script to generate the final atlas. 
