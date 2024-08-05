import pandas as pd
import numpy as np
import pathlib
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import os
import argparse
import sys

def fetch_reads(bamfile, region):
    """
    Fetches reads from a BAM file at a specific region.
    Parameters:
    - bamfile (str): Path to the BAM file.
    - region (str): Genomic region to fetch reads from. This should be taken from the file CpG_clusters_whole_genome_radius_100.bed
    Returns:
    - readdf (pandas.DataFrame): DataFrame containing the fetched reads with the following columns:
        - chrom: Chromosome name
        - start: Start position of the read
        - cigar: CIGAR string
        - flen: Length of the read
        - seq: Sequence of the read
        - methyl_string: Methyl string
        - XR: XR tag value
        - XG: XG tag value
        - sample: Sample name
        - region: Genomic region
    Example: 
    bamfile = "/Users/hieunguyen/data/bam_files/highdepth_cancer_WGBS_bismark.bam"
    region = "1:1103264-1103363"
    """
    # Function implementation goes here
    pass
    
    bamfile_obj = pysam.AlignmentFile(bamfile).fetch(region = region)
    reads = []
    for read in bamfile_obj:
        reads.append(read)
    readdf = pd.DataFrame()
    readdf["chrom"] = [read.to_dict()["ref_name"] for read in reads]
    readdf["start"] = [read.to_dict()["ref_pos"] for read in reads]
    readdf["cigar"] = [read.to_dict()["cigar"] for read in reads]
    readdf["flen"] = [read.to_dict()["length"] for read in reads]
    readdf["seq"] = [read.to_dict()["seq"] for read in reads]
    readdf["methyl_string"] = [read.to_dict()["tags"][2] for read in reads]
    readdf["XR"] = [read.to_dict()["tags"][3] for read in reads]
    readdf["XG"] = [read.to_dict()["tags"][4] for read in reads]
    readdf["sample"] = str(bamfile).split("/")[-1].split(".")[1]
    readdf["region"] = region
    return readdf

def fetch_reads_multiple_bams(bamfolder, region):
    """
    Fetches reads from multiple BAM files at a specific region.
    Parameters:
    - bamfolder (str): Path to the folder containing BAM files. Each bam file should be sorted and indexed with `samtools index bam.file`
    - region (str): Genomic region to fetch reads from. This should be taken from the file CpG_clusters_whole_genome_radius_100.bed
    Returns:
    - readdf (pandas.DataFrame): DataFrame containing the fetched reads with the following columns:
        - chrom: Chromosome name
        - start: Start position of the read
        - cigar: CIGAR string
        - flen: Length of the read
        - seq: Sequence of the read
        - methyl_string: Methyl string
        - XR: XR tag value
        - XG: XG tag value
        - sample: Sample name
        - region: Genomic region
    Example: 
    bamfolder = "/Users/hieunguyen/data/bam_files"
    region = "1:1103264-1103363"
    """
    # Function implementation goes here
    pass
    
    all_bam_files = [item for item in pathlib.Path(bamfolder).glob("*.bam")]
    output_readdf = pd.DataFrame()
    for bamfile in all_bam_files:
        readdf = fetch_reads(bamfile, region)
        output_readdf = pd.concat([output_readdf, readdf], axis = 0)
    return output_readdf

