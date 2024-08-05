import pandas as pd
import numpy as np
import pathlib
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import os
import argparse
import pyfaidx
import sys

##### helper functions for processing input BAM files
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

##### helper functions for processing CpG status in reads#
def get_refseq(path_to_all_fa, chrom, start, end):
    """
    Retrieves the reference sequence from a given FASTA file.
    Args:
        path_to_all_fa (str): The path to the directory containing all the FASTA files.
        chrom (str): The chromosome identifier.
        start (int): The starting position of the sequence.
        end (int): The ending position of the sequence.
    Returns:
        str: The uppercase reference sequence.
    Raises:
        FileNotFoundError: If the FASTA file for the specified chromosome is not found.
    """
    refseq = pyfaidx.Fasta(os.path.join(path_to_all_fa, "chr{}.fa".format(chrom)))
    return(str.upper(refseq.get_seq(name = "chr{}".format(chrom), start = start, end = end).seq))
            
def get_CpG_status(read_start, read_end, read, cpg_pos, mode = "string"):    
    """
    Returns the CpG status of a given genomic position within a read.

    Parameters:
    - read_start (int): The starting position of the read.
    - read_end (int): The ending position of the read.
    - read (str): The sequence of the read.
    - cpg_pos (int): The genomic position of interest.
    - mode (str): The mode of the output. Can be "string" or "num". Defaults to "string".

    Returns:
    - str or int: The CpG status of the genomic position. If mode is "string", returns the sequence at the position. If mode is "num", returns -1 for "not_covered", 1 for "CG", and 0 for other sequences.
    """
    if (read_start <= cpg_pos) and (read_end >= cpg_pos):
        seq = read[cpg_pos - read_start: cpg_pos + 1 - read_start + 1]
    else:
        seq = "not_covered"
    if mode == "string":
        return seq
    elif mode == "num":
        if seq == "not_covered":
            seq = -1
        elif seq == "CG":
            seq = 1
        else:
            seq = 0
        return seq
