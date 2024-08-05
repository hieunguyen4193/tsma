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
import re

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
    readdf["sample"] = str(bamfile).split("/")[-1].split(".")[0]
    readdf["region"] = region
    return readdf

def fetch_reads_multiple_bams(all_bam_files, region):
    """
    Fetches reads from multiple BAM files at a specific region.
    Parameters:
    - all_bam_files (str): List of path to bam files. Each bam file should be sorted and indexed with `samtools index bam.file`
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
    
    # all_bam_files = [item for item in pathlib.Path(bamfolder).glob("*.bam")]
    output_readdf = pd.DataFrame()
    for bamfile in tqdm(all_bam_files):
        print("Fetching reads from file {}".format(str(bamfile).split("/")[-1]))
        readdf = fetch_reads(bamfile, region)
        output_readdf = pd.concat([output_readdf, readdf], axis = 0)
    return output_readdf

def generate_betadf(bamfiles, region, path_to_all_fa, outputdir, save_file = False):
    # bamdir = "/mnt/archiving/DATA_HIEUNHO/GW-bismark_cfDNA-highdepth_pair-end/04_bismark_deduplicated_sorted"
    # region = "1:1103264-1103363"
    # path_to_all_fa = "/datassd/hieunguyen/hg19"
    
    # metadata = pd.read_csv("/mnt/DATASM14/hieunho/hieu_project/metadata/ECD/metadata_GW_highdepth.csv")
    # bamfiles = [item for item in pathlib.Path(bamdir).glob("*.bam") if str(item) in metadata.Bam_file.unique()]
    # print("Number of bam files in this location: {}".format(len(bamfiles)))    
    # outputdir = "./outputs"
    
    os.system("mkdir -p {}".format(outputdir))
    
    output_readdf = fetch_reads_multiple_bams(bamfiles, region)

    # save the file
    if save_file:
        output_readdf.to_csv(os.path.join(outputdir, "{}.all_samples.csv".format(region)))
    
    region = region.replace(":", "_").replace("-", "_")
        
    region_chrom = region.split("_")[0]
    region_start = int(region.split("_")[1])
    region_end = int(region.split("_")[2])
    
    refseq_at_cluster = get_refseq(path_to_all_fa = path_to_all_fa, 
                                    chrom = region_chrom, 
                                    start = region_start, 
                                    end = region_end + 1)
    all_cpg_in_cluster = [m.start(0) for m in re.finditer("CG", refseq_at_cluster)]
    cpg_coords = [item + region_start for item in all_cpg_in_cluster]
    
    df = output_readdf.copy()
    
    pattern = re.compile("^[1-9][0-9]M")
    df["check_cigar"] = df["cigar"].apply(lambda x: bool(pattern.fullmatch(x)))
    df = df[df["check_cigar"] == True]
                
    df["end"] = df[["start", "cigar"]].apply(lambda x: int(x[0]) + int(x[1].replace("M", "")), axis = 1)
    
    for cpg_pos in cpg_coords:
        df[cpg_pos] = df[["start", "end", "seq"]].apply(lambda x: get_CpG_status(x[0], x[1], x[2], cpg_pos, mode = "num"), axis = 1)
        
    betadf = pd.DataFrame(data = df["sample"].unique())
    betadf.columns = ["sample"]
    for cpg_pos in tqdm(cpg_coords):    
        tmpdf = df[["sample", cpg_pos]].copy()
        tmpcountdf = tmpdf.groupby('sample')[cpg_pos].apply(lambda x: (x == 1).sum()/((x == 0).sum() + (x == 1).sum()) ).reset_index(name= "meth_level_{}".format(cpg_pos))
        betadf = betadf.merge(tmpcountdf[["sample", "meth_level_{}".format(cpg_pos)]], right_on = "sample", left_on = "sample", how = "outer")
            
    betadf["avg_beta"] = betadf[[item for item in betadf.columns if item != "sample"]].apply(lambda x: np.mean([item for item in x if np.isnan(item) == False]), axis = 1)
    if save_file:
        betadf.to_csv(os.path.join(outputdir, "{}_beta_values.all_samples.csv".format(region)))
    return betadf

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
    read_start = int(read_start)
    read_end  = int(read_end)
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
