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

'''
Use this function to collect reads at pre-defined CpG regions. 
There are three options for the sample_type: Tissue, WBC or cfDNA. 
See the attached metadata file for more details on sample_type.
'''
##### EXAMPLE INPUT
# bamdir = "/mnt/DATASM14/hieunho/hieu_data/GW-bismark_cfDNA-highdepth_pair-end/04_bismark_deduplicated_filtered_Q30"

def main(args):
    ##### INPUT ARGS
    region = args.region
    outputdir = args.output
    path_to_metadata = args.metadata
    bamdir = args.bamdir
    sample_type = args.sample_type
    os.system("mkdir -p {}".format(outputdir))
    
    ##### Assertion
    assert sample_type in ["Tissue", "cfDNA", "WBC", "Control"]
    region_name = region.replace(":", "_").replace("-", "_")
    ##### MAIN
    if os.path.isfile(os.path.join(outputdir, "read_at_region_{}.csv".format(region_name))) == False:
        metadata = pd.read_csv(path_to_metadata)    
        all_bam_files = [item for item in pathlib.Path(bamdir).glob("*.bam")]
        chosen_bam_files = [item for item in all_bam_files if item.name.split(".")[0] in metadata[metadata["Sample type"] == sample_type].filename.unique()]
        output_readdf = pd.DataFrame()
        for bamfile in chosen_bam_files:
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
            readdf["sample"] = bamfile.name.split(".")[0]
            readdf["region"] = region
            output_readdf = pd.concat([output_readdf, readdf], axis = 0)
        output_readdf.to_csv(os.path.join(outputdir, "read_at_region_{}.csv".format(region_name)))
    else:
        print("region {} is already in the output dir {} of sample_type {}".format(region, outputdir, sample_type))
        
def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--region', type=str,
        help='names of the region to collect reads from BAM file', action='store')
    parser.add_argument('--output', type=str, action = "store",
        help='Path to output file')
    parser.add_argument('--metadata', type=str, action = "store",
        help='Path to metadata file')
    parser.add_argument('--bamdir', type=str, action = "store",
        help='Path to BAM directory')
    parser.add_argument('--sample_type', type=str, action = "store",
        help='Choose the sample type, this can only be Tissue, WBC or cfDNA')
    
    return parser.parse_args(argv)

if __name__ == '__main__':
    main(parse_arguments(sys.argv[1:]))