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
from helper_functions import *

# outputdir = "./outputdir_02102023"
outputdir = "./outputdir_23102023_hypo_hyper"
topK = 500
atlas_sample_types = "Tissue,WBC"

path_to_03_output = os.path.join(outputdir, "03_output_noFDR")

path_to_10_output = os.path.join(outputdir, "10_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
os.system("mkdir -p {}".format(path_to_10_output))

atlas = pd.read_csv(os.path.join(path_to_03_output, "top{}_atlas_{}.final.csv".format(topK, atlas_sample_types)), index_col =[0])

atlas = atlas[[item for item in atlas.columns if "_y" not in item ]]
atlas.columns = [item.replace("_x", "") for item in atlas.columns]

avg_atlas = atlas.set_index("sample").fillna(0).groupby("label").mean()
avg_atlas = avg_atlas.loc[atlas_labels[atlas_sample_types]]
atlas_regions = [item for item in atlas.columns if item not in ["sample", "label"]]


bamdir = os.path.join("outputdir_02102023", "cLOD_samples", "BAM", "sorted_bam")
all_bam_files = [item for item in pathlib.Path(bamdir).glob("*.sorted.bam")]

# bamfile = os.path.join(bamdir, "1-LODCRC20R2CT496W_M561-M761.deduplicated.sorted.bam")
for bamfile in tqdm(all_bam_files):
    # samplename = bamfile.split("/")[-1]
    samplename = bamfile.name
    
    if (os.path.isfile(os.path.join(path_to_10_output, "Sample_{}.deconvo.csv".format(samplename))) == False):
        output_readdf = pd.DataFrame()
        for region in tqdm(atlas_regions):
            region = "{}:{}-{}".format(region.split("_")[0], region.split("_")[1], region.split("_")[2])
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
            readdf["sample"] = samplename
            readdf["region"] = region
            output_readdf = pd.concat([output_readdf, readdf], axis = 0)
        output_readdf.to_csv(os.path.join(path_to_10_output, "Sample_{}.csv".format(samplename)))
        
        
        all_betadf = pd.DataFrame(data = [samplename], columns = ["sample"])
        for region in tqdm(atlas_regions):
            df = output_readdf[output_readdf["region"] == "{}:{}-{}".format(region.split("_")[0],
                                                             region.split("_")[1],
                                                             region.split("_")[2])]
            region_chrom = region.split("_")[0]
            region_start = int(region.split("_")[1])
            region_end = int(region.split("_")[2])
            refseq_at_cluster = get_refseq(path_to_all_fa = path_to_all_fa, 
                                                chrom = region_chrom, 
                                                start = region_start, 
                                                end = region_end + 1)
            all_cpg_in_cluster = [m.start(0) for m in re.finditer("CG", refseq_at_cluster)]
            cpg_coords = [item + region_start for item in all_cpg_in_cluster]
            
            df["check_cigar"] = df["cigar"].apply(lambda x: bool(pattern.fullmatch(x)))
            df = df[df["check_cigar"] == True]
            
            if df.shape[0] != 0:
                df["end"] = df[["start", "cigar"]].apply(lambda x: int(x[0]) + int(x[1].replace("M", "")), axis = 1)
                df["start"] = df["start"].astype(int)
                df["end"] = df["end"].astype(int)
                for cpg_pos in cpg_coords:
                    df[cpg_pos] = df[["start", "end", "seq"]].apply(lambda x: get_CpG_status(x[0], x[1], x[2], cpg_pos, mode = "num"), axis = 1)
                
                betadf = pd.DataFrame(data = [samplename], columns = ["sample"])
                
                for cpg_pos in cpg_coords:    
                    tmpdf = df[["sample", cpg_pos]].copy()
                    tmpcountdf = tmpdf.groupby('sample')[cpg_pos].apply(lambda x: (x == 1).sum()/((x == 0).sum() + (x == 1).sum()) ).reset_index(name= "meth_level_{}".format(cpg_pos))
                    betadf = betadf.merge(tmpcountdf[["sample", "meth_level_{}".format(cpg_pos)]], right_on = "sample", left_on = "sample", how = "outer")
                
                betadf["avg_beta"] = betadf[[item for item in betadf.columns if item != "sample"]].apply(lambda x: np.mean([item for item in x if np.isnan(item) == False]), axis = 1)
                betadf = betadf[["sample", "avg_beta"]]
                betadf.columns = ["sample", region]
                all_betadf = all_betadf.merge(betadf[["sample", region]], right_on = "sample", left_on = "sample")
        deconvo_res_cfdna = deconvo(all_betadf, avg_atlas, atlas_sample_types)
        
        deconvo_res_cfdna.to_csv(os.path.join(path_to_10_output, "Sample_{}.deconvo.csv".format(samplename)))