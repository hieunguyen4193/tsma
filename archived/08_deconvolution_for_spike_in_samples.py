import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm
import pathlib
import os

from helper_functions import *
metadata = pd.read_csv("metadata.csv")
metadata_cfdna = metadata[metadata["Sample type"].isin(["cfDNA", "Control"])]
import warnings
warnings.filterwarnings("ignore")
import random

for topK in [500, 1000]:
    # for atlas_sample_types in ["Tissue", "Tissue,WBC", "Tissue,WBC,Control", "Tissue,Control"]:
    for atlas_sample_types in ["Tissue,WBC", "Tissue,WBC,Control", "Tissue,Control"]:
        datadir = "./datadir"
        outputdir = "./outputdir_02102023"
        
        path_to_03_output = os.path.join(outputdir, "03_output_noFDR")
        path_to_04_output = os.path.join(outputdir, "04_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
        path_to_05_output = os.path.join(outputdir, "05_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
        path_to_06_output = os.path.join(outputdir, "06_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
        path_to_07_output = os.path.join(outputdir, "07_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
        path_to_08_output = os.path.join(outputdir, "08_output_noFDR", "top{}_{}".format(topK, atlas_sample_types.replace(",", "_and_")))
        os.system("mkdir -p {}".format(path_to_08_output))
        os.system("mkdir -p {}".format(os.path.join(path_to_08_output, "tmp")))
        atlas = pd.read_csv(os.path.join(path_to_03_output, "top{}_atlas_{}.final.csv".format(topK, atlas_sample_types)), index_col =[0])
        
        atlas = atlas[[item for item in atlas.columns if "_y" not in item ]]
        atlas.columns = [item.replace("_x", "") for item in atlas.columns]
        
        avg_atlas = atlas.set_index("sample").fillna(0).groupby("label").mean()
        avg_atlas = avg_atlas.loc[atlas_labels[atlas_sample_types]]
        atlas_regions = [item for item in atlas.columns if item not in ["sample", "label"]]
        
        
        resdf = pd.DataFrame()
        all_control_samples = [item for item in pathlib.Path(path_to_07_output).glob("*")]
        
        for sample in all_control_samples:
            selected_metadata = metadata[metadata["filename"] == sample.name]
            for spike_in_rate in [0.01, 0.05, 0.1, 0.15, 0.2, 0.3]:
                for label in ["Liver", "Gastric", "Breast", "Lung", "CRC"]:
                    all_spike_in_samples = [item for item in pathlib.Path(
                        os.path.join(str(sample), "spike_in_rate_{}".format(spike_in_rate), label)).glob("*.csv")]
                    for file in all_spike_in_samples[0:4]:
                        rep = file.name.split("_")[-1].strip(".csv")     
                        if os.path.exists(os.path.join(path_to_08_output, "tmp", "deconvolution_results_{}_{}_{}_{}.csv".format(sample.name, label, spike_in_rate, rep))) == False:
                            dffull = pd.read_csv(file, index_col = [0])
                            all_betadf = pd.DataFrame(data = [item for item in selected_metadata["filename"].to_list()], columns = ["sample"])
                                            
                            for region in tqdm(atlas_regions):
                                df = dffull[dffull["region"] == "{}:{}-{}".format(region.split("_")[0],
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
                                
                                df["label"] = df["sample"].apply(lambda x: metadata[metadata["filename"] == x]["Label"].unique()[0])
                                
                                df["check_cigar"] = df["cigar"].apply(lambda x: bool(pattern.fullmatch(x)))
                                df = df[df["check_cigar"] == True]
                                
                                df["Label"] = df["sample"].apply(lambda x: metadata[metadata["filename"] == x]["Label"].unique()[0])
                                
                                if df.shape[0] != 0:
                                    df["end"] = df[["start", "cigar"]].apply(lambda x: int(x[0]) + int(x[1].replace("M", "")), axis = 1)
                                    for cpg_pos in cpg_coords:
                                        df[cpg_pos] = df[["start", "end", "seq"]].apply(lambda x: get_CpG_status(x[0], x[1], x[2], cpg_pos, mode = "num"), axis = 1)
                                    
                                    betadf = selected_metadata[["filename"]].copy()
                                    betadf.columns = ["sample"]
                                    
                                    for cpg_pos in cpg_coords:    
                                        tmpdf = df[["sample", cpg_pos, "Label"]].copy()
                                        tmpcountdf = tmpdf.groupby('sample')[cpg_pos].apply(lambda x: (x == 1).sum()/((x == 0).sum() + (x == 1).sum()) ).reset_index(name= "meth_level_{}".format(cpg_pos))
                                        betadf = betadf.merge(tmpcountdf[["sample", "meth_level_{}".format(cpg_pos)]], right_on = "sample", left_on = "sample", how = "outer")
                                    
                                    betadf["avg_beta"] = betadf[[item for item in betadf.columns if item != "sample"]].apply(lambda x: np.mean([item for item in x if np.isnan(item) == False]), axis = 1)
                                    betadf["label"] = betadf["sample"].apply(lambda x: metadata[metadata["filename"] == x]["Label"].unique()[0])
                                    betadf = betadf[["sample", "avg_beta"]]
                                    betadf.columns = ["sample", region]
                                    all_betadf = all_betadf.merge(betadf[["sample", region]], right_on = "sample", left_on = "sample")
                            deconvo_res_cfdna = deconvo(all_betadf, avg_atlas, atlas_sample_types)
                            tmpresdf = deconvo_res_cfdna.set_index("TOO").T.copy()
                            tmpresdf["spike_in_rate"] = spike_in_rate
                            tmpresdf["spike_in_sample"] = label
                            # rep = file.name.split("_")[-1].strip(".csv")
                            tmpresdf["rep"] = rep
                            resdf = pd.concat([resdf, tmpresdf], axis = 0)
                            tmpresdf.to_csv(os.path.join(path_to_08_output, "tmp", "deconvolution_results_{}_{}_{}_{}.csv".format(sample.name, label, spike_in_rate, rep)))
                        else:
                            tmpresdf = pd.read_csv(os.path.join(path_to_08_output, "tmp", "deconvolution_results_{}_{}_{}_{}.csv".format(sample.name, label, spike_in_rate, rep)), index_col = [0])
        resdf.to_csv(os.path.join(path_to_08_output, "deconvolution_results_for_spike_in_samples.csv"))
