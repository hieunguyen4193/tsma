import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
import sys
import pathlib
from tqdm import tqdm
import re
pattern = re.compile("^[1-9][0-9]M")
import warnings
warnings.filterwarnings("ignore")
from helper_functions import *

##### EXAMPLE OF INPUT ARGS
# path_to_metadata = "metadata.csv"
# metadata = pd.read_csv(path_to_metadata)
# region = "10:100992305-100992405"
# atlas_sample_types = "Tissue"
# path_to_main_output = "."    
# LOO_sample = "not_given"

def main(args):
    ##### INPUT ARGS
    region = args.region
    path_to_main_output = args.output
    path_to_metadata = args.metadata
    metadata = pd.read_csv(path_to_metadata)
    atlas_sample_types = args.atlas_sample_types
    LOO_sample = args.LOO_sample
    
    ##### MAIN
    path_to_02_output = os.path.join(path_to_main_output, "02_output", atlas_sample_types.replace(",", "_and_"))
    if (LOO_sample != "not_given" and LOO_sample in metadata.filename.unique()):
        path_to_02_output = os.path.join(path_to_main_output, "02_output_LOO", atlas_sample_types.replace(",", "_and_"))
    os.system("mkdir -p {}".format(path_to_02_output))

    summ_test_file = os.path.join(path_to_02_output, "summary_of_all_ttest_results.csv")
    if os.path.isfile(summ_test_file) == False:
        os.system("touch {}".format(summ_test_file))
        # os.system("echo {} >> {}".format(summ_file_columns[atlas_sample_types], summ_test_file))
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

    if os.path.exists(os.path.join(path_to_02_output, "ttest_res_tissue_and_wbc_samples_{}.csv".format(region))) == False:
        df = pd.DataFrame()
        if "," in atlas_sample_types:
            selected_sample_types = atlas_sample_types.split(",")
            selected_metadata = metadata[metadata["Sample type"].isin(selected_sample_types)]
            for sampletype in selected_sample_types:
                assert sampletype in ["Tissue", "WBC", "cfDNA", "Control"]
                tmpdf = pd.read_csv(os.path.join(readdf_dir[sampletype], "read_at_region_{}.csv".format(region)), index_col = [0])
                df = pd.concat([df, tmpdf], axis = 0)
        elif (atlas_sample_types == "cfDNA"):
            selected_sample_types = ["cfDNA", "Control"]
            selected_metadata = metadata[metadata["Sample type"].isin(selected_sample_types)]
            for sampletype in selected_sample_types:
                assert sampletype in ["Tissue", "WBC", "cfDNA", "Control"]
                tmpdf = pd.read_csv(os.path.join(readdf_dir[sampletype], "read_at_region_{}.csv".format(region)), index_col = [0])
                df = pd.concat([df, tmpdf], axis = 0)
        else:
            assert atlas_sample_types in ["Tissue", "WBC", "cfDNA", "Control"]
            df = pd.read_csv(os.path.join(readdf_dir[atlas_sample_types], "read_at_region_{}.csv".format(region)), index_col = [0])
            selected_metadata = metadata[metadata["Sample type"] == atlas_sample_types]
        
        if (LOO_sample != "not_given" and LOO_sample in metadata.filename.unique()):
            df = df[df["Sample"] != LOO_sample]
        
        if df.shape[0] != 0:
            df["label"] = df["sample"].apply(lambda x: metadata[metadata["filename"] == x]["Label"].unique()[0])
            
            df["check_cigar"] = df["cigar"].apply(lambda x: bool(pattern.fullmatch(x)))
            df = df[df["check_cigar"] == True]
            
            df["Label"] = df["sample"].apply(lambda x: metadata[metadata["filename"] == x]["Label"].unique()[0])
                        
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
            import scipy
            tmpresdf = pd.DataFrame(data = [region], columns = ["region"])
            
            for group in atlas_labels[atlas_sample_types]:
                target_group = betadf[betadf["label"] == group].avg_beta.to_numpy()
                background_group = betadf[betadf["label"] != group].avg_beta.to_numpy()
                target_group = [item for item in target_group if np.isnan(item) == False]
                background_group = [item for item in background_group if np.isnan(item) == False]
                pval = scipy.stats.ttest_ind(target_group, 
                                     background_group).pvalue
                logFC = np.log2(0.0001 + (np.mean(target_group)/np.mean(background_group)))
                tmpresdf["{}_pval".format(group)] = pval
                tmpresdf["{}_logFC".format(group)] = logFC
            tmpresdf.to_csv(os.path.join(path_to_02_output, "ttest_res_tissue_and_wbc_samples_{}.csv".format(region)))
            outputstring = ",".join([str(item) for item in tmpresdf.loc[0].values])
            os.system("echo {} >> {}".format(outputstring, summ_test_file))

    else:
        print("region {} is already in the output dir {} of sample_type {}".format(region, path_to_02_output, atlas_sample_types))

def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', type=str,
        help='names of the region to collect reads from BAM file', action='store')
    parser.add_argument('--output', type=str, action = "store",
        help='Path to output file')
    parser.add_argument('--metadata', type=str, action = "store",
        help='Path to metadata file')
    parser.add_argument('--LOO_sample', type=str, action = "store",
        help='Choose a sample to be held out in LOOCV')
    parser.add_argument('--atlas_sample_types', type=str, action = "store",
        help='Choose the sample type, this can only be Tissue, WBC or cfDNA')
    return parser.parse_args(argv)

if __name__ == '__main__':
    main(parse_arguments(sys.argv[1:]))