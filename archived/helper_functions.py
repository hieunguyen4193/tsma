import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pathlib
from tqdm import tqdm
import re
pattern = re.compile("^[1-9][0-9]M")
import warnings
warnings.filterwarnings("ignore")
import pyfaidx

#####--------------------------------------------------------------------------#####
##### GENERAL CONFIGURATIONS
#####--------------------------------------------------------------------------#####

readdf_dir = dict()
readdf_dir["Tissue"] = "./datadir/readdf_tumor"
readdf_dir["WBC"] = "./datadir/readdf_WBC"
readdf_dir["cfDNA"] = "./datadir/readdf_cfDNA"
readdf_dir["Control"] = "./datadir/readdf_Control"

atlas_labels = dict()
atlas_labels["Tissue"] = ["Liver", "Breast", "Gastric", "Lung", "CRC"]
atlas_labels["cfDNA"] = ["Liver", "Breast", "Gastric", "Lung", "CRC", "Control"]
atlas_labels["Tissue,WBC"] = ["Liver", "Breast", "Gastric", "Lung", "CRC", "WBC"]
atlas_labels["Tissue,Control"] = ["Liver", "Breast", "Gastric", "Lung", "CRC", "Control"]
atlas_labels["Tissue,WBC,Control"] = ["Liver", "Breast", "Gastric", "Lung", "CRC", "Control", "WBC"]

path_to_all_fa = "/datassd/hieunguyen/ECD/storage/resources/hg19"

summ_file_columns = dict()
summ_file_columns["Tissue"] = "region,Liver_pval,Liver_logFC,Breast_pval,Breast_logFC,Gastric_pval,Gastric_logFC,Lung_pval,Lung_logFC,CRC_pval,CRC_logFC"
summ_file_columns["cfDNA"] = "region,Liver_pval,Liver_logFC,Breast_pval,Breast_logFC,Gastric_pval,Gastric_logFC,Lung_pval,Lung_logFC,CRC_pval,CRC_logFC,Control_pval,Control_logFC"
summ_file_columns["Tissue,WBC"] = "region,Liver_pval,Liver_logFC,Breast_pval,Breast_logFC,Gastric_pval,Gastric_logFC,Lung_pval,Lung_logFC,CRC_pval,CRC_logFC,WBC_pval,WBC_logFC"
summ_file_columns["Tissue,Control"] = "region,Liver_pval,Liver_logFC,Breast_pval,Breast_logFC,Gastric_pval,Gastric_logFC,Lung_pval,Lung_logFC,CRC_pval,CRC_logFC,Control_pval,Control_logFC"
summ_file_columns["Tissue,WBC,Control"] = "region,Liver_pval,Liver_logFC,Breast_pval,Breast_logFC,Gastric_pval,Gastric_logFC,Lung_pval,Lung_logFC,CRC_pval,CRC_logFC,Control_pval,Control_logFC,WBC_pval,WBC_logFC"

def get_refseq(path_to_all_fa, chrom, start, end):
    
    refseq = pyfaidx.Fasta(os.path.join(path_to_all_fa, "chr{}.fa".format(chrom)))
    return(str.upper(refseq.get_seq(name = "chr{}".format(chrom), start = start, end = end).seq))
            
def get_CpG_status(read_start, read_end, read, cpg_pos, mode = "string"):    
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

def deconvo(inputdf, input_atlas, atlas_sample_types):
    resdf = pd.DataFrame(data = atlas_labels[atlas_sample_types], columns = ["TOO"])
    for samplename in inputdf["sample"].unique():
        sampledf = inputdf.set_index("sample").loc[samplename].reset_index()
        sampledf.columns = ["region"] + [samplename]
        
        deconv_atlas = input_atlas.T.reset_index()
        deconv_atlas.columns = ["region"] + list(deconv_atlas.columns[1:])
        deconv_atlas = deconv_atlas.merge(sampledf, right_on = "region", left_on = "region")
        deconv_atlas = deconv_atlas[~deconv_atlas[samplename].isna()]
        if deconv_atlas.shape[0] != 0:
            from scipy import optimize
            mixture, residual = optimize.nnls(deconv_atlas[atlas_labels[atlas_sample_types]].to_numpy(), deconv_atlas[samplename].to_numpy())
            mixture /= np.sum(mixture)
            tmpresdf = pd.DataFrame(data = atlas_labels[atlas_sample_types], columns = ["TOO"])
            tmpresdf[samplename] = mixture
        else:
            tmpresdf = pd.DataFrame(data = atlas_labels[atlas_sample_types], columns = ["TOO"])
            tmpresdf[samplename] = ["NaN deconv." for item in range(6)]
        resdf = resdf.merge(tmpresdf, right_on = "TOO", left_on = "TOO")
    return resdf