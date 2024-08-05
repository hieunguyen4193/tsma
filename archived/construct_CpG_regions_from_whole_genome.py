import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import pathlib
import pysam
import pyfaidx

path_to_all_fa = "/datassd/hieunguyen/ECD/storage/resources/hg19"

all_fa_files = [file for file in pathlib.Path(path_to_all_fa).glob("*.fa")]

import numpy as np
radius = 100
output = open("CpG_clusters_whole_genome_radius_{}.bed".format(radius), 'w')
    
for file in sorted(all_fa_files):
    filename = file.name.strip(".fa")
    
    print("Working on {}".format(file))
    fasta_file = pyfaidx.Fasta(file)

    ref_seq = fasta_file.get_seq(name = file.name.strip(".fa"), start = 1, end = fasta_file[file.name.strip(".fa")][-1].start)
    ref_seq = ref_seq.seq

    def find_CpG(ref_seq):
        import re
        all_C_positions = [m.start(0) for m in re.finditer("CG", ref_seq)]
        all_c_positions = [m.start(0) for m in re.finditer("cg", ref_seq)]
        all_cG_positions = [m.start(0) for m in re.finditer("cG", ref_seq)]
        return all_C_positions + all_c_positions + all_cG_positions

    all_cpg_positions = sorted(find_CpG(ref_seq))

    all_pos = all_cpg_positions.copy()

    anchor = all_pos[0]

    tmp_cpg_clusters = dict()
    cpg_cluster_idx = 0
    tmp_cpg_clusters[cpg_cluster_idx] = [anchor]

    for item in tqdm(all_pos[1:]):
        if (item - anchor) <= radius:
            tmp_cpg_clusters[cpg_cluster_idx].append(item)
        else:
            cpg_cluster_idx += 1
            anchor = item        
            tmp_cpg_clusters[cpg_cluster_idx] = [anchor]

            all_cpg_cluster = dict()

    for key in tmp_cpg_clusters.keys():
        if (len(tmp_cpg_clusters[key]) >= 4):
            all_cpg_cluster[key] = tmp_cpg_clusters[key]
            
    for key in all_cpg_cluster.keys():
        cluster = all_cpg_cluster[key]
        start = np.min(cluster) - 1
        end = np.max(cluster) + 1
        chrom = filename.strip("chr")
        num_cpg = len(cluster)
        region_name = "{}:{}-{}".format(chrom, start, end)

        output.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + "{}_bin_{}".format(chrom, key) + "\t" + str(num_cpg) + "\t" + region_name + '\n')
