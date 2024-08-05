all_file_splits=$(ls ./datadir/CpG_clusters_whole_genome_radius_100_gt5_split_regions/*.split.regions);
for region_file_split in $all_file_splits;do \
regions=$(cat ${region_file_split} | cut -d, -f7);
echo -e "working on region split number " $region_file_split "\n";
parallel -j 250 'python 02_conduct_ttest_one_class_versus_others.py --region {} --output ./outputdir --metadata metadata.csv --LOO_sample not_given --atlas_sample_types Tissue' ::: $regions;done;
