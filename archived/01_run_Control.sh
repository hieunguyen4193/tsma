all_file_splits=$(ls ./datadir/CpG_clusters_whole_genome_radius_100_gt5_split_regions/*.split.regions);
bamdir="/mnt/DATASM14/hieunho/hieu_data/GW-bismark_cfDNA-highdepth_pair-end/04_bismark_deduplicated_filtered_Q30";
for region_file_split in $all_file_splits;do \
regions=$(cat ${region_file_split} | cut -d, -f7);
echo -e "working on region split number " $region_file_split "\n";
parallel -j 100 'python 01_collect_reads_from_BAM.py --region {} --output ./datadir/readdf_Control --metadata metadata.csv --sample_type Control --bamdir '$bamdir ::: $regions;
echo $region_file_split >> finished_split_regions_Control.txt;done