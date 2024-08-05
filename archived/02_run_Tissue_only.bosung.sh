regions=$(cat missing_regions_tissue_only.txt);
echo -e "working on region split number " $region_file_split "\n";
parallel -j 250 'python 02_conduct_ttest_one_class_versus_others.py --region {} --output ./outputdir --metadata metadata.csv --LOO_sample not_given --atlas_sample_types Tissue' ::: $regions;
