#!/bin/sh

# download data from https://osf.io/rg6nw
wget --content-disposition https://osf.io/rg6nw/download

# move the downloaded .zip file into the 'data' directory and unzip
mv Pyramids_Data_Codes_NHB.zip data
cd data
unzip Pyramids_Data_Codes_NHB.zip

# Copy the relevant files into the base data directory
mv Pyramids_analysis_NHB_revision/Studies_5ab_6ab_7ab_XPs_with_infants/data_pyramids_infants.csv .
mv Pyramids_analysis_NHB_revision/Studies_S1ab_2_3_4_XPs_with_adults/Data_pyramids_adults.csv .
mv Pyramids_analysis_NHB_revision/Study\ 1a/dom.* .
mv Pyramids_analysis_NHB_revision/Study\ 1a/Metadata.csv .
mv Pyramids_analysis_NHB_revision/Study\ 1a/list_children_dominance_matrix.Rdata .
mv Pyramids_analysis_NHB_revision/Study_1b/Primate_Tree.nex .

# Delete remaining files
rm -rf Pyramids_analysis_NHB_revision
rm Pyramids_Data_Codes_NHB.zip
