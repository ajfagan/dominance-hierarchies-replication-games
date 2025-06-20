#!/bin/sh

# Install libraries
Rscript install_relevant_libraries_I4R.R

# Run each Rscript
cd Pyramids_analysis_NHB_revision/Study\ 1a
Rscript ST1a_animals.R
Rscript ST1a_children.R
Rscript ST1a_combine.R
Rscript ST1a_plots.R

cd ../Study\ 1a\ -\ exact_pvalues
Rscript ST1a_animals-exact_pvalues.R
Rscript ST1a_children-exact_pvalues.R
Rscript combined_exact_pvalues.R

cd ../Study\ 1a\ -\ reordered\ removal
Rscript ST1a_animals.R
Rscript ST1a_children.R
Rscript ST1a_combine.R
Rscript ST1a_plots.R


cd ../Study\ 1a\ -\ Strictness\ Check
Rscript ST1a_animals.R
Rscript ST1a_children.R
Rscript ST1a_combine.R
Rscript ST1a_plots.R
Rscript ST1a_animals_strict.R
Rscript ST1a_children_strict.R
Rscript ST1a_combine_strict.R
Rscript ST1a_plots_strict.R
Rscript ST1a_animals_not_strict.R
Rscript ST1a_children_not_strict.R
Rscript ST1a_combine_not_strict.R
Rscript ST1a_plots_not_strict.R


cd ../Study\ 1a\ -\ Norm\ Fix
Rscript ST1a_animals.R
Rscript ST1a_children.R
Rscript ST1a_combine.R
Rscript ST1a_plots.R
Rscript ST1a_animals_all.R
Rscript ST1a_children_all.R
Rscript ST1a_combine_all.R
Rscript ST1a_plots_all.R

cd ../Study_1b
Rscript Study_1b.R

cd ../Studies_S1ab_2_3_4_XPs_with_adults
Rscript Studies_S1ab_2_3_4.R

cd ../Studies_5ab_6ab_7ab_XPs_with_infants
Rscript Studies_5ab_6ab_7ab.R


