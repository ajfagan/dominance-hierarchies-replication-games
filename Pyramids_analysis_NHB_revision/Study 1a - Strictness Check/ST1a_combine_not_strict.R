# Study 1a - Combine human data with non-human animal data

rm(list = ls(all = TRUE))

library(tidyverse)
library(statnet)

# load CHILDREN data
load(file="dom_mat_SP_df_long_children_not_strict.Rdata")
load(file="dom_mat_SP_df_children_not_strict.Rdata")

# load ANIMAL data
load(file="dom_mat_SP_df_long_animals.Rdata")
load(file="dom_mat_SP_df_animals.Rdata")

# Combine the two datasets
dom_mat_SP_df_animals <-dom_mat_SP_df_animals%>%
  select(fileid,order,species,starts_with("X"),number_individuals,ttm,tpm)%>%
  rename(File_name=fileid)%>%
  rename(Group=order)%>%
  rename(Species=species)%>%
  rename(N = number_individuals)

dom_mat_SP_df_children<-dom_mat_SP_df_children%>%
  select(File_name,starts_with("X"),N,ttm,tpm)%>%
  mutate(Group="Children",Species="Homo sapiens")

combined_dom_mat_SP <- full_join(dom_mat_SP_df_children,dom_mat_SP_df_animals)
glimpse(combined_dom_mat_SP)
nrow(combined_dom_mat_SP)

save(combined_dom_mat_SP,file="combined_dom_mat_SP_not_strict.Rdata")
