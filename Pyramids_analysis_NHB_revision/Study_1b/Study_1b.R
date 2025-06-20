# Phylogenetic analysis of dominance patterns

# The first part is an analysis of primate data
# The second part is the same analysis with the full tree of life

rm(list = ls(all = TRUE))

# Packages
library(rotl)
library(datelife)
library(tidyverse)
library(MCMCglmm)
library(ggpubr)
library(phytools)
library(phylobase)
library(phylosignal)

# Data
load(file="../Study 1a/combined_dom_mat_SP.Rdata")
glimpse(combined_dom_mat_SP)

####################################################!
### Study of Primates only  (with branch length) ###
####################################################!

# https://10ktrees.nunn-lab.org/

# Subset data with primates only
primates_dom_mat_SP<-combined_dom_mat_SP%>%
  filter(Group=="Primates" | Group=="Children")
length(unique(primates_dom_mat_SP$Species))
sort(unique(primates_dom_mat_SP$Species))

# Import 10k tree
primate_tree<-read.nexus("../../data/Primate_Tree.nex") # This is the corresponding tree

# change tip labels
primate_tree$tip.label <- gsub("_"," ",primate_tree$tip.label)
sort(primate_tree$tip.label)
sort(unique(primates_dom_mat_SP$Species[!gsub("_"," ",primates_dom_mat_SP$Species) %in% sort(primate_tree$tip.label)]))

# Change names to make them match between data and tree
primate_tree$node.label <- NULL # remove node names
primates_dom_mat_SP=primates_dom_mat_SP%>%
  mutate(tree_names=gsub("_"," ",Species))%>%
  mutate(tree_names = replace(tree_names, tree_names == "Gorilla beringei beringei","Gorilla beringei"))%>%
  mutate(tree_names = replace(tree_names, tree_names == "Papio hamadryas hamadryas","Papio hamadryas"))%>%
  mutate(tree_names = replace(tree_names, tree_names == "Macaca fuscata fuscata","Macaca fuscata"))%>%
  mutate(tree_names = replace(tree_names, tree_names == "Eulemur fulvus","Eulemur fulvus fulvus"))%>%
  mutate(tree_names = replace(tree_names, tree_names == "Pan troglodytes","Pan troglodytes troglodytes"))%>% # We lump together Pan troglodytes troglodytes and Pan troglodytes
  mutate(tree_names = replace(tree_names, tree_names == "Gorilla gorilla","Gorilla gorilla gorilla")) # We lump together Gorilla gorilla gorilla and Gorilla gorilla 

primates_dom_mat_SP_tr <- primates_dom_mat_SP[primates_dom_mat_SP$tree_names %in% primate_tree$tip.label, ]
setdiff(primates_dom_mat_SP$tree_names,primates_dom_mat_SP_tr$tree_names) # loosing 1 species, no tree data on Cercocebus sanjei
setdiff(primate_tree$tip.label,primates_dom_mat_SP_tr$tree_names)
glimpse(primates_dom_mat_SP_tr)

# Plot tree
pdf(file = "primate_tree.pdf",width=10,height=10)
plot.phylo(primate_tree,cex = 0.8)
axisPhylo() # Add the time axis
mtext("Time (myrs)", side = 1, line = 2, at = max(get("last_plot.phylo",envir = .PlotPhyloEnv)$xx) * 0.5) # add the axis name
dev.off()

# Model comparison to evaluate importance of phylogeny and group size

# Make the tree ultrametric for analysis
is.ultrametric(primate_tree)
primate_tree<-phytools::force.ultrametric(primate_tree,method = "nnls") # transforms the tree to make it ultrametric
is.ultrametric(primate_tree)

# Prepare for MCMCglmm models
primates_dom_mat_SP_tr<-data.frame(primates_dom_mat_SP_tr)
inv.phylo<-inverseA(primate_tree,nodes="TIPS")
pr <- list( # Some non committed priors
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 0.002))
)

set.seed(123) # for reproducibility

# Full model with phylogeny and group size
primate_phylo_N_model <- MCMCglmm(tpm ~ N,
                                  random = ~tree_names,
                                  family="gaussian",
                                  ginverse=list(tree_names=inv.phylo$Ainv),
                                  data = primates_dom_mat_SP_tr,
                                  prior = pr,
                                  verbose = FALSE
)
summary(primate_phylo_N_model)

# Only phylogeny
primate_phylo_model <- MCMCglmm(tpm ~ 1,
                                  random = ~tree_names,
                                  family="gaussian",
                                  ginverse=list(tree_names=inv.phylo$Ainv),
                                  data = primates_dom_mat_SP_tr,
                                  prior = pr,
                                  verbose = FALSE
)
summary(primate_phylo_model)

# Only group size
primate_N_model <- MCMCglmm(tpm ~ N,
                            random = ~tree_names,
                            family="gaussian",
                            data = primates_dom_mat_SP_tr,
                            prior = pr,
                            verbose = FALSE
)
summary(primate_N_model)

# Neither group size nor phylogeny
primate_model <- MCMCglmm(tpm ~ 1,
                          random = ~tree_names,
                          family="gaussian",
                          data = primates_dom_mat_SP_tr,
                          prior = pr,
                          verbose = FALSE
)
summary(primate_model)

# Model comparison based on DIC
summary(primate_phylo_N_model)$DIC
summary(primate_phylo_model)$DIC
summary(primate_N_model)$DIC
summary(primate_model)$DIC
summary(primate_phylo_N_model)$DIC-summary(primate_model)$DIC
summary(primate_phylo_model)$DIC-summary(primate_model)$DIC
summary(primate_N_model)$DIC-summary(primate_model)$DIC
# Conclusion: no detectable effect of network size (N) or phylogeny

# This is confirmed by calculating Pagel's lambda
lambda <- primate_phylo_model$VCV[,'tree_names']/(primate_phylo_model$VCV[,'tree_names']+primate_phylo_model$VCV[,'units'])
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)
# Above 0 (no signal) but still close to it and far from 1 (max signal)

# Plot the mean effect on the phylogeny (this ignores multiple values per species)

# compute mean values for each species
avg_tpm_species<-primates_dom_mat_SP_tr%>%
  group_by(tree_names)%>%
  summarise(mean_tpm=mean(tpm),mean_ttm=mean(ttm),mean_X021D=mean(X021D),
            mean_X021U=mean(X021U),mean_X021C=mean(X021C),mean_X030T=mean(X030T),mean_X030C=mean(X030C))%>%
  arrange(match(tree_names,primate_tree$tip.label))

# Plot tree with average value
primate_tree_p4d <- phylo4d(primate_tree, tip.data = avg_tpm_species$mean_tpm) 
barplot.phylo4d(primate_tree_p4d,trait.labels=c("Average triadic pyramidal metric"),tree.type = "phylo", tree.ladderize = TRUE,center = FALSE , scale = FALSE, tip.cex= 1.2, trait.cex= 0.9)
pdf(file = "primate_tree_tpm.pdf",width=10,height=10)
barplot.phylo4d(primate_tree_p4d,trait.labels=c("Average triadic pyramidal metric"),tree.type = "phylo", tree.ladderize = TRUE,center = FALSE , scale = FALSE, tip.cex= 1.2, trait.cex = 0.9)
dev.off()

#################################################################!
### Tree of life phylogeny for all species (no branch length) ###
#################################################################!

# Match names
combined_dom_mat_SP[76,11] <- "Pan_troglodytes_schweinfurthii" # info from Foerster_2016
combined_dom_mat_SP[77,11] <- "Pan_troglodytes_schweinfurthii" # info from Foerster_2016
combined_dom_mat_SP[317,11] <- "Pan_troglodytes_verus" # info from Wittig_2003

combined_dom_mat_SP_fulltree<-combined_dom_mat_SP%>%
  mutate(Match_names=gsub("_", " ", tolower(Species)))%>%
  mutate(Match_names = replace(Match_names, Match_names == "diacamma sp","diacamma"))%>%
  mutate(Match_names = replace(Match_names, Match_names == "rangifer tarandus","rangifer tarandus caribou")) # info from Barette_1986
  
# Find the tree
dom_taxa <- tnrs_match_names(unique(combined_dom_mat_SP_fulltree$Match_names), context = "Animals")
head(dom_taxa) # 100% match
dom_tr <- tol_induced_subtree(ott_id(dom_taxa)[is_in_tree(ott_id(dom_taxa))])
plot(dom_tr, show.tip.label = FALSE)
dom_tr 

# Tip labels are OTT Ids, change for dataset names
dom_tr$tip.label <- strip_ott_ids(dom_tr$tip.label, remove_underscores = TRUE)
sort(dom_tr$tip.label) 

# Remove node names for analysis
dom_tr$node.label <- NULL

# Add the data from matched names to the dataframe
combined_dom_mat_SP_fulltree<-combined_dom_mat_SP_fulltree%>%
  left_join(select(dom_taxa,search_string,unique_name,ott_id), by = c('Match_names' = 'search_string'))

# 5 species have no tip label and will not be analysed
setdiff(combined_dom_mat_SP_fulltree$unique_name,dom_tr$tip.label) # Tip labels == unique_name
combined_dom_mat_SP_fulltree=combined_dom_mat_SP_fulltree%>%
  filter(unique_name%in%dom_tr$tip.label)

# Now perform model comparison
# Warning! MCMCglmm does not work with tibble...
combined_dom_mat_SP_fulltree<-data.frame(combined_dom_mat_SP_fulltree)

# Prepare for MCMCglmm 
pr <- list( # Some non committed priors
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 0.002))
)
inv.phylo<-inverseA(dom_tr,nodes="TIPS",scale=TRUE) # Warning because there are no branch length for this analysis

set.seed(123) # for reproducibility

# Full model with phylogeny and group size
fulltree_phylo_N_model <- MCMCglmm(tpm ~ N,
                                   random = ~unique_name,
                                   family="gaussian",
                                   ginverse=list(unique_name=inv.phylo$Ainv),
                                   data = combined_dom_mat_SP_fulltree,
                                   prior = pr,
                                   verbose = FALSE
)
summary(fulltree_phylo_N_model)

# Only phylogeny
fulltree_phylo_model <- MCMCglmm(tpm ~ 1,
                                   random = ~unique_name,
                                   family="gaussian",
                                   ginverse=list(unique_name=inv.phylo$Ainv),
                                   data = combined_dom_mat_SP_fulltree,
                                   prior = pr,
                                   verbose = FALSE
)
summary(fulltree_phylo_model)

# Only group size
fulltree_N_model <- MCMCglmm(tpm ~ N,
                                   random = ~unique_name,
                                   family="gaussian",
                                   data = combined_dom_mat_SP_fulltree,
                                   prior = pr,
                                   verbose = FALSE
)
summary(fulltree_N_model)

# Neither group size nor phylogeny
fulltree_model <- MCMCglmm(tpm ~ 1,
                                   random = ~unique_name,
                                   family="gaussian",
                                   data = combined_dom_mat_SP_fulltree,
                                   prior = pr,
                                   verbose = FALSE
)
summary(fulltree_model)

# Model comparison based on DIC
summary(fulltree_phylo_N_model)$DIC
summary(fulltree_phylo_model)$DIC
summary(fulltree_N_model)$DIC
summary(fulltree_model)$DIC

summary(fulltree_phylo_N_model)$DIC-summary(fulltree_model)$DIC
summary(fulltree_phylo_model)$DIC-summary(fulltree_model)$DIC
summary(fulltree_N_model)$DIC-summary(fulltree_model)$DIC
# Conclusion: no effect of network size (N), no effect of phylogeny


# This is confirmed by calculating Pagel's lambda
lambda <- fulltree_phylo_model$VCV[,'unique_name']/(fulltree_phylo_model$VCV[,'unique_name']+fulltree_phylo_model$VCV[,'units'])
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)

# Plot the mean effect on the phylogeny (this ignores multiple values per species)

# compute mean values for each species
avg_tpm_species<-combined_dom_mat_SP_fulltree%>%
  group_by(unique_name)%>%
  summarise(mean_tpm=mean(tpm))%>%
  arrange(match(unique_name,dom_tr$tip.label))

# Plot tree with average value
full_tree_p4d <- phylo4d(dom_tr, tip.data = avg_tpm_species$mean_tpm) 
barplot.phylo4d(full_tree_p4d,trait.labels=c("Average triadic pyramidal metric"),tree.type = "phylo", tree.ladderize = TRUE,center = FALSE , scale = FALSE, tip.cex= 1.5, trait.cex= 1.2)
pdf(file = "full_tree_tpm.pdf",width=10,height=20)
barplot.phylo4d(full_tree_p4d,trait.labels=c("Average triadic pyramidal metric"),tree.type = "phylo", tree.ladderize = TRUE,center = FALSE , scale = FALSE, tip.cex= 1.5, trait.cex= 1.2)
dev.off()

