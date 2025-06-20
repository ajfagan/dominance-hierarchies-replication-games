# The origins of pyramidal hierarchies #
# Study 1a animal data #
rm(list = ls(all = TRUE))

library(tidyverse)
library(statnet)

set.seed(123) # for reproducibility

# Get the data
# load data: a list with all the studies containing dataframe with contest matrix and metadata
load(file="../../data/dom.data.RData")
load(file="../../data/dom.metadata.RData")
dom.data[[1]] # for instance

# First add matrix if data in edgelist format
library(igraph) # pb compatibility with statnet

# function to add matrix if data are in edgelist format 
contest_matrix <-function(dom.data.study){ # For a study
  if (!is.null(dom.data.study$edgelist)){ # If there is an edgelist make matrix
    simple_edgelist<- data.frame(winner = dom.data.study$edgelist$winner, loser = dom.data.study$edgelist$loser)
    mygraph <- graph.data.frame(simple_edgelist)
    temp_mat <- get.adjacency(mygraph, sparse = FALSE)
    dom.data.study$matrix<-temp_mat # Add matrix to the study
    print(dom.data.study$metadata$species)
  }
  dom.data.study # Return the study
}
dom.data<-lapply(dom.data,contest_matrix) # Apply to all studies
detach("package:igraph", unload=TRUE) # Now detach package igraph to avoid conflicts

# Find dominance matrix from contests
dom_matrix<-function(contest_mat){
  N = ncol(contest_mat)
  temp_mat<-matrix(rep(0,N^2),N,N)
  for (i in 1:N){
    for (j in 1:N){
      if (contest_mat[i,j]>0 & contest_mat[i,j]>=contest_mat[j,i]){
        temp_mat[i,j]=1
      }
      if (i==j){temp_mat[i,j]=0}
    }
  }
  temp_mat
}
# For instance
dom_matrix(matrix(1:9,3,3))

# Find dominance matrix from contests under 
#  the more general definition of dominance 
#  determined by the odds of winning a contest.
dom_matrix_odds <-function(contest_mat, odds = 1){
  N = ncol(contest_mat)
  temp_mat<-matrix(rep(0,N^2),N,N)
  for (i in 1:N){
    for (j in 1:N){
      #prop_win = contest_mat[i,j]/(contest_mat[i,j] + contest_mat[j,i])
      if (contest_mat[i,j]>0 & ( contest_mat[i,j] >= contest_mat[j,i]/odds)){
        temp_mat[i,j]=1
      }
      if (i==j){temp_mat[i,j]=0}
    }
  }
  temp_mat
}

# Apply to dataset

dom_matrix(matrix(1:9,3,3))


# Apply to dataset
# dom_data_noNA<-lapply(dom.data, function(mat) {mat$matrix[is.na(mat$matrix)] <- 0;mat}) # replace NA by 0
# dom_mat<-lapply(dom_data_noNA,function(dat) {dom_matrix(dat$matrix)})
odds <- 2 # 1:1 dominance
dom_data_noNA<-lapply(dom.data, function(mat) {mat$matrix[is.na(mat$matrix)] <- 0;mat}) # replace NA by 0
dom_mat<-lapply(dom_data_noNA,function(dat) {dom_matrix_odds(dat$matrix, odds = odds)})

# Null model, randomise direction of links

# Randomise matrix
rand_dom_mat<-function(dom_mat){ 
  N = ncol(dom_mat)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (runif(1)>0.5){
        m1=dom_mat[i,j]
        m2=dom_mat[j,i]
        dom_mat[i,j]=m2
        dom_mat[j,i]=m1
        if ((dom_mat[i,j] == 1) && (dom_mat[j,i] == 1)){
          if (runif(1)>0.5) { # in 1/2 of cases the relation isn't symmetric
            if (runif(1)>0.5) { # in 1/4 of cases the relation becomes A dominates B
              dom_mat[i,j]=0
              dom_mat[j,i]=1}
            else { # in 1/4 of cases the relation becomes B dominates A
              dom_mat[i,j]=1
              dom_mat[j,i]=0}
          }
        }
      }
    }
  }
  dom_mat
}

###########################################################################################
## Previously consider triad_z function which generates z-scores given sampling under null
###########################################################################################

# Returns the z-score 
triad_z <- function(dom_mat,nrep=10){
  obs_val <- triad.census(dom_mat)
  sim_val<-replicate(nrep,triad.census(rand_dom_mat(dom_mat)),simplify =TRUE)
  mean_val<-apply(sim_val,1,mean)
  sd_val<-apply(sim_val,1,sd)
  z_score<-(obs_val-mean_val)/sd_val
  z_score
  # adding triadic indexes: triadic pyramidal metric (TPM) and triadic transitivity metric (TTM)
  obs_val_df<-data.frame(obs_val)
  tpm = obs_val_df$X021D/(obs_val_df$X021D+obs_val_df$X021U)
  ttm = obs_val_df$X030T/(obs_val_df$X030T+obs_val_df$X030C)
  
  return(cbind(z_score,ttm,tpm))
}
# For instance
triad_z(matrix(c(0,1,1,0,0,1,0,0,0),3,3))
# useful below
triad.census(matrix(c(0,1,1,0,0,1,0,0,0),3,3))

# Returns the z-score 
num_larger <- function(x) {
  
}

###########################################################################################
## New triad_z_exact function which generates p-values based on number of observed motifs
## under the null.
###########################################################################################

triad_z_exact <- function(dom_mat,nrep=10){
  obs_val <- triad.census(dom_mat)
  #print(obs_val)
  sim_val<-replicate(nrep,triad.census(rand_dom_mat(dom_mat)),simplify =TRUE)
  #print(dim(sim_val))

  pvalues <-  numeric(length = length(obs_val))
  signs <-  numeric(length = length(obs_val))
  for (i in 1:length(obs_val)) {
    sim_mvec <- sim_val[i,]
    obs_mval <- obs_val[i]
    num_larger <- sum(obs_mval > sim_mvec) + 0.5 * sum(obs_mval == sim_mvec)
    #print(num_larger)
    signs[i] <- sign(nrep/2 - num_larger)
    pvalues[i] <- 2 * min(c(nrep - num_larger, num_larger))/(nrep+1) + 1/(nrep + 1)
  }
  out_df <- data.frame(signs = signs, pvalues = pvalues)
  row.names(out_df) <- colnames(obs_val)
  return(out_df)
}
# For instance
triad_z_exact(matrix(c(0,1,1,0,0,1,0,0,0),3,3),nrep = 1000)
# useful below

motifs_names = colnames(triad_z(matrix(c(0,1,1,0,0,1,0,0,0),3,3)))


# Example of applying random generation
#sim_val<-replicate(100,triad.census(rand_dom_mat(dom_mat[[1]])),simplify =TRUE)

# Example of the frequency metric for pyramids:trees 
#  (4 corresponds to pyramid)
#sim_tpm <- sim_val[4,]/(sim_val[4,] + sim_val[5,])


###########################################################################################
## New triad_prop_exact function which generates p-values based on proportion of pyramids
## to trees under the null.
###########################################################################################

triad_prop_exact <- function(dom_mat,nrep=10){
  obs_val <- triad.census(dom_mat)
  #print(obs_val)
  sim_val<-replicate(nrep,triad.census(rand_dom_mat(dom_mat)),simplify =TRUE)
  #print(dim(sim_val))
  obs_val_df <- data.frame(obs_val)
  #obs_tpm = obs_val_df$X021D/(obs_val_df$X021D+obs_val_df$X021U)
  
  if (obs_val_df$X021D+obs_val_df$X021U > 0) {
    obs_tpm = obs_val_df$X021D/(obs_val_df$X021D+obs_val_df$X021U)
  } else {
    obs_tpm = 0
  }
  
  index_null <- which(sim_val[4,] + sim_val[5,] == 0)
  sim_tpm <- sim_val[4,]/(sim_val[4,] + sim_val[5,])
  sim_tpm[index_null] <- 0
  
  #ttm = obs_val_df$X030T/(obs_val_df$X030T+obs_val_df$X030C)
  
  #sim_tpm <- sim_val[4,]/(sim_val[4,] + sim_val[5,])
  
  num_larger <- sum(obs_tpm > sim_tpm) + 0.5 * sum(obs_tpm == sim_tpm)
  #print(num_larger)
  sign <- sign(nrep/2 - num_larger)
  pvalue <- 2 * min(c(nrep - num_larger, num_larger))/(nrep+1) + 1/(nrep + 1)


  out_df <- data.frame(signs = sign, tpm = obs_tpm, pvalues = pvalue)
  #row.names(out_df) <- colnames(obs_val)
  #print(out_df)
  return(out_df)
}

###########################################################################################
## Generating test statistics
###########################################################################################

# t1<-Sys.time()
dom_mat_z_score<-lapply(dom_mat,triad_z,nrep=1000)
# t2<-Sys.time()
# t2-t1

## Exact pvalues by sampling of triad counts by sampling network under the null of random dominance directions 
# t1<-Sys.time()
dom_mat_exact<-lapply(dom_mat,triad_z_exact,nrep=1000)
# t2<-Sys.time()
# t2-t1

## Exact pvalues by sampling of tpm proportion by sampling network under the null of random dominance directions 
# t1<-Sys.time()
dom_mat_prop_exact<-lapply(dom_mat,triad_prop_exact,nrep=1000)
# t2<-Sys.time()
# t2-t1

###########################################################################################
## Un-comment to save triad_z_exact data
###########################################################################################

orders <- sapply(dom.data, FUN = function(x){x$metadata$order})
#study_names <- names(dom_mat_exact)
study <- dom_mat_exact[1]
study_name <- names(study)
study_mat <- study[[1]]
n_motif <- dim(study_mat)[1]
temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep(orders[1], n_motif), motif = rownames(study_mat),  study_mat)
template <- temp_df
for (i in 2:length(dom_mat_exact)) {
  study <- dom_mat_exact[i]
  study_name <- names(study)
  study_mat <- study[[1]]
  n_motif <- dim(study_mat)[1]
  temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep(orders[i], n_motif), motif = rownames(study_mat),  study_mat)
  template <- rbind(template, temp_df)
}
rownames(template) <- NULL
write_csv(template, file = "ST1a_animals_exact_pvalues.csv")


###########################################################################################
## Un-comment to save triad_prop_exact data
###########################################################################################

#study_names <- names(dom_mat_exact)

orders <- sapply(dom.data, FUN = function(x){x$metadata$order})
study <- dom_mat_prop_exact[1]
study_name <- names(study)
study_mat <- study[[1]]
n_motif <- dim(study_mat)[1]
temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep(orders[1], n_motif),  study_mat)
template <- temp_df
for (i in 2:length(dom_mat_prop_exact)) {
  study <- dom_mat_prop_exact[i]
  study_name <- names(study)
  study_mat <- study[[1]]
  n_motif <- dim(study_mat)[1]
  temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep(orders[i], n_motif),  study_mat)
  template <- rbind(template, temp_df)
}
rownames(template) <- NULL
write_csv(template, file = "ST1a_animals_exact_prop_pvalues.csv")

###########################################################################################
## Previous significance pattern for triad_z data
###########################################################################################

rownames(dom_mat_z_score)=motifs_names
dom_mat_z_score[,1:6]

# Find significance profile (SP)
norm_factor<-apply(dom_mat_z_score[1:16,]^2,2,function(x) sqrt(sum(x,na.rm=T))) # calculate sqrt of SS
dom_mat_signif_prof<-t(t(dom_mat_z_score) * (1/norm_factor))
dom_mat_signif_prof<-dom_mat_signif_prof[c(4,5,6,9,10),] # only the 5 motifs we are interested in
dom_mat_signif_prof[,1:10]

# Plot SP
# Add metadata and create useful tibble
dom.metadata<-dom.metadata%>%
  select(-proportion_unknown,-dci,-ttri,-ds_steepness,-modified_landaus_h) # remove columns that are not used and do not exist for edgelist data
glimpse(dom.metadata)
n_occur <- data.frame(table(dom.metadata$groupid)) # added to report the number of networks with repeated measures
repeated_measures <- n_occur[n_occur$Freq > 1,]
length (repeated_measures$Var1) # repeated measures were reported for 48 groups, we keep only the first measure in these cases
dom_mat_signif_prof_df=data.frame(t(dom_mat_signif_prof))
dom_mat_signif_prof_df$fileid=rownames(dom_mat_signif_prof_df)
dom_mat_signif_prof_df$ttm = dom_mat_z_score["ttm",]
dom_mat_signif_prof_df$tpm = dom_mat_z_score["tpm",]
dom_mat_signif_prof_df=dom_mat_signif_prof_df[complete.cases(dom_mat_signif_prof_df),] # filter out complete structures
dom_mat_SP_df<-inner_join(dom.metadata,dom_mat_signif_prof_df,by="fileid")
dom_mat_SP_df=dom_mat_SP_df[match(unique(dom_mat_SP_df[,13]),dom_mat_SP_df[,13]),] # keep only the first measure for each network
dom_mat_SP_df_long = dom_mat_SP_df %>% # wide to long
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Dom_pattern",
    values_to = "norm_Z_Score")

# Reorder variables
# create a vector with the desired order
dom_pattern_order <- c("X021D", "X021U", "X021C","X030T","X030C")
dom_mat_SP_df_long=dom_mat_SP_df_long%>%
  mutate(order =  factor(order)) %>%
  mutate(Dom_pattern =  factor(Dom_pattern, levels = dom_pattern_order)) %>%
  arrange(order,species,Dom_pattern)

# Now plot SP
dom_mat_SP_plot<-ggplot(dom_mat_SP_df_long,aes(x = Dom_pattern,y=norm_Z_Score,colour=captivity,group=fileid))+
  geom_line(size=1)+
  labs(x="Dominance pattern",y="Normalised Z-score")+
  theme_bw()+
  ylim(-1,1)+
  facet_wrap(~order,ncol=2)+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12+2))
dom_mat_SP_plot
ggsave(dom_mat_SP_plot, filename="dom_mat_SP_plot.jpeg",width = 15, height=30, units = "cm")

# Save outcomes
dom_mat_SP_df_animals<-dom_mat_SP_df
dom_mat_SP_df_long_animals<-dom_mat_SP_df_long
save(dom_mat_SP_df_long_animals,file="dom_mat_SP_df_long_animals_odds2.Rdata")
save(dom_mat_SP_df_animals,file="dom_mat_SP_df_animals_odds2.Rdata")

