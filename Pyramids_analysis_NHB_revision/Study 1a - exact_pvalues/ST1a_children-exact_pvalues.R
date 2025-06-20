
# The origins of pyramidal hierarchies  
# Study 1a - human data

rm(list = ls(all = TRUE))

library(tidyverse)
library(statnet)
library(parallel)
library(rstatix)
library(broom)

set.seed(123) # for reproducibility

# get metadata
children_dominance_matadata=read_csv(file="../../data/Metadata.csv")

#read_csv()
# load data
load(file="../../data/list_children_dominance_matrix.Rdata")

# Find dominance matrix from contests
dom_matrix<-function(contest_mat){
  N = ncol(contest_mat)
  temp_mat<-matrix(rep(0,N^2),N,N)
  for (i in 1:N){
    for (j in 1:N){
      if (contest_mat[i,j]>0 & contest_mat[i,j]>contest_mat[j,i]){ # Here if the link is non zero symetric it is transformed into 0-0
        temp_mat[i,j]=1
      }
      if (i==j){temp_mat[i,j]=0} # Just to be sure...
    }
  }
  temp_mat
}
# For instance
dom_matrix(matrix(1:9,3,3))
dom_matrix(matrix(rep(1,9),3,3))


# Find weighted dominance matrix from contests
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
odds = 1
children_dom_mat<-lapply(children_dom_mat,dom_matrix_odds, odds = odds)
children_dom_mat[[9]] # For instance

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
      }
    }
  }
  dom_mat
}
rand_dom_mat(matrix(1:9,3,3)) # For instance

###########################################################################################
## Previously consider triad_z function which generates z-scores given sampling under null
###########################################################################################
triad_z <- function(dom_mat,nrep=10){
  obs_val <- triad.census(dom_mat)
  sim_val<-replicate(nrep,triad.census(rand_dom_mat(dom_mat)),simplify =TRUE)

  # Normalisation
  mean_val<-apply(sim_val,1,mean)
  sd_val<-apply(sim_val,1,sd)
  z_score<-(obs_val-mean_val)/sd_val
  
  # adding triadic indexes: triadic pyramidal metric (TPM) and triadic transitivity metric (TTM)
  obs_val_df<-data.frame(obs_val)
  tpm = obs_val_df$X021D/(obs_val_df$X021D+obs_val_df$X021U)
  ttm = obs_val_df$X030T/(obs_val_df$X030T+obs_val_df$X030C)
  
  return(cbind(z_score,ttm,tpm))
  
}
# For instance
triad_z(matrix(c(0,1,1,0,0,1,0,0,0),3,3))

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

motifs_names = colnames(triad_z(matrix(c(0,1,1,0,0,1,0,0,0),3,3))) # useful below



sim_val<-replicate(100,triad.census(rand_dom_mat(children_dom_mat[[1]])),simplify =TRUE)
sim_tpm <- sim_val[4,]/(sim_val[4,] + sim_val[5,])

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
children_z_score<-sapply(children_dom_mat,triad_z,nrep=1000)
# t2<-Sys.time()
# t2-t1

# t1<-Sys.time()
children_mat_exact<-lapply(children_dom_mat,triad_z_exact,nrep=1000)
# t2<-Sys.time()
# t2-t1

# t1<-Sys.time()
children_mat_prop_exact<-lapply(children_dom_mat,triad_prop_exact,nrep=1000)
# t2<-Sys.time()
# t2-t1


###########################################################################################
## Un-comment to save triad_z_exact data
###########################################################################################

study <- children_mat_exact[1]
study_name <- names(study)
study_mat <- study[[1]]
n_motif <- dim(study_mat)[1]
temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep("Children", n_motif), motif = rownames(study_mat),  study_mat)
template <- temp_df
for (i in 2:length(children_mat_exact)) {
  study <- children_mat_exact[i]
  study_name <- names(study)
  study_mat <- study[[1]]
  n_motif <- dim(study_mat)[1]
  temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep("Children", n_motif), motif = rownames(study_mat),  study_mat)
  template <- rbind(template, temp_df)
}
rownames(template) <- NULL
write_csv(template, file = "ST1a_children_exact_pvalues.csv")


###########################################################################################
## Un-comment to save triad_prop_exact data
###########################################################################################

study <- children_mat_prop_exact[1]
study_name <- names(study)
study_mat <- study[[1]]
n_motif <- dim(study_mat)[1]
temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep("Children", n_motif),  study_mat)
template <- temp_df
for (i in 2:length(children_mat_prop_exact)) {
  study <- children_mat_prop_exact[i]
  study_name <- names(study)
  study_mat <- study[[1]]
  n_motif <- dim(study_mat)[1]
  temp_df <- cbind(study_name = rep(study_name, n_motif), order = rep("Children", n_motif),  study_mat)
  template <- rbind(template, temp_df)
}
rownames(template) <- NULL
write_csv(template, file = "ST1a_children_exact_prop_pvalues.csv")


###########################################################################################
## Previous significance pattern for triad_z data
###########################################################################################

rownames(children_z_score)=motifs_names
children_z_score[1:nrow(children_z_score),1:6]

norm_factor<-apply(children_z_score[1:16,]^2,2,function(x) sqrt(sum(x,na.rm=T))) # calculate sqrt of SS
children_signif_prof<-t(t(children_z_score) * (1/norm_factor))
children_signif_prof<-children_signif_prof[c(4,5,6,9,10),] # finally we select the 5 motifs we are analysing
children_signif_prof[,1:6]

# Plot SP
# Add metadata and create useful tibble
glimpse(children_dominance_matadata)
children_signif_prof_df=data.frame(t(children_signif_prof))
children_signif_prof_df$ttm = children_z_score["ttm",]
children_signif_prof_df$tpm = children_z_score["tpm",]
children_signif_prof_df$File_name=rownames(children_signif_prof_df)
children_SP_df<-inner_join(children_dominance_matadata[,1:4],children_signif_prof_df,by="File_name")
children_SP_df_long = children_SP_df %>% # wide to long
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Dom_pattern",
    values_to = "norm_Z_Score")

# Reorder variables
dom_pattern_order <- c("X021D", "X021U", "X021C","X030T","X030C") # create a vector with the desired order of motifs
children_SP_df_long=children_SP_df_long%>%
  mutate(Dom_pattern =  factor(Dom_pattern, levels = dom_pattern_order))

# Now plot SP
children_SP_plot<-ggplot(children_SP_df_long,aes(x = Dom_pattern,y=norm_Z_Score,group=File_name))+
  geom_line(size=1)+
  labs(x="Dominance pattern",y="Normalised Z-score")+
  ylim(-1,1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12+2))
children_SP_plot
ggsave(children_SP_plot, filename="children_SP_plot.jpeg",width = 15, height=30, units = "cm")

# Save outcomes for comparison with animals
dom_mat_SP_df_children<-children_SP_df
dom_mat_SP_df_long_children<-children_SP_df_long
save(dom_mat_SP_df_long_children,file="dom_mat_SP_df_long_children_odds2.Rdata")
save(dom_mat_SP_df_children,file="dom_mat_SP_df_children_odds2.Rdata")


# Complementary analysis: Effect of Age
# NB: This analysis has been added at the request of a reviewer.

# Prepare data
children_SP_df<-inner_join(children_dominance_matadata[,1:6],children_signif_prof_df,by="File_name") # get age values
children_SP_df$Mean_age <- scale(children_SP_df$Mean_age, scale = FALSE) # mean center age
children_SP_df$Mean_age_above_Mdn <- as.factor(children_SP_df$Mean_age_above_Mdn)
str(children_SP_df)

# Analysis 1: effect of Mean-Age as a continuous var (mean-centered, in month)
Age.tpm.lm <- lm(tpm ~ Mean_age, data = children_SP_df)
summary(Age.tpm.lm) # no effect of Age
tidy(Age.tpm.lm, conf.int = TRUE)

# same analysis but with standardized variables (for standardized coefficients)
Age.tpm.lm2 <- lm(scale(tpm) ~ scale(Mean_age), data = children_SP_df)
summary(Age.tpm.lm2) # no effect of Age
tidy(Age.tpm.lm2, conf.int = TRUE)


# Analysis 2: effect of Mean-Age as a categorical var
# descriptive statistics
tapply(children_SP_df$tpm, children_SP_df$Mean_age_above_Mdn, summary) # Median and mean values are close across age groups

# analysiq
Age_catego.tpm.lm <- lm(tpm ~ Mean_age_above_Mdn, data = children_SP_df)
summary(Age_catego.tpm.lm) # no effect of Age
tidy(Age_catego.tpm.lm, conf.int = TRUE)

# same analysis but with standardized variables (for standardized coefficients)
Age_catego.tpm.lm2 <- lm(scale(tpm) ~ Mean_age_above_Mdn, data = children_SP_df)
summary(Age_catego.tpm.lm2) # no effect of Age
tidy(Age_catego.tpm.lm2, conf.int = TRUE)

# Analysis 3:separate analyses per age group (younger vs. older)
# compare the triadic pyramidal metric to chance in each sub-group
tpm.per.age.group <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_test(tpm ~ 1, paired = FALSE, mu = 0.5, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))

tpm.per.age.group # mean tpms are above chance in both age groups

# computes effect sizes, takes time because of bootstrapping
effect.size.tpm <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_effsize(tpm ~ 1, paired = FALSE, mu = 0.5, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size.tpm

# Analyses on z-scores
# Pairwise comparisons between z-scores for triadic pyramids and triadic trees

# Convert to long format
children_SP_df_long <- gather(children_SP_df, pattern, norm_z_score, c("X021D","X021U","X021C","X030T", "X030C"), factor_key=TRUE)

# filter data
SP_long_for_compa <-dplyr::filter(children_SP_df_long, pattern == "X021D" | pattern == "X021U")

# run the analysis
# compute p-values & corrected p-values
Ztree.per.age.group <- SP_long_for_compa %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_test(norm_z_score ~ pattern, paired = TRUE, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
Ztree.per.age.group

# computes effect sizes
# Warning: takes time because of bootstrapping
effect.size.pairwise <- SP_long_for_compa %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_effsize(norm_z_score ~ pattern, paired = TRUE, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size.pairwise

# compare Z-scores to chance for triadic pyramids in each sub-group
# descriptive statistics for pyramids
tapply(children_SP_df$X021D, children_SP_df$Mean_age_above_Mdn, summary) # Median and mean values

# analysis
Zpyramid.per.age.group <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_test(X021D ~ 1, paired = FALSE, mu = 0, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))

Zpyramid.per.age.group # mean z-scores for triadic pyramids are above chance in both age groups

# computes effect sizes, takes time because of bootstrapping
effect.size.Zpyramid <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_effsize(X021D ~ 1, paired = FALSE, mu = 0, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size.Zpyramid

# descriptive statistics for trees
tapply(children_SP_df$X021U, children_SP_df$Mean_age_above_Mdn, summary) # Median and mean values

# compare Z-scores for triadic trees to chance in each sub-group
Ztree.per.age.group <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_test(X021U ~ 1, paired = FALSE, mu = 0, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))

Ztree.per.age.group # mean z-scores for triadic trees do not differ from chance in any age group

# computes effect sizes, takes time because of bootstrapping
effect.size.Ztree <- children_SP_df %>%
  group_by(Mean_age_above_Mdn) %>%
  wilcox_effsize(X021U ~ 1, paired = FALSE, mu = 0, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size.Ztree

