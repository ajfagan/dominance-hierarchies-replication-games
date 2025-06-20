# Analysis of adults' data, pyramids paper (Studies S1ab, 2-4)
rm(list = ls(all = TRUE))

## load required packages
ipak <- function (pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Install devtools package if necessary
if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

packages <- c("tidyverse","devtools", "ggthemes", "cowplot", "rcompanion", "emmeans", "afex", "bruceR","rstatix","pairwiseCI","ggpubr","effectsize", "viridis")

ipak(packages)

set.seed(123) # for reproducibility

## Load data

pyramid_adults.long <- read.csv("../../data/Data_pyramids_adults.csv", header = TRUE, sep =";")
pyramid_adults.long$study <- as.factor(pyramid_adults.long$study)
pyramid_adults.long$gender <- as.factor(pyramid_adults.long$gender)
pyramid_adults.long$names <- as.factor(pyramid_adults.long$names)
pyramid_adults.long$test_type <- as.factor(pyramid_adults.long$test_type)
pyramid_adults.long$Inference <- as.factor(pyramid_adults.long$Inference)
pyramid_adults.long$age <- as.numeric(pyramid_adults.long$age)

str(pyramid_adults.long)

## demographics
# age
tapply(pyramid_adults.long$age, pyramid_adults.long$xp, summary)

#gender
# NB: repeated measures in long format for ST2, values have to be divided by 2 for this study
tapply(pyramid_adults.long$gender, pyramid_adults.long$xp, summary)

# set the effect size for ANOVAs to partial eta squared
afex_options(es_aov = "pes")

# function performing all the analyses for studies S1ab, 3-4
analysis.adults <- function(the_study) {

# ANOVA followed up by t-tests
# ANOVA
exp1 <-dplyr::filter(pyramid_adults.long, study == the_study)
summary(exp1)

a1 <- aov_car(scoreInf ~ Inference*names + Error(id), exp1)
print("anova")
print(a1)


#Follow-up tests

#Comparison of inference scores to chance in each Study
stat.test <- exp1 %>% 
  group_by(Inference) %>%
  t_test(scoreInf ~ 1, paired = FALSE, mu = 0.5, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))

print("t-tests")
print(stat.test)

# effect sizes
effect.size <- exp1 %>% 
  group_by(Inference) %>%
  cohens_d(scoreInf ~ 1, paired = FALSE, mu = 0.5, ci = TRUE, ci.type = "perc", nboot = 1000)

print("effect sizes")
print(effect.size)

# create table summarizing pairwise comparisons
 stats.t_test.chance <- cbind(stat.test,effect.size[,4], effect.size[,7:9])
 
# write csv
 stats.t_test.chance <- apply(stats.t_test.chance,2,as.character)
 
 write.csv(stats.t_test.chance, file = paste(the_study, ".one_sample_t_test.csv", sep=""))
}

# Study S1
analysis.adults("STS1")
F_to_eta2(f = c(31.29,.41,.25), df = c(1,1,1), df_error = c(36,36,36))# added to get 95%CI for partial eta-squared

#Study 2

# ANOVA ST2
exp1 <-dplyr::filter(pyramid_adults.long, study == "ST2")
summary(exp1)

a1 <- aov_car(scoreInf ~ names*test_type*lingformat + Error(id/(test_type)), exp1)
print("anova for ST2")
print(a1)
F_to_eta2(f = c(6.19,2.43,.37,4.45,6.80,3.84,0.06), df = c(1,1,1,1,1,1,1), df_error = c(36,36,36,36,36,36,36))# added to get 95%CI for partial eta-squared


# Follow-up t-tests ST2
stat.test <- exp1 %>% 
  group_by(Inference) %>%
  t_test(scoreInf ~ 1, paired = FALSE, mu = 0.5, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
stat.test

effect.size <- exp1 %>% 
  group_by(Inference) %>%
  cohens_d(scoreInf ~ 1, paired = FALSE, mu = 0.5, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size

# create table summarizing pairwise comparisons
stats.t_test.chance <- cbind(stat.test,effect.size[,4], effect.size[,7:9])

# write csv
stats.t_test.chance <- apply(stats.t_test.chance,2,as.character)
write.csv(stats.t_test.chance, file = paste("ST2.one_sample_t_test.csv", sep=""))

# Study 3 
analysis.adults("ST3")
F_to_eta2(f = c(13.99,.58,1.51), df = c(1,1,1), df_error = c(36,36,36))# added to get 95%CI for partial eta-squared


# Study 4
analysis.adults("ST4")
F_to_eta2(f = c(14.51,.63,.61), df = c(1,1,1), df_error = c(36,36,36))# added to get 95%CI for partial eta-squared


#### Plot S1ab
pyramid_adults.graph <- read.csv("../../data/Data_pyramids_adults.csv", header = TRUE, sep =";")
pyramid_adults.graph <-dplyr::filter(pyramid_adults.graph, study == "STS1")

pyramid_adults.graph$Inference <- factor(pyramid_adults.graph$Inference, levels=c("test", "ctrl"))

p <- ggplot(pyramid_adults.graph,aes(x=Inference, y=scoreInf, fill =Inference,palette="npg")) +
  geom_boxplot(width=0.25, outlier.alpha = 0.5, outlier.size = 0.7, color="black", lwd = 0.5) +
  geom_jitter(height = 0.03, width=0, alpha=0.5, size = 0.7) +
  stat_summary(geom = "point",
               fun.y = "mean",
               col = "red",
               size = 5,
               shape = 20,
               fill = "red") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, size = 1, colour = "red") +
  facet_grid(cols=vars(study), scales = "free_x")+
  ylab("Inference Score")+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size=12, face="bold"))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, option = "D")


p +theme(legend.position = "None")
ggsave("Plot_adults_S1.pdf", device = "pdf", width = 8, height = 10, units = c("cm"), dpi = 900)

# Plot Studies 2-4
pyramid_adults.graph <- read.csv("../../data/Data_pyramids_adults.csv", header = TRUE, sep =";")
pyramid_adults.graph <-dplyr::filter(pyramid_adults.graph, study != "STS1")

p <- ggplot(pyramid_adults.graph,aes(x=Inference, y=scoreInf, fill =Inference,palette="npg")) +
  geom_boxplot(width=0.25, outlier.alpha = 0.5, outlier.size = 0.7, color="black", lwd = 0.5) +
  geom_jitter(height = 0.03, width=0, alpha=0.5, size = 0.7) +
  stat_summary(geom = "point",
               fun = "mean",
               col = "red",
               size = 5,
               shape = 20,
               fill = "red") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.1, size = 1, colour = "red") +
  facet_grid(cols=vars(study), scales = "free_x")+
  ylab("Inference Score")+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size=12, face="bold"))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, option = "D")

p +theme(legend.position = "None")
ggsave("Plot_adults_S2-4.pdf", device = "pdf", width = 18, height = 10, units = c("cm"), dpi = 900)

