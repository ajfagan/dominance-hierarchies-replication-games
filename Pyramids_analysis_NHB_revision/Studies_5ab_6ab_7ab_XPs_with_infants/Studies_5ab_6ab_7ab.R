# Analysis of infants' data, pyramids paper (Studies 5ab, 6ab, 7ab)
rm(list = ls(all = TRUE))

# load required packages
ipak <- function (pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Install devtools package if necessary
if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")


packages <- c("tidyverse","devtools", "ggthemes", "cowplot", "rcompanion", "emmeans", "afex", "bruceR","rstatix","pairwiseCI","ggpubr", "viridis")

ipak(packages)

set.seed(123) # for reproducibility

## Load data

pyramid_infants.wide <- read.csv("../../data/data_pyramids_infants.csv", header = TRUE, sep =";")

# convert to long format
pyramid_infants.long <- gather(pyramid_infants.wide, test_type, LogLT, logcoh:logincoh) 
pyramid_infants.long.raw <- gather(pyramid_infants.wide, test_type, rawLT, coh:incoh)

pyramid_infants.long$first_test <- as.factor(pyramid_infants.long$first_test)
pyramid_infants.long$test_type <- as.factor(pyramid_infants.long$test_type)
pyramid_infants.long$group <- as.factor(pyramid_infants.long$group)
pyramid_infants.long$study <- as.factor(pyramid_infants.long$study)
pyramid_infants.long$gender <- as.factor(pyramid_infants.long$gender)

pyramid_infants.long.raw$first_test <- as.factor(pyramid_infants.long.raw$first_test)
pyramid_infants.long.raw$test_type <- as.factor(pyramid_infants.long.raw$test_type)
pyramid_infants.long.raw$group <- as.factor(pyramid_infants.long.raw$group)
pyramid_infants.long.raw$study <- as.factor(pyramid_infants.long.raw$study)

str(pyramid_infants.long)
str(pyramid_infants.long.raw)

# demographics
# age
tapply(pyramid_infants.long$age, pyramid_infants.long$xp, summary)

# gender (repeated measures in long format; values need to be divided by 2)
tapply(pyramid_infants.long$gender, pyramid_infants.long$xp, summary)


# set the effect size for ANOVAs to partial eta squared
afex_options(es_aov = "pes")


# function performing all the analyses for each infant study
analysis.infants <- function(the_study) {
# ANOVA followed up by t-tests

# ANOVA
exp1 <-dplyr::filter(pyramid_infants.long, study == the_study)
summary(exp1)

a1 <- aov_car(LogLT ~ group*test_type + Error(id/(test_type)), exp1)
print("anova")
print(a1)

# follow-up on LogLT
stat.test <- exp1 %>% 
  group_by(group) %>%
  t_test(LogLT ~ test_type, paired = TRUE, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))

print("t-tests")
print(stat.test)

# effect sizes
effect.size <- exp1 %>% 
  group_by(group) %>%
  cohens_d(LogLT ~ test_type, paired = TRUE,ci = TRUE, ci.type = "perc", nboot = 1000)

print("effect sizes")
print(effect.size)

# create table summarizing pairwise comparisons
stats.t_test.chance <- cbind(stat.test,effect.size[,4], effect.size[,7:9])

# write csv
stats.t_test.chance <- apply(stats.t_test.chance,2,as.character)
write.csv(stats.t_test.chance, file = paste(the_study, ".t_tests.csv", sep=""))

# follow-up on raw data
# Wilcoxon tests on raw data

# test condition
# prepare raw data
exp1.experimental.raw <-dplyr::filter(pyramid_infants.long.raw, study == the_study & group == "test")

# wilcoxon test for matched pairs
#wilcoxon.experimental <- wilcox.test(rawLT ~ test_type, data=exp1.experimental.raw, paired=TRUE, exact=FALSE)
wilcoxon.experimental <- wilcox.test(
  exp1.experimental.raw$rawLT[1:20],  
  exp1.experimental.raw$rawLT[21:40], 
  paired=TRUE, exact=FALSE
)

# compute effect size r & 95% CI using the percentile bootstrapping method
ES.wilcoxon.experimental <- wilcoxonPairedR(x = exp1.experimental.raw$rawLT,g = exp1.experimental.raw$test_type,
                ci = TRUE,conf = 0.95, R = 10000)

print("Wilcoxon test for experimental condition (raw data)")
print(wilcoxon.experimental)
print("Wilcoxon effect size")
print(ES.wilcoxon.experimental)

# control condition
# prepare raw data
exp1.control.raw <-dplyr::filter(pyramid_infants.long.raw, study == the_study & group == "ctrl")

# wilcoxon test for matched pairs
#wilcoxon.control <- wilcox.test(rawLT ~ test_type, data=exp1.control.raw, paired=TRUE, exact=FALSE)
wilcoxon.control <- wilcox.test(
  exp1.control.raw$rawLT[1:20],
  exp1.control.raw$rawLT[21:40], 
  paired=TRUE, exact=FALSE
)

# compute effect size r & 95% CI using the percentile bootstrapping method
ES.wilcoxon.control <-wilcoxonPairedR(x = exp1.control.raw$rawLT,g = exp1.control.raw$test_type,
                ci = TRUE,conf = 0.95, R = 10000)

print("Wilcoxon test for control condition (raw data)")
print(wilcoxon.control)
print("Wilcoxon effect size")
print(ES.wilcoxon.control)
}


# Study 5ab
analysis.infants("ST5")
F_to_eta2(f = c(.33,4.80,5.93), df = c(1,1,1), df_error = c(38,38,38))# added to get 95%CI for partial eta-squared

# Study 6ab
analysis.infants("ST6")
F_to_eta2(f = c(1.32,.52,6.03), df = c(1,1,1), df_error = c(38,38,38))# added to get 95%CI for partial eta-squared

# Study 7ab
analysis.infants("ST7")
F_to_eta2(f = c(.33,0,.03), df = c(1,1,1), df_error = c(38,38,38))# added to get 95%CI for partial eta-squared


# Plots

# rename levels and reorder factors to facilitate plotting 
levels(pyramid_infants.long.raw$test_type)[levels(pyramid_infants.long.raw$test_type)=="coh"] <- "coherent"
levels(pyramid_infants.long.raw$test_type)[levels(pyramid_infants.long.raw$test_type)=="incoh"] <- "incoherent"

pyramid_infants.long.raw$study = factor(pyramid_infants.long.raw$study, levels=c('ST6','ST7','ST5'))
pyramid_infants.long.raw$group = factor(pyramid_infants.long.raw$group, levels=c('test','ctrl'))

# graph for paired data

# Plots
# ST5ab
data.plot.5ab <-dplyr::filter(pyramid_infants.long.raw, study == "ST5")
p <- ggpaired(data.plot.5ab, 
              x="test_type", 
              y="rawLT", 
              width = 0.3,
              ylim = c(0, 60),
              color="black",
              line.color="gray",
              line.size = 0.3,
              point.size = 0.2,
              outlier.shape = NA,
              fill = "test_type",
              ylab = "LT (s)",
              palette="npg",
              panel.labs = list(group = c("ST5a", "ST5b")),
              facet.by=c("group"))+
  stat_summary(geom = "point",
               fun.y = "mean",
               col = "red",
               size = 3,
               shape = 20,
               fill = "red")+ 
  stat_summary(fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = .1,
               colour = "red")+
  theme(axis.text = element_text(size = 9), 
        axis.title.y = element_text(size=10, 
                                    face="bold"))+
  scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, option = "D")
p

ggsave("Plot5a-5b.pdf", device = "pdf", width = 8, height = 8, units = c("cm"), dpi = 900)

# ST6ab
data.plot.6ab <-dplyr::filter(pyramid_infants.long.raw, study == "ST6")
p <- ggpaired(data.plot.6ab, 
              x="test_type", 
              y="rawLT", 
              width = 0.3,
              ylim = c(0, 60),
              color="black",
              line.color="gray",
              line.size = 0.3,
              point.size = 0.2,
              outlier.shape = NA,
              fill = "test_type",
              ylab = "LT (s)",
              palette="npg",
              panel.labs = list(group = c("ST6a", "ST6b")),
              facet.by=c("group"))+
              stat_summary(geom = "point",
                            fun.y = "mean",
                              col = "red",
                             size = 3,
                            shape = 20,
                             fill = "red")+ 
              stat_summary(fun.data = mean_cl_boot, 
                             geom = "errorbar", 
                             width = .1,
                             colour = "red")+
              theme(axis.text = element_text(size = 9), 
                      axis.title.y = element_text(size=10, 
                      face="bold"))+
              scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, option = "D")
  #+ geom_point(alpha = 1/10, size = .6)
p

ggsave("Plot6a-6b.pdf", device = "pdf", width = 8, height = 8, units = c("cm"), dpi = 900)

# ST7ab
data.plot.7ab <-dplyr::filter(pyramid_infants.long.raw, study == "ST7")
p <- ggpaired(data.plot.7ab, 
              x="test_type", 
              y="rawLT", 
              width = 0.3,
              ylim = c(0, 60),
              color="black",
              line.color="gray",
              line.size = 0.3,
              point.size = 0.2,
              outlier.shape = NA,
              fill = "test_type",
              ylab = "LT (s)",
              palette="npg",
              panel.labs = list(group = c("ST7a", "ST7b")),
              facet.by=c("group"))+
  stat_summary(geom = "point",
               fun.y = "mean",
               col = "red",
               size = 3,
               shape = 20,
               fill = "red")+ 
  stat_summary(fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = .1,
               colour = "red")+
  theme(axis.text = element_text(size = 9), 
        axis.title.y = element_text(size=10, 
                                    face="bold"))+
  scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, option = "D")
p

ggsave("Plot7a-7b.pdf", device = "pdf", width = 8, height = 8, units = c("cm"), dpi = 900)

