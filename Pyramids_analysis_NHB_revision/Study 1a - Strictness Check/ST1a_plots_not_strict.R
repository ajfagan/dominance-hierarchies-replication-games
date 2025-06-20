# Study 1a Plots and statistical comparisons

rm(list=ls())

# Load packages
library(rstatix)
library(ggpubr)
library(viridis)
library(tidyverse)
library(statnet)

set.seed(123) # for reproducibility
options(digits = 2)

# Load data
load(file="combined_dom_mat_SP_not_strict.Rdata")

# recode & rename the groups (order) by putting together all the orders with <= 16 networks in the "Other" group
combined_dom_mat_SP$Group = recode(combined_dom_mat_SP$Group, Accipitriformes = "Other", Anseriformes = "Other", Carcharhiniformes = "Other", 
         Cichliformes= "Other",Decapoda = "Other",Diprotodontia = "Other",Galliformes= "Other",
         Lagomorpha= "Other",Osmeriformes= "Other",Perissodactyla= "Other",Proboscidea= "Other", 
         Psittaciformes= "Other",Salmoniformes= "Other",Sphenisciformes= "Other",Squamata= "Other",
         Testudines= "Other",Passeriformes ="Birds (Passeriformes)", Artiodactyla = "Ungulates (Artiodactyla)", 
         Carnivora = "Carnivores (Carnivora)", Children = "Human Children (Primates)", 
         Hymenoptera = "Social insects (Hymenoptera)",Primates = "Non-human primates (Primates)",
         Rodentia= "Rodents (Rodentia)")

# normalize by group
combined_dom_mat_SP = 
  do.call(data.frame,
          combined_dom_mat_SP %>%
            group_by(Group) %>%
            mutate_at(c("X021D","X021U","X021C","X030T","X030C"),
                      function(x){x/sqrt(sum(x^2))}))

# convert to long format
combined_dom_mat_SP_long <- gather(combined_dom_mat_SP, pattern, norm_z_score, c("X021D","X021U","X021C","X030T", "X030C"), factor_key=TRUE)

# arrange values order
combined_dom_mat_SP_long$Group <- factor(combined_dom_mat_SP_long$Group, levels=c("Other","Social insects (Hymenoptera)","Birds (Passeriformes)","Rodents (Rodentia)","Carnivores (Carnivora)","Ungulates (Artiodactyla)","Non-human primates (Primates)", "Human Children (Primates)"))

# Plot triadic pyramidal metric for all groups

# filter data
combined_dom_mat_SP_for_tpm <- dplyr::filter(combined_dom_mat_SP_long, pattern == "X021D")

# Calculate wilcox_test p-values, with holm correction for multiple comparisons
stat.test.tpm <- combined_dom_mat_SP_for_tpm %>% 
  group_by(Group) %>%
  wilcox_test(tpm ~ 1, paired = FALSE, mu = 0.5, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))%>%
  add_xy_position(x = "Group") 
stat.test.tpm

#add information about taxonomic group size
groups_info <- c(expression(paste("Other (", italic("n"), " = 50)")),
                 expression(paste("Social insects (Hymenoptera, ", italic("n"), " = 17)")),
                 expression(paste("Birds (Passeriformes, ", italic("n"), " = 29)")),
                 expression(paste("Rodents (Rodentia, ", italic("n"), " = 53)")),
                 expression(paste("Carnivores (Carnivora, ", italic("n"), " = 30)")),
                 expression(paste("Ungulates (Artiodactyla, ", italic("n"), " = 38)")),
                 expression(paste("Non-human primates (Primates, ", italic("n"), " = 81)")),
                 expression(paste("Human Children (Primates, ", italic("n"), " = 20)"))
                 )

# plot tpm
tpm_group_plot <- ggplot(combined_dom_mat_SP_for_tpm,aes(x=Group, y=tpm, fill =Group)) +
  coord_flip()+
  geom_hline(yintercept = 0.5, linetype = 2, color = "darkgray")+ 
  geom_boxplot(width=0.5, outlier.alpha = 0.2, outlier.size = 0.05, color="black", lwd = 0.5) +
  geom_jitter(height = 0, width=0, alpha=0.2, size = 0.05) +
  stat_summary(geom = "point",
               fun = "mean",
               col = "red",
               size = 1,
               shape = 20,
               fill = "red") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 0.5, color = "red") +
  stat_pvalue_manual(stat.test.tpm, x = "Group", label = "p.adj.signif", size = 5, y.position = 0,vjust=0.3)+
  ylab("Triadic pyramidal metric")+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
#  scale_x_discrete(limits=rev)+
  scale_x_discrete(labels=groups_info)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE)+
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(size=12, face="bold"),axis.title.y = element_blank())+
  theme(legend.position = "none")
tpm_group_plot

ggsave(tpm_group_plot, device = "eps" ,filename="Global_tpm_plot_not_strict.eps",width = 15, height=10, units = "cm", dpi = 900)
ggsave(tpm_group_plot, device = "jpg" ,filename="Global_tpm_plot_not_strict.jpg",width = 15, height=10, units = "cm", dpi = 900)

# computes effect sizes, takes time because of bootstrapping
effect.size.tpm <- combined_dom_mat_SP_for_tpm %>% 
  group_by(Group) %>%
  wilcox_effsize(tpm ~ 1, paired = FALSE, mu = 0.5, ci = TRUE, ci.type = "perc", nboot = 1000)
effect.size.tpm 

# create table summarizing tpm comparisons
stats.tpm <- cbind(stat.test.tpm [,1:2],stat.test.tpm [,5:9],effect.size.tpm[,4], effect.size.tpm[,7:9])
stats.tpm
# write csv
write.csv(stats.tpm,"tableS1.csv")


# Plot normalized z-scores for all groups

# add p-values, with holm correction for multiple comparisons
stat.test <- combined_dom_mat_SP_long %>% 
  group_by(Group,pattern) %>%
  wilcox_test(norm_z_score ~ 1, paired = FALSE, mu = 0, alternative = "two.sided") %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))%>% 
  add_xy_position(x = "Group") 
stat.test

# plot
norm_z_score_group_plot <- ggplot(combined_dom_mat_SP_long,aes(x=Group, y=norm_z_score, fill =Group)) +
  coord_flip()+
#  scale_x_discrete(limits=rev)+
  scale_x_discrete(labels=groups_info)+
  geom_hline(yintercept = 0, linetype = 2, color = "darkgray")+
  geom_boxplot(width=0.5, outlier.alpha = 0.2, outlier.size = 0.05, color="black", lwd = 0.5) +
  geom_jitter(height = 0, width=0, alpha=0.2, size = 0.05) +
  stat_summary(geom = "point",
               fun = "mean",
               col = "red",
               size = 1,
               shape = 20,
               fill = "red") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 0.5, color = "red") +
  stat_pvalue_manual(stat.test, x = "Group", label = "p.adj.signif",size = 4, y.position = -0.95,vjust=0.2)+
  facet_grid(cols=vars(pattern), scales = "free_x")+
  ylab("Normalized z-score")+
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size=12, face="bold"))+
  scale_y_continuous(limits = c(-1, 1), oob = scales::squish)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE)+
  theme(axis.title.x = element_text(size=12, face="bold"),axis.title.y = element_blank())+
  theme(legend.position = "none")
norm_z_score_group_plot

# save plot
ggsave(norm_z_score_group_plot, device = "eps" ,filename="Global_SP_plot_not_strict.eps",width = 25, height=10, units = "cm", dpi = 900)
ggsave(norm_z_score_group_plot, device = "jpg" ,filename="Global_SP_plot_not_strict.jpg",width = 25, height=10, units = "cm", dpi = 900)

# # compute effect sizes
# # Warning: takes quite some time because of bootstrapping
# effect.size <- combined_dom_mat_SP_long %>% 
#   group_by(Group, pattern) %>%
#   wilcox_effsize(norm_z_score ~ 1, paired = FALSE, mu = 0, ci = TRUE, ci.type = "perc", nboot = 1000)
# 
# # create table summarizing pairwise comparisons
# stats.z.score.chance <- cbind(stat.test[,1:3],stat.test[,6:10],effect.size[,4], effect.size[,8:10])
# 
# # write csv
# stats.z.score.chance <- apply(stats.z.score.chance,2,as.character)
# write.csv(stats.z.score.chance, "tableS3.csv")
# 
# # Pairwise comparisons between z-scores for triadic pyramids and triadic trees
# 
# # filter data
# SP_long_for_compa <-dplyr::filter(combined_dom_mat_SP_long, pattern == "X021D" | pattern == "X021U")
# 
# # run the analysis
# # compute p-values & corrected p-values
# stat.test.pairwise <- SP_long_for_compa %>% 
#   group_by(Group) %>%
#   wilcox_test(norm_z_score ~ pattern, paired = TRUE, alternative = "two.sided") %>%
#   adjust_pvalue(method = "holm") %>%
#   add_significance("p.adj",cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
# stat.test.pairwise
# 
# # computes effect sizes
# # Warning: takes time because of bootstrapping
# effect.size.pairwise <- SP_long_for_compa %>% 
#   group_by(Group) %>%
#   wilcox_effsize(norm_z_score ~ pattern, paired = TRUE, ci = TRUE, ci.type = "perc", nboot = 1000)
# effect.size.pairwise 
# 
# # create table summarizing pairwise comparisons
# stats.z.score.pairwise <-cbind(stat.test.pairwise,effect.size.pairwise[,4], effect.size.pairwise[,8:10])
# 
# # write csv
# write.csv(stats.z.score.pairwise, "tableS2.csv")

