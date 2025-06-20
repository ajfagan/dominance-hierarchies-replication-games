packages = c("ggplot2",
             "dplyr",
             "tidyverse",
             "statnet",
             "igraph",
             "rstatix",
             "coin",
             "ggpubr",
             "viridis",
             "Hmisc",
             "rotl",
             #"datelife",
             "MCMCglmm",
             "phytools",
             "phylobase",
             "phylosignal",
             "devtools",
             "ggthemes",
             "cowplot",
             "rcompanion",
             "emmeans",
             "afex",
             "bruceR",
             "pairwiseCI",
             "effectsize",
             "DHARMa"
)
#datelife not avaialble on this version of R
## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)


packages[!(packages %in% (.packages())) ]
