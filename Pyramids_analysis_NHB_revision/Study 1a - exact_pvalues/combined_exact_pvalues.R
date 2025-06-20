# Study 1a - Combine human data with non-human animal data

rm(list = ls(all = TRUE))

library(tidyverse)
library(MASS)
library(Deriv)
library(numDeriv)
library(DHARMa)

# Uncomment one of the two sets below to combine the data and generate relevant plots

#(1) Exact pvalues for the number of observed motifs across studies
ST1a_animals_exact_pvalues <- read.csv("ST1a_animals_exact_pvalues.csv")
ST1a_children_exact_pvalues <- read.csv("ST1a_children_exact_pvalues.csv")

#(1) Exact pvalues for the ratio (number of triadic pyramids)/ (number of triadic pyramids and triadic trees)
# ST1a_animals_exact_pvalues <- read.csv("ST1a_animals_exact_prop_pvalues.csv")
# ST1a_children_exact_pvalues <- read.csv("ST1a_children_exact_prop_pvalues.csv")



combined_exact_pvalues <- rbind(ST1a_animals_exact_pvalues, ST1a_children_exact_pvalues)

# recode & rename the groups (order) by putting together all the orders with <= 16 networks in the "Other" group
combined_exact_pvalues$order = recode(combined_exact_pvalues$order, Accipitriformes = "Other", Anseriformes = "Other", Carcharhiniformes = "Other", 
                                   Cichliformes= "Other",Decapoda = "Other",Diprotodontia = "Other",Galliformes= "Other",
                                   Lagomorpha= "Other",Osmeriformes= "Other",Perissodactyla= "Other",Proboscidea= "Other", 
                                   Psittaciformes= "Other",Salmoniformes= "Other",Sphenisciformes= "Other",Squamata= "Other",
                                   Testudines= "Other",Passeriformes ="Birds (Passeriformes)", Artiodactyla = "Ungulates (Artiodactyla)", 
                                   Carnivora = "Carnivores (Carnivora)", Children = "Human Children (Primates)", 
                                   Hymenoptera = "Social insects (Hymenoptera)",Primates = "Non-human primates (Primates)",
                                   Rodentia= "Rodents (Rodentia)")


# To select specific motifs of interest for the count tests
#combined_exact_pvalues <-  combined_exact_pvalues %>% filter(motif == "021D")
# pyr: 021D
# tre: 021U




# An example of a specific order's set of pvalues
#pv_rodents <- combined_exact_pvalues %>% filter(order == "Rodents (Rodentia)")



# Function to fit a beta distribution to data
fit_beta <- function(xvec) {
  # Validate input data
  if(any(xvec <= 0 | xvec >= 1)) {
    warning("Some values are not strictly between 0 and 1. Results may be unstable.")
    # Adjust values that are exactly 0 or 1 to be slightly inside the interval
    xvec[xvec <= 0] <- .Machine$double.eps
    xvec[xvec >= 1] <- 1 - .Machine$double.eps
  }
  
  # Method 1: Using method of moments
  mean_x <- mean(xvec)
  var_x <- var(xvec)
  
  # Calculate alpha and beta using method of moments
  mm_alpha <- mean_x * (mean_x * (1 - mean_x) / var_x - 1)
  mm_beta <- (1 - mean_x) * (mean_x * (1 - mean_x) / var_x - 1)
  
  # Method 2: Maximum likelihood estimation
  # Log-likelihood function for the beta distribution
  beta_loglik <- function(params) {
    alpha <- params[1]
    beta <- params[2]
    
    if(alpha <= 0 || beta <= 0) return(1e10)  # Return large value for invalid parameters
    
    # Negative log-likelihood
    -sum(dbeta(xvec, shape1 = alpha, shape2 = beta, log = TRUE))
  }
  
  # Initial parameter estimates from method of moments
  init_params <- c(mm_alpha, mm_beta)
  if(any(is.na(init_params) | is.infinite(init_params))) {
    init_params <- c(1, 1)  # Default to uniform distribution if MM fails
  }
  
  # Perform optimization
  mle_result <- optim(init_params, beta_loglik, method = "BFGS")
  
  # Extract the MLE parameters
  mle_alpha <- mle_result$par[1]
  mle_beta <- mle_result$par[2]
  
  # Return results
  result <- list(
    method_of_moments = list(alpha = mm_alpha, beta = mm_beta),
    maximum_likelihood = list(alpha = mle_alpha, beta = mle_beta, 
                              convergence = mle_result$convergence,
                              log_likelihood = -mle_result$value)
  )
  
  return(result)
}


# Example of fitting a beta distribution the pvalues
#fit_pv_rodents <- fit_beta(pv_rodents$pvalues)
#alpha_pv_rodents <- fit_pv_rodents$maximum_likelihood$alpha
#beta_pv_rodents <- fit_pv_rodents$maximum_likelihood$beta


# Plot of pvalue distribution across orders with fitted data density overlayed
orders <- unique(combined_exact_pvalues$order) # 12 orders
par(mfrow = c(3,4), oma = c(0, 0, 3, 0), mar = c(4, 4, 2, 1))

df_ord_beta_fits <- data.frame(order = c("test"), alpha = c(1), beta = c(2))
for (ord in orders) {
  order_study <- combined_exact_pvalues %>% filter(order == ord)
  
  order_study_pv <- order_study$pvalues
  fit_pv <- fit_beta(order_study_pv)
  
  # order_study_tpm <- order_study$tpm
  # fit_pv <- fit_beta(order_study_tpm)
  
  alpha_pv <- fit_pv$method_of_moments$alpha
  beta_pv <- fit_pv$method_of_moments$beta
  
  
  #order_study_tpm <- order_study$tpm
  
  hist(order_study_pv, freq = FALSE, main = ord, breaks = 12, xlab = "p-value")
  # hist(order_study_tpm, freq = FALSE, main = ord, breaks = 12)
  lines(seq(0,1, 0.01), dbeta(seq(0,1, 0.01), shape1 = alpha_pv, shape2 = beta_pv), col = "red")
  #df_ord_beta_fits <- rbind(df_ord_beta_fits, data.frame(order = c(ord), alpha = c(alpha_pv), beta = c(beta_pv)))
  
}
mtext("Distribution of p-values for NP across Orders", outer = TRUE, cex = 1.0, line = 1)
# plot.new()  # Create a blank plot for legend
# legend("center", 
#        legend = c("Data", "Beta fit"), 
#        col = c("black", "red"), 
#        lwd = 2, 
#        title = "Legend")
par(fig = c(0, 1, 0, 1), new = TRUE)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
legend("bottom", 
        legend = c("Data", "Beta fit"), 
        col = c("black", "red"), 
        lwd = 2, 
        xpd = TRUE,
        horiz = TRUE,
        inset = c(0, -0.1))  # Adjust inset to position below the plots
df_ord_beta_fits <- df_ord_beta_fits[-1,]
par(mfrow = c(1,1))

# hist(pv_rodents$pvalues, freq = FALSE)
# lines(seq(0,1, 0.01), dbeta(seq(0,1, 0.01), shape1 = alpha_pv_rodents, shape2 = beta_pv_rodents))



# data <- combined_exact_pvalues
# beta_params <- data %>%
#   group_by(order) %>%
#   summarise(
#     params = list(fit_beta(pvalues)$method_of_moments),
#     alpha = sapply(params, `[[`, "alpha"),
#     beta = sapply(params, `[[`, "beta")
#   )# %>% filter(!is.na(alpha), !is.na(beta)) %>% filter(!is.infinite(alpha), !is.infinite(beta))
# 
# 

