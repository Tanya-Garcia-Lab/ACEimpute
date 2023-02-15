# clear workspace
rm(list = ls())

# set working and library directories
setwd("~/Documents/GitHub/ACEimpute/Manuscripts_Simulations/Correctly_Specified_Imputation_Model/")

# Load packages
library(tidyverse)
library(lme4)

# read in simulated datasets
simulated_datasets = readr::read_rds("simulated_data/simulated_datasets_25perc_cens.rds")

# number of simulations
max_sim <- max(simulated_datasets$replicate)
num_sim = 5

param_names = c("beta1", "alpha", "sigma2")
n_param = 3
method_names = c("fd_reml")
n_method = 1

save_res = matrix(data = NA, 
                  nrow = num_sim, 
                  ncol = 2 + 2*n_param*n_method) %>%
  as.data.frame()

colnames(save_res) = c("sim", "censored", paste0(rep(paste0(c("", "se_"), rep(param_names, each = 2)), n_method), 
                                                 "_",
                                                 rep(method_names, each = 2*n_param)))

save_res$sim = 1:num_sim

head(save_res)

for (s in 1:num_sim) {
  ## Filter to replicate == s only
  data_s <- simulated_datasets %>%
    filter(replicate == s)
  
  ## use REML when full data is available
  # fit mixed effects model with lmer()
  fd_reml_sum <- lmer(formula = y ~ x1 + time_to_event + (z1 - 1 | id) - 1, 
                      data = data_s, REML = T) %>%
    summary()
  
  # save parameter estimates
  save_res[s, paste0(param_names, "_fd_reml")] <- c(coef(fd_reml_sum)[, 1], fd_reml_sum$sigma^2)
  # save standard error estimates
  save_res[s, paste0("se_", param_names, "_fd_reml")] <- c(sqrt(diag(fd_reml_sum$vcov)), NA)
  
  if (s %% 25 == 0) print(paste("Simulation", s, "complete! :D"))
  
  write.csv(x = save_res, file = "sim_results/full_data_reml_estimates.csv", row.names = F)
}