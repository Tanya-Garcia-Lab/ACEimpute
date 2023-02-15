rm(list = ls())

set.seed(95)

setwd("~/Documents/GitHub/ACEimpute/Manuscripts_Simulations/Correctly_Specified_Imputation_Model/")

# # Run once: 
# install.packages("devtools")
# devtools::install_github(repo = "sarahlotspeich/imputeCensRd", ref = "main")
# devtools::install_github(repo = "Tanya-Garcia-Lab/ACEimpute/ACEimpute")

# Load packages
library(tidyverse)
library(imputeCensRd)
library(lme4)
library(geex)
library(ACEimpute)

# number of resamples
M = 15

# # parameter and method names
param_names = c("beta1", "alpha", "sigma2")
n_param = 3
method_names = c("cmi_mi_reml")
n_method = 1

for (f in list.files(path = "sim_data/")) {
  # # for testing
  # f = list.files(path = "simulated_data/")[1]
  
  # read in simulated dataset
  simulated_datasets = readr::read_rds(paste0("simulated_data/", f))
  cens = str_extract(string = f, pattern = "[:digit:]+")
  
  # number of simulations
  max_sim <- max(simulated_datasets$replicate)
  # change this for testing
  num_sim = 5
  
  # data.frame for storing results
  save_res = matrix(data = NA,
                    nrow = num_sim,
                    ncol = 2 + 2*n_param*n_method) %>%
    as.data.frame()

  # name columns
  colnames(save_res) = c("sim", "censored", paste0(rep(paste0(c("", "se_"), rep(param_names, each = 2)), n_method),
                                                   "_",
                                                   rep(method_names, each = 2*n_param)))

  save_res$sim = 1:num_sim
  
  for (s in 1:num_sim) {
    # # for testing
    # s = 1
    
    # grab data from one replicate
    data_s <- simulated_datasets %>%
      filter(replicate == s)
    
    # impute using only one row per subject
    imp_data <- data_s[which(data_s$visit == 1), ]
    
    # get list of M imputed datasets with sampled lambda
    imputed_data_list <- cmi_sp_resample(W = "w", Delta = "delta", Z = c("X_t1", "X_t2"), data = imp_data,
                                         M = M, trapezoidal_rule = T, surv_beyond = "d")
    
    # left join each dataset in list to longitudinal data
    imputed_data_list <- lapply(X = imputed_data_list,
                                FUN = function(x)  {
                                  x %>%
                                    dplyr::select(id, imp) %>%
                                    dplyr::left_join(data_s, ., by = "id") %>%
                                    # recalculate time_to_event using imp
                                    mutate(time_to_event_imp = visit - 1 - imp)
                                }
    )

    # collect M sets of parameter estimates
    coef_matrix <- matrix(0, nrow = M, ncol = n_param)
    var_matrix <- matrix(0, nrow = M, ncol = n_param)
    for (i in 1:M) {
      model_sum = imputed_data_list[[i]] %>%
        lmer(formula = y ~ x1 + time_to_event_imp + (z1 - 1 | id) - 1, data = .) %>%
        summary()
      coef_matrix[i, ] <- c(model_sum$coefficients[, 1], model_sum$sigma^2)
      var_matrix[i, ] <- c(model_sum$coefficients[, 2], NA)
    }
    
    # pool parameter estimates using Rubin's rules
    save_res[s, paste0(param_names, "_cmi_mi_reml")] <- colMeans(coef_matrix)
    save_res[s, paste0("se_", param_names, "_cmi_mi_reml")] <- colMeans(var_matrix) + (1 + 1/M)*apply(coef_matrix, 2, var)
    
    # show progress
    if (s %% 25 == 0) print(paste("Simulation", s, "complete! :D"))
    
    # continually save results
    write.csv(x = save_res, file = paste0("sim_results/cmi_mi_estimates_", cens, "perc.csv"), row.names = F)
  }
}