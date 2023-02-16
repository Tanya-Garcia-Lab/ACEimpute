# clear workspace
rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscripts_Simulations/Correctly_Specified_Imputation_Model")

# # Run once: 
# install.packages("devtools")
# devtools::install_github(repo = "sarahlotspeich/imputeCensRd", ref = "main")
# devtools::install_github(repo = "Tanya-Garcia-Lab/ACEimpute/ACEimpute")

# Load packages
library(tidyverse)
library(lme4)
library(geex)
library(imputeCensRd)
library(ACEimpute)

# parameter and method names
param_names = c("beta1", "alpha", "sigma2")
n_param = 3
method_names = c("cmi_correct_ace")
n_method = 1

for (f in list.files(path = "sim_data/")) {
  # # for testing
  # f = list.files(path = "simulated_data/")[1]
  
  # read in simulated datasets
  simulated_datasets = readr::read_rds(paste0("sim_data/", f))
  cens = str_extract(string = f, pattern = "[:digit:]+")
  
  # number of simulations
  max_sim <- max(simulated_datasets$replicate)
  # change this for testing
  num_sim = max_sim
  
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
    
    ## get initial estimates from complete case analysis
    # use complete case analysis with REML
    cc_reml_sum <- data_s %>%
      dplyr::filter(delta == 1) %>%
      lmer(formula = y ~ x1 + time_to_event + (z1 - 1 | id) - 1, 
           data = .) %>%
      summary()
    
    cc_reml_est_s = c(coef(cc_reml_sum)[, 1], cc_reml_sum$sigma^2)
    
    # impute using only one row per subject
    imp_data <- data_s[which(data_s$visit == 1), ]
    
    # impute censored covariate t
    # imputation model linear predictor: x_t1 + x_t2 (matches DGM)
    imp_data <- cmi_sp(W = "w", Delta = "delta", Z = c("X_t1", "X_t2"), data = imp_data, 
                       trapezoidal_rule = T, surv_beyond = "d")
    
    # assign imputed values based on subject id
    data_s <- imp_data[[1]] %>%
      dplyr::select(id, imp) %>%
      dplyr::left_join(data_s, ., by = "id") %>%
      # recalculate time_to_event using imp
      mutate(time_to_event_imp = visit - 1 - imp)
    
    # use m_estimate() to solve estimating equation with summand given by impeRfect::eff_score_vec()
    cmi_ace_fit <- m_estimate(estFUN = eff_score_vec,
                              data = data_s, units = "id",
                              root_control = setup_root_control(start = cc_reml_est_s),
                              outer_args = list(response = "y",
                                                X.names = c("x1", "time_to_event_imp"),
                                                Z.names = c("z1"),
                                                # delta tells eff_score_vec which subjects are censored
                                                cens = "delta"))
    # save parameter estimates
    save_res[s, paste0(param_names, "_cmi_correct_ace")] <- coef(cmi_ace_fit)
    # save standard error estimates
    save_res[s,  paste0("se_", param_names, "_cmi_correct_ace")] <- sqrt(diag(vcov(cmi_ace_fit)))
    
    # show progress
    if (s %% 25 == 0) print(paste("Simulation", s, "complete! :D"))
    
    # continually save results
    write.csv(x = save_res, file = paste0("sim_results/ace_estimates_", cens, "perc.csv"), row.names = F)
  }
}