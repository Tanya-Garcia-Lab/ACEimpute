rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscripts_Simulations/Correctly_Specified_Imputation_Model/")

library(tidyverse)
library(latex2exp)
library(xtable)

# create table to display simulation results
create_table <- function(sim_results, true_params) {
  bias_df <- sim_results %>%
    dplyr::select(-dplyr::starts_with("se_")) %>%
    tidyr::gather("param_method", "est", -1) %>% 
    dplyr::mutate(data = ifelse(test = grepl(pattern = "_fd_", x = param_method), 
                                yes = "Full Data", 
                                no = ifelse(test = grepl(pattern = "_cc_", x = param_method), 
                                            yes = "Complete Case",
                                            no = ifelse(test = grepl(pattern = "_cmi_mi_", x = param_method),
                                                        yes = "CMI-MI",
                                                        no = "CMI"))),
                  method = paste0(toupper(sub(".*_", "", param_method)), " (", data, ")"),
                  param = sub("_.*", "", param_method),
                  truth = unlist(true_params[param]),
                  bias = est - truth)
  
  see_df <- sim_results %>%
    dplyr::select(sim, starts_with("se_")) %>%
    dplyr::mutate(sim = 1:dplyr::n()) %>% 
    tidyr::gather("param_method", "se_est", dplyr::starts_with("se_")) %>%
    dplyr::mutate(data = ifelse(test = grepl(pattern = "_fd_", x = param_method), 
                                yes = "Full Data", 
                                no = ifelse(test = grepl(pattern = "_cc_", x = param_method), 
                                            yes = "Complete Case",
                                            no = ifelse(test = grepl(pattern = "_cmi_mi_", x = param_method),
                                                        yes = "CMI-MI",
                                                        no = "CMI"))),
                  param_method = sub("se_", "", param_method))
  
  mse_df <- bias_df %>%
    left_join(see_df, by = c("sim", "data", "param_method")) %>%
    mutate(mse = bias^2 + se_est^2,
           cover = I(est + qnorm(0.025)*se_est < truth & est + qnorm(0.975)*se_est > truth))
  
  mse_df <- mse_df %>%
    dplyr::filter(!(param == "sigma2" & est < 0)) %>%
    dplyr::mutate(
      data = factor(x = data, 
                    levels = c("CMI", "CMI-MI", "Complete Case", "Full Data"),
                    labels = c("CMI", "CMI-MI", "Complete Case", "Full Data")),
      Method = factor(x = method,
                      levels = c("ACE (CMI)", "REML (CMI)", "REML (CMI-MI)", "REML (Complete Case)", "REML (Full Data)"),
                      labels = c("ACE", "CMI", "CMI-MI", "CCA", "Oracle")),
      Parameter = factor(param, 
                         levels = c("alpha", "beta1", "sigma2"),
                         labels = c("alpha", "beta1", "sigma2")))
  
  summary_df = mse_df %>% 
    group_by(Parameter, Method) %>%
    summarize(Bias = mean(bias),
              SEE = mean(se_est, na.rm = T),
              ESE = sd(est),
              MSE = mean(bias^2),
              CPr = mean(cover))
  
  return(summary_df)
}

# true parameter values
truth = list(alpha = 1, beta1 = 1, sigma2 = 1)

# create light censoring results table
results_25perc = read.csv("sim_results/full_data_reml_estimates.csv") %>%
  select(-censored)
results_25perc = read.csv("sim_results/ace_estimates_25perc.csv") %>%
  select(-censored) %>%
  left_join(results_25perc, ., by = "sim")
results_25perc = read.csv("sim_results/cmi_mi_estimates_25perc.csv") %>%
  select(-censored) %>%
  left_join(results_25perc, ., by = "sim")
table_25 = create_table(sim_results = results_25perc, true_params = truth)

# create medium censoring results table
results_50perc = read.csv("sim_results/full_data_reml_estimates.csv") %>%
  select(-censored)
results_50perc = read.csv("sim_results/ace_estimates_50perc.csv") %>%
  select(-censored) %>%
  left_join(results_50perc, ., by = "sim")
results_50perc = read.csv("sim_results/cmi_mi_estimates_50perc.csv") %>%
  select(-censored) %>%
  left_join(results_50perc, ., by = "sim")
table_50 = create_table(sim_results = results_50perc, true_params = truth) 

# create heavy censoring results table
results_75perc = read.csv("sim_results/full_data_reml_estimates.csv") %>%
  select(-censored)
results_75perc = read.csv("sim_results/ace_estimates_75perc.csv") %>%
  select(-censored) %>%
  left_join(results_75perc, ., by = "sim")
results_75perc = read.csv("sim_results/cmi_mi_estimates_75perc.csv") %>%
  select(-censored) %>%
  left_join(results_75perc, ., by = "sim")
table_75 = create_table(sim_results = results_75perc, true_params = truth) 

sink("tables/medium_and_heavy_xtable.txt")
# combine medium and heavy censoring results
rbind(table_50, table_75) %>%
  cbind(data.frame(Censoring = c("Medium", rep("", 8), 
                                 "Heavy", rep("", 8))), .) %>%
  mutate(Parameter = rep(c("alpha", rep("", 2),
                           "beta", rep("", 2),
                           "sigma2", rep("", 2)), 2)) %>%
  # print xtable for LaTeX use
  xtable(x = ., digits = 3) %>%
  print(include.rownames = F)
sink()

sink("tables/light_xtable.txt")
table_25 %>%
  cbind(data.frame(Censoring = c("Light", rep("", 8))), .) %>%
  mutate(Parameter = c("alpha", rep("", 2),
                       "beta", rep("", 2),
                       "sigma2", rep("", 2))) %>%
  # print xtable for LaTeX use
  xtable(x = ., digits = 3) %>%
  print(include.rownames = F)
sink()