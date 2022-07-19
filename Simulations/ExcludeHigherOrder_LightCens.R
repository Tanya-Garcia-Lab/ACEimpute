# libFolder = "/nas/longleaf/home/lotspeic/R/"
# library(impeRfect, lib.loc = libFolder)
library(lme4)
library(geex)
library(magrittr)
library(tidyverse)

# Generate data set for n clusters each of size m
# with a response y, subject index id, as well as
# continuous and time-dependent covariates X
# (associated with fixed effects) and Z
# (associated with random effects)
DGM_1_random_slope = function(n = 1000, m = 3, beta, sigma = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  
  # random errors
  epsilon = rnorm(n = N, mean = 0, sd = sigma)
  
  # collect data
  sample.data = data.frame(id = rep(1:n, each = m),
                           visit = rep(1:m, n),
                           # begin constructing response with intercept and random error
                           y = beta[1] + epsilon)
  
  # fixed effects design matrix X
  X = rnorm(n = N*(p - 1), mean = 0, sd = 1) %>%
    matrix(nrow = N)
  
  # add X %*% beta[2:p] to response and add X to data
  sample.data$y = sample.data$y + X %*% beta[2:p]
  sample.data = cbind(sample.data, X)
  
  # random effects design matrix Z
  Z = rnorm(n = N, mean = 5, sd = 1) %>%
    matrix(nrow = N)
  # random slopes; 1 random slope per subject
  b = rnorm(n, 0, 1) %>%
    matrix(nrow = n) %>%
    # repeat each random slope m times using kronecker product
    kronecker(rep(1, m))
  
  # add Z * b to response and add Z to data
  sample.data$y = sample.data$y + Z * b
  sample.data = cbind(sample.data, Z)
  
  # designed to handle p fixed effects and 1 random effect
  colnames(sample.data) = c("id", "visit", "y", paste0("x", 1:(p - 1)), "z1")
  
  return(sample.data)
}

################################### UPDATED TO MISSPECIFY THE COX MODEL LATER ###################################
misspecified_cox_simulation = function(n, logHR, A, lambda = 1) {
  # Generate n iid unif(0, 1) random variable
  U = runif(n = n, min = 0, max = 1)
  
  # Create a quadratic term from the column of A
  quadA <- A*A
  
  # Calculate the linear predictor (including the interaction)
  lp <- A %*% logHR - 0.25 * quadA
  
  # Calculate H_0(T) according to Bender (2005)
  H_0.T = - log(U) * exp(- lp)
  
  # Return simulated response based on distribution character provided
  return(H_0.T / lambda) 
}
################################### UPDATED TO MISSPECIFY THE COX MODEL LATER ###################################

# Generate data set for n clusters each of size m
# with a response y, subject index id, continuous
# and time-dependent covariates X and Z, continuous and
# time-independent covariates X_t, censored covariate t, 
# as well as data (w, c, and delta) which describe the censoring
DGM_2_random_slopes = function(n = 1000, m = 3, b = NULL, 
                               beta, alpha, sigma = 1,
                               logHR, rate.t = 1, rate.c = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  
  ## STEP 1: generate X_t, time-independent covariates used to generate event time t
  # number of covariates in X_t
  p.logHR <- length(logHR)
  
  # since X_t are time-independent, the number of X_t for each subject is p.logHR for a total of n*p.logHR
  # simulate X_t as iid from Normal(0, 1)
  X_t <- rnorm(n = n*p.logHR, mean = 0, sd = 1) %>%
    matrix(ncol = p.logHR)
  colnames(X_t) = paste0("X_t", 1:p.logHR)
  
  ################################### UPDATED TO MISSPECIFY THE COX MODEL LATER ###################################
  ## STEP 2: generate t from Cox simulation using covariates X_t
  t <- misspecified_cox_simulation(n = n, logHR = logHR, A = X_t, lambda = rate.t)
  ################################### UPDATED TO MISSPECIFY THE COX MODEL LATER ###################################
  
  ## STEP 2: generate t from Cox simulation using covariates X_t
  #t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = X_t, dist = "Exponential", lambda = rate.t)
  
  ## STEP 3: generate y, X, and Z using DGM_1_random_slope()
  long.data = DGM_1_random_slope(n = n, m = m, beta = beta, sigma = sigma)
  
  # add t to long.data
  long.data$t <- rep(t, each = m)
  # time_to_event variable
  long.data$time_to_event <- with(long.data, visit - 1 - t)
  
  # add alpha * time_to_event to y
  long.data$y <- with(long.data, y + alpha * time_to_event)
  
  # add X_t to data
  long.data <- cbind(long.data, kronecker(X_t, rep(1, m)))
  
  # I HAD TROUBLE LETTING THIS LINE VARY WITH
  # DIFFERENT NUMBER OF X_t - right now it's
  # hard-coded for exactly 2 covariates X_t
  colnames(long.data)[8] = paste0("x_t", 1:p.logHR)
  
  ## STEP 4:
  # censoring mechanism
  c <- rexp(n = n, rate = rate.c)
  # censoring indicator
  delta <- as.numeric(t <= c)
  # censored value
  w <- ifelse(delta, t, c)
  
  # add censoring information to data
  long.data$c <- rep(c, each = m)
  long.data$delta <- rep(delta, each = m)
  long.data$w <- rep(w, each = m)
  
  return(long.data)
}

set.seed(114)
# number of subjects
n <- 1000
# number of visits per subject
m <- 3

## REGRESSION PARAMETERS
# log hazard ratio on X_t for simulation of t, the survival outcome
logHR <- c(1)
rate.t <- 0.5
rate.c <- 0.1 # rate.c = 0.1 --> ~20% censoring, rate.c = 0.3 --> ~40% censoring
# intercept and coefficient on x1 (in that order)
beta <- c(0, 1)
# coefficient on time_to_event for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error
sigma <- 1

long.data <- DGM_2_random_slopes(n = n, m = m, b = NULL,
                                 beta = beta, alpha = alpha, sigma = sigma,
                                 logHR = logHR, rate.t = rate.t, rate.c = rate.c)

# check censoring rate
long.data$delta %>% mean()

test.data = long.data %>%
  mutate(x_t1_2 = x_t1^2)

# check Cox coefficient estimates
survival::coxph(formula = survival::Surv(t, delta) ~ x_t1 + x_t1_2, data = test.data)

# number of simulations
num_sim <- 5
# save results in this data frame
save_res <- data.frame(sim = 1:num_sim, censored = NA,
                       beta1_fd_reml = NA, se_beta1_fd_reml = NA, alpha_fd_reml = NA, se_alpha_fd_reml = NA, sigma2_fd_reml = NA,
                       beta1_cc_reml = NA, se_beta1_cc_reml = NA, alpha_cc_reml = NA, se_alpha_cc_reml = NA, sigma2_cc_reml = NA,
                       beta1_cmi_ace = NA, se_beta1_cmi_ace = NA, alpha_cmi_ace = NA, se_alpha_cmi_ace = NA, sigma2_cmi_ace = NA, se_sigma2_cmi_ace = NA,
                       beta1_cmi_reml = NA, se_beta1_cmi_reml = NA, alpha_cmi_reml = NA, se_alpha_cmi_reml = NA, sigma2_cmi_reml = NA)

for (s in 1:num_sim) {
  ## Generate longitudinal data with 1 random slope per subject and no random intercept
  long_data <- DGM_2_random_slopes(n = n, m = m, b = NULL,
                                   beta = beta, alpha = alpha, sigma = sigma,
                                   logHR = logHR, rate.t = rate.t, rate.c = rate.c)
  
  # check censoring rate
  save_res[s, "censored"] <- 1 - mean(long_data$delta)
  
  ## DATASET 1: use REML when full data is available
  # fit mixed effects model with lmer()
  fd_reml_sum <- lmer(formula = y ~ x1 + time_to_event + (z1 - 1 | id) - 1, 
                      data = long_data) %>%
    summary()
  # save parameter estimates
  save_res[s, c("beta1_fd_reml", "alpha_fd_reml", "sigma2_fd_reml")] <- c(coef(fd_reml_sum)[, 1], fd_reml_sum$sigma^2)
  # save standard error estimates
  save_res[s, c("se_beta1_fd_reml", "se_alpha_fd_reml")] <- sqrt(diag(fd_reml_sum$vcov))
  
  ## DATASET 2: use REML with complete cases only
  # use complete case analysis
  cc_reml_sum <- long_data %>%
    dplyr::filter(delta == 1) %>%
    lmer(formula = y ~ x1 + time_to_event + (z1 - 1 | id) - 1, data = .) %>%
    summary()
  
  # save parameter estimates
  cc_reml_est <- save_res[s, c("beta1_cc_reml", "alpha_cc_reml", "sigma2_cc_reml")] <- c(coef(cc_reml_sum)[, 1], cc_reml_sum$sigma^2)
  # save standard error estimates
  save_res[s, c("se_beta1_cc_reml", "se_alpha_cc_reml")] <- sqrt(diag(cc_reml_sum$vcov))
  
  ## DATASET 3: use CMI to impute censored values of T
  # impute using only one row per subject
  imp_data <- long_data[which(long_data$visit == 1), ]
  
  # impute censored covariate t using imputeCensoRd::condl_mean_impute()
  imp_mod <- survival::coxph(formula = survival::Surv(w, delta) ~ x_t1, data = imp_data)
  # imputation model linear predictor: x_t1 + x_t2 (matches DGM)
  imp_data <- imputeCensoRd::condl_mean_impute(fit = imp_mod, obs = "w", event = "delta", addl_covar = c("x_t1"), data = imp_data)
  
  # assign imputed values based on subject id
  long_data <- imp_data %>%
    dplyr::select(id, imp) %>%
    dplyr::left_join(long_data, ., by = "id") %>%
    # recalculate time_to_event using imp
    mutate(time_to_event_imp = visit - 1 - imp)
  
  # fit mixed effects model with lmer() using imputed values of T
  cmi_reml_sum <- lmer(formula = y ~ x1 + time_to_event_imp + (z1 - 1 | id) - 1, 
                       data = long_data) %>%
    summary()
  
  # save parameter estimates
  save_res[s, c("beta1_cmi_reml", "alpha_cmi_reml", "sigma2_cmi_reml")] <- c(coef(cmi_reml_sum)[, 1], cmi_reml_sum$sigma^2)
  # save standard error estimates
  save_res[s, c("se_beta1_cmi_reml", "se_alpha_cmi_reml")] <- sqrt(diag(cmi_reml_sum$vcov))
  
  # ACE Imp estimating equation to fit mixed effects model while correcting for imputation error
  # use m_estimate() to solve estimating equation with summand given by impeRfect::eff_score_vec()
  cmi_ace_fit <- m_estimate(estFUN = eff_score_vec,
                            data = long_data, units = "id",
                            root_control = setup_root_control(start = cc_reml_est),
                            outer_args = list(response = "y",
                                              X.names = c("x1", "time_to_event_imp"),
                                              Z.names = c("z1"),
                                              # delta tells eff_score_vec which subjects are censored
                                              cens = "delta"))
  # save parameter estimates
  save_res[s, c("beta1_cmi_ace", "alpha_cmi_ace", "sigma2_cmi_ace")] <- coef(cmi_ace_fit)
  # save standard error estimates
  save_res[s, c("se_beta1_cmi_ace", "se_alpha_cmi_ace", "se_sigma2_cmi_ace")] <- sqrt(diag(vcov(cmi_ace_fit)))
  
  # show sim progress
  if (s %% 25 == 0) print(paste("Simulation", s, "complete! :D"))
  
  write.csv(save_res, "data/ExcludeHigherOrder_LightCens.csv", row.names = F)
}

# display results
save_res %>%
  dplyr::select(-sim) %>%
  dplyr::summarize_all(mean)
