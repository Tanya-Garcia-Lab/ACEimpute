# clear workspace
rm(list = ls())

setwd("/Users/kylegrosser/Documents/GitHub/ACEimpute/Manuscript_Simulations/Correctly_Specified_Imputation_Model")

library(tidyverse)
library(lme4)

# Generate data set for n clusters each of size m
# with a response y, subject index id, as well as
# continuous and time-dependent covariates 
# x1 (associated with fixed effects, beta) and 
# z1 (associated with random effects)
DGM_1_random_slope = function(n = 1000, m = 3, beta, sigma = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  
  # random errors
  epsilon = rnorm(n = N, mean = 0, sd = sigma)
  
  # generate data frame with id, visit number, and response
  sample.data = data.frame(id = rep(1:n, each = m),
                           visit = rep(1:m, n),
                           # begin constructing response with intercept and random error
                           y = beta[1] + epsilon)
  
  # fixed effects design matrix X
  X = rnorm(n = N, mean = 0, sd = 1) %>%
    matrix(nrow = N)
  
  # add X %*% beta[2:p] to response and append X to data
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
  
  # add Z * b to response and append Z to data
  sample.data$y = sample.data$y + Z * b
  sample.data = cbind(sample.data, Z)
  
  # name columns
  colnames(sample.data) = c("id", "visit", "y", paste0("x", 1:(p - 1)), "z1")
  
  return(sample.data)
}

# Generate data set for n clusters each of size m
# with a response y, subject index id, continuous
# and time-dependent covariates x1 and z1, continuous and
# time-independent covariates (X_t1, X_t2), censored covariate t, 
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
  X_t <- rnorm(n = n*(p.logHR), mean = 0, sd = 1) %>%
    matrix(ncol = p.logHR)
  colnames(X_t) = paste0("X_t", 1:p.logHR)
  
  ## STEP 2: generate t from Cox simulation using covariates X_t
  t <- ACEimpute::cox_simulation(n = n, logHR = logHR, A = X_t, dist = "Exponential", lambda = rate.t)
  
  ## STEP 3: generate y, X, and Z using DGM_1_random_slope()
  long.data = DGM_1_random_slope(n = n, m = m, beta = beta, sigma = sigma)
  
  # append t to long.data
  long.data$t <- rep(t, each = m)
  # time_to_event variable
  long.data$time_to_event <- with(long.data, visit - 1 - t)
  
  # add alpha * time_to_event to y
  long.data$y <- with(long.data, y + alpha * time_to_event)
  
  # append X_t to data
  long.data <- cbind(long.data, kronecker(X_t, rep(1, m)))
  
  # name columns
  colnames(long.data)[8:9] = paste0("X_t", 1:p.logHR)
  
  ## STEP 4:
  # censoring variable v
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

# number of subjects
n <- 1000
# number of visits per subject
m <- 3

## REGRESSION PARAMETERS
# log hazard ratios on X_t for simulation of t, the censored covariate
logHR <- c(1, -0.5)
rate.t <- 0.5
# rate.c = 0.125 --> ~25% censoring, 0.5 --> ~50% censoring, and 2 --> 75% censoring

# (intercept, coefficient on x1)
beta <- c(0, 1)
# coefficient on time_to_event for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error
sigma <- 1

# number of simulations = 5 for demonstration
num_sim <- 5

# set seed
set.seed(114)
# generate data for n*num_sim subjects with 25% censoring
simulated_datasets_25perc = DGM_2_random_slopes(n = n*num_sim, m = m, b = NULL,
                                                beta = beta, alpha = alpha, sigma = sigma,
                                                logHR = logHR, rate.t = rate.t, rate.c = 0.125)

# reset seed
set.seed(114)
# generate data for n*num_sim subjects with 50% censoring
simulated_datasets_50perc = DGM_2_random_slopes(n = n*num_sim, m = m, b = NULL,
                                                beta = beta, alpha = alpha, sigma = sigma,
                                                logHR = logHR, rate.t = rate.t, rate.c = 0.5)

# reset seed
set.seed(114)
# generate data for n*num_sim subjects with 75% censoring
simulated_datasets_75perc = DGM_2_random_slopes(n = n*num_sim, m = m, b = NULL,
                                                beta = beta, alpha = alpha, sigma = sigma,
                                                logHR = logHR, rate.t = rate.t, rate.c = 2)

# assign replicate IDs
simulated_datasets_25perc$replicate = rep(1:num_sim, each = n*m)
simulated_datasets_50perc$replicate = rep(1:num_sim, each = n*m)
simulated_datasets_75perc$replicate = rep(1:num_sim, each = n*m)

# check parameter estimates of simulated data
# simulated_datasets_25perc %>%
#   filter(visit == 1) %>%
#   survival::coxph(formula = survival::Surv(w, delta) ~ X_t1 + X_t2, data = .)
# 
# simulated_datasets_25perc %>%
#   filter(replicate == 1) %>%
#   lmer(formula = y ~ x1 + time_to_event + (z1 - 1 | id) - 1, data = .)

# see summary
summary(simulated_datasets_25perc)
summary(simulated_datasets_50perc)
summary(simulated_datasets_75perc)

# save simulated datasets
saveRDS(object = simulated_datasets_25perc, file = "sim_data/simulated_datasets_25perc_cens.rds")
saveRDS(object = simulated_datasets_50perc, file = "sim_data/simulated_datasets_50perc_cens.rds")
saveRDS(object = simulated_datasets_75perc, file = "sim_data/simulated_datasets_75perc_cens.rds")