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
  X_t <- rnorm(n = n*p.logHR, mean = 5, sd = 1) %>%
    matrix(ncol = p.logHR)
  colnames(X_t) = paste0("X_t", 1:p.logHR)
  
  ## STEP 2: generate t from Cox simulation using covariates X_t
  t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = X_t, dist = "Exponential", lambda = rate.t)
  
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
  colnames(long.data)[8:9] = paste0("x_t", 1:p.logHR)
  
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
logHR <- c(1, -0.5)
rate.t <- 0.5
rate.c <- 5 # rate.c = 1 --> ~20% censoring, rate.c = 5 --> ~45% censoring
# intercept and coefficient on x1 (in that order)
beta <- c(1, 0.25)
# coefficient on time_to_event for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error
sigma <- 1

# number of simulations
num_sim <- 1000
# save results in this data frame
save_res <- data.frame(sim = 1:num_sim, censored = NA,
                       beta1_ee = NA, se_beta1_ee = NA, alpha_ee = NA, se_alpha_ee = NA, sigma2_ee = NA,
                       beta1_lme = NA, se_beta1_lme = NA, alpha_lme = NA, se_alpha_lme = NA, sigma2_lme = NA)

for (s in 1:num_sim) {
  ## Generate longitudinal data with 1 random slope per subject and no random intercept
  long.data <- DGM_2_random_slopes(n = n, m = m, b = NULL,
                                   beta = beta, alpha = alpha, sigma = sigma,
                                   logHR = logHR, rate.t = rate.t, rate.c = rate.c)
  
  # check censoring rate
  save_res[s, "censored"] <- 1 - mean(long.data$delta)
  
  # impute using only one row per subject
  imp_data <- long.data[which(long.data$visit == 1), ]
  
  # impute censored covariate t using imputeCensoRd::condl_mean_impute()
  imp_mod <- survival::coxph(formula = survival::Surv(w, delta) ~ x_t1 + x_t2, data = imp_data)
  imp_data <- imputeCensoRd::condl_mean_impute(fit = imp_mod, obs = "w", event = "delta", addl_covar = c("x_t1", "x_t2"), data = imp_data)
  
  # repeat imputed covariate imp m times for each subject
  long.data$imp = rep(imp_data$imp, each = m)
  
  # recalculate time_to_event using imp
  long.data$time_to_event_imp <- with(long.data, (visit - 1) - imp)
  
  # get initial parameter estimates with lmer()
  lme.fit <- lmer(formula = y ~ x1 + time_to_event_imp + (z1 - 1 | id), data = long.data)
  lme.sum <- summary(lme.fit)
  # save parameter estimates (not intercept)
  lme.est <- save_res[s, c("beta1_lme", "alpha_lme", "sigma2_lme")] <- c(coef(lme.sum)[2:3, 1], lme.sum$sigma^2)
  # save standard error estimates (not intercept)
  save_res[s, c("se_beta1_lme", "se_alpha_lme")] <- sqrt(diag(vcov(lme.fit))[2:3])
  
  # use m_estimate() to solve estimating equation with summand given by impeRfect::eff_score_vec()
  ee.fit <- m_estimate(estFUN = eff_score_vec,
                       data = long.data, units = "id",
                       root_control = setup_root_control(start = lme.est),
                       outer_args = list(response = "y",
                                         X.names = c("x1", "time_to_event_imp"),
                                         Z.names = c("z1"),
                                         # This line is the only real change needed for censored data
                                         cens = "delta"))
  # save parameter estimates
  save_res[s, c("beta1_ee", "alpha_ee", "sigma2_ee")] <- coef(ee.fit)
  # save standard error estimates
  save_res[s, c("se_beta1_ee", "se_alpha_ee")] <- sqrt(diag(vcov(ee.fit)[-3, -3]))
  
  if (s %% 25 == 0) print(paste("Simulation", s, "complete!"))
  
  write.csv(save_res, "~ImputeT_RandomSlopes_HeavyCens.csv", row.names = F)
}

save_res %>%
  dplyr::select(-sim) %>%
  dplyr::summarize_all(mean)