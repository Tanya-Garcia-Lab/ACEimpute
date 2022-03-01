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
  
  # andom errors
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

  # add Z_r * b to response and add Z_r to data
  sample.data$y = sample.data$y + Z * b
  sample.data = cbind(sample.data, Z)

  colnames(sample.data) = c("id", "visit", "y", paste0("x", 1:(p - 1)), "z1")
  
  return(sample.data)
}

# Generate data set for n clusters each of size m
# with a response y, subject index id, continuous
# and time-dependent covariates X and Z, continuous and
# time-independent covariates X_t, subject age age,
# censored covariate t, as well as data (w,
# c, and delta) which describe the censoring
DGM_2_random_slopes = function(n = 1000, m = 3, b = NULL, 
                               beta, alpha, sigma = 1,
                               logHR, rate.t = 1, rate.c = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  
  ## STEP 1: generate X_t, time-independent covariates
  # number of covariates in X_t
  p.logHR <- length(logHR)
  
  # since X_t are time-independent, the number of X_t
  #  for each subject is p.logHR for a total of n*p.logHR
  # simulate X_t as iid from Normal(0, 1)
  X_t <- rnorm(n = n*p.logHR, mean = 5, sd = 1) %>%
    matrix(ncol = p.logHR)
  colnames(X_t) = paste0("X_t", 1:p.logHR)
  
  ## STEP 2: generate t from Cox simulation and age s
  t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = X_t, dist = "Exponential", lambda = rate.t)
  
  ## STEP 3: generate y, X, and Z
  # intercept and coefficients on X
  p.beta = length(beta)
  long.data = DGM_1_random_slope(n = n, m = m, beta = beta, sigma = sigma)
  
  long.data$t <- rep(t, each = m)
  # time_to_event variable
  long.data$time_to_event <- with(long.data, visit - 1 - t)
  
  # add alpha * time_to_event to y
  long.data$y <- long.data$y + alpha * (long.data$time_to_event)
  
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

# This function calculates the efficient score vector proposed by
# Grosser et al. (2022)
#
# Estimates the parameters of a linear mixed model when we do
# not want to specify a distribution for the random effects and
# when one covariate is (possibly) subject to random right-censoring
# Takes in dataframe 'data', which corresponds to data from a
# single cluster, then returns the value of efficient score vector for
# theta = (alpha, beta, sigma^2) corresponding to that cluster.
# Intended for use with m_estimate() function
# At the moment, this function only allows for one random slope and no random intercept
eff_score_vec = function(data, response,
                         X.names, Z.names = NULL, cens = NULL) {
  # number of observations in cluster
  m = nrow(data)
  
  # cens.bool = 1 if cens != NULL and data[1, cens] == 0
  cens.bool = !is.null(cens)
  if (cens.bool) { cens.bool = I(data[1, cens] == 0) }
  
  # number of fixed effects
  p = length(X.names)
  
  # longitudinal response
  y = as.matrix(data[, response])
  # fixed effects design matrix
  X = as.matrix(data[, X.names])
  
  # random effects design matrix
  if (is.null(Z.names)) {
    # if Z.names not provided, Z = 1_m
    Z = matrix(1, m, 1)
    # number of random effects is 1
    q = 1
  } else {
    # otherwise, use random effects design matrix specified by Z.names
    Z = as.matrix(data[, Z.names])
    # number of random effects is number of columns specified by Z.names
    q = length(Z.names)
  }
  
  # add leading column of 1_m if censored
  if (cens.bool) {
    Z = cbind(1, Z)
    q = q + 1
  }
  
  # The following lines are written so that this function is
  #  compatible with the function m_estimate() from the geex package
  function(theta) {
    # variance of random error
    sigma2 = theta[p+1]
    
    # calculate zeta, the fixed effects component of linear predictor
    zeta = X %*% theta[1:p]
    
    # conditional expectation of y given w, x, s, t, z
    cond.expect.y <- zeta + Z %*% MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% (y -  zeta)
    
    # conditional expectation of y'y given w, x, s, t, z
    cond.expect.y2 <- sigma2*(m - q) + sum(cond.expect.y^2)
    
    # efficient score vector for beta
    eff.score.beta <- t(X) %*% (y - cond.expect.y)/sigma2
    
    # efficient score for sigma^2
    eff.score.sigma2 <- (0.5*(sum(y^2) - cond.expect.y2) - t(zeta) %*% (y - cond.expect.y))/(sigma2^2)
    
    # return efficient vector for theta = c(beta, sigma^2)
    c(eff.score.beta, eff.score.sigma2)
  }
}

set.seed(114)
# number of subjects
n <- 1000
# number of visits per subject
m <- 3

## REGRESSION PARAMETERS
# log hazard ratio on z.t for simulation of t, the survival outcome
logHR <- c(1, -0.5)
rate.t <- 0.5
rate.c <- 5 # rate.c = 1 --> ~20% censoring, rate.c = 5 --> ~45% censoring
# intercept and coefficients on z.y (in that order)
beta <- c(1, 0.25)
# coefficient on t for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error
sigma <- 1

num_sim <- 25
save_res <- data.frame(sim = 1:num_sim, censored = NA,
                       beta1_ee = NA, se_beta1_ee = NA, alpha_ee = NA, se_alpha_ee = NA, sigma2_ee = NA,
                       beta1_lme = NA, se_beta1_lme = NA, alpha_lme = NA, se_alpha_lme = NA, sigma2_lme = NA)

for (s in 1:num_sim) {
  ## GENERATE LONGITUDINDAL DATA USING impeRfect::DGM_2
  long.data <- DGM_2_random_slopes(n = n, m = m, b = NULL,
                                   beta = beta, alpha = alpha, sigma = sigma,
                                   logHR = logHR, rate.t = rate.t, rate.c = rate.c)

  # check censoring rate
  save_res[s, "censored"] <- 1 - mean(long.data$delta)

  # impute using one row per subject
  imp_data <- long.data[which(long.data$visit == 1), ]
  
  # impute using imputeCensoRd::condl_mean_impute()
  imp_mod <- survival::coxph(formula = survival::Surv(w, delta) ~ x_t1 + x_t2, data = imp_data)

  # Impute censored covariate t
  imp_data <- imputeCensoRd::condl_mean_impute(fit = imp_mod, obs = "w", event = "delta", addl_covar = c("x_t1", "x_t2"), data = imp_data)
  long.data$imp = rep(imp_data$imp, each = m)
  
  # Recalculate time_to_event using imp
  long.data$time_to_event_imp <- with(long.data, (visit - 1) - imp)

  # analyze with impeRfect::eff_score_vector(), as done in the following code
  # get initial parameter estimates with lmer()
  lme.fit <- lmer(formula = y ~ x1 + time_to_event_imp + (z1 - 1 | id), data = long.data)
  lme.sum <- summary(lme.fit)
  lme.est <- save_res[s, c("beta1_lme", "alpha_lme", "sigma2_lme")] <- c(coef(lme.sum)[2:3, 1], lme.sum$sigma^2)
  save_res[s, c("se_beta1_lme", "se_alpha_lme")] <- sqrt(diag(vcov(lme.fit))[2:3])
  # use m_estimate() to solve estimating equations defined above
  ee.fit <- m_estimate(estFUN = eff_score_vec,
                       data = long.data, units = "id",
                       root_control = setup_root_control(start = lme.est),
                       outer_args = list(response = "y",
                                         X.names = c("x1", "time_to_event_imp"),
                                         Z.names = c("z1"),
                                         # This line is the only real change needed for censored data
                                         cens = "delta"))
  save_res[s, c("beta1_ee", "alpha_ee", "sigma2_ee")] <- coef(ee.fit)
  save_res[s, c("se_beta1_ee", "se_alpha_ee")] <- sqrt(diag(vcov(ee.fit)[-3, -3]))

  if (s %% 5 == 0) print(paste("Simulation", s, "complete!"))

  # write.csv(save_res, "~kylegrosser/Documents/GitHub/Random-Error-Imputation/Simulations/Censoring/Impute_T/Results/IT_PhRMASetting_LightCens.csv", row.names = F)
}

save_res %>%
  dplyr::select(-sim) %>%
  dplyr::summarize_all(mean)