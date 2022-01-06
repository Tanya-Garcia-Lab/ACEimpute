#' Generate longitudinal dataset with one censored covariate
#'
#' Generate data set for \code{n} clusters each of size \code{m}
#' with a response \code{y}, subject index \code{id}, continuous
#' and time-dependent covariates \code{z_y}, continuous and
#' time-independent covariates \code{z_t}, subject age \code{age},
#' censored covariate \code{t}, as well as data (\code{w},
#' \code{c}, and \code{delta}) which describe the censoring
#'
#' @param n sample size; default is 1000
#' @param m number of visits per subject; default is 3
#' @param b \code{n}-dimensional vector of random intercepts; if NULL, random intercept is generated from N(0, 1)
#' @param logHR vector of log hazard ratios on \code{z_t}, used for generation of \code{t}
#' @param beta vector of regression coefficients on \code{z_y}, used for generation of \code{y}
#' @param alpha regression coefficient on (\code{age} - \code{t}) used for generation of \code{y}
#' @param sigma standard deviation of the random error, epsilon; default is 1
#' @param min.age minimum for simulation of age, i.e., age ~ Uniform(min.age, max.age)
#' @param max.age maximum for simulation of age, i.e., age ~ Uniform(min.age, max.age)
#'
#' @export
DGM_2 = function(n = 1000, m = 3, b = NULL,
                 logHR, lambda, beta, alpha, sigma = 1,
                 min.c = 0, max.c = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)

  ## STEP 1: generate z.t, time-independent covariates

  # number of covariates in z.t
  p.logHR <- length(logHR)
  # since z.t are time-independent, the number of z.t
  # for each subject is p.logHR for a total of n*p.logHR

  # simulate z.t as iid from Normal(0, 1)
  z.t <- matrix(data = rnorm(n = n*p.logHR, mean = 5, sd = 1),
                ncol = p.logHR)
  colnames(z.t) = paste0("z_t", 1:p.logHR)

  ## STEP 2: generate t from Cox simulation and age s
  t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = z.t, dist = "Exponential", lambda = lambda)
  age.init <- runif(n = n, min = 30, max = 50)
  # max.age <- age.init + m - 1

  ## STEP 3: generate y, the longitudinal response, and z.y, the time-dependent covariates
  # intercept and coefficients on z.y
  p.beta = length(beta)
  long.data = impeRfect::DGM_1(n = n, m = m, b = b, beta = beta, sigma = sigma)

  # increase age by 1 for each visit
  long.data$age <- rep(age.init, each = m) + (long.data$visit - 1)
  # long.data$max.age <- rep(max.age, each = m)
  long.data$t <- rep(t, each = m)

  # add alpha * (age - t) to Y
  long.data$y <- long.data$y + alpha * (long.data$age - long.data$t)

  # add time-independent covariates to data
  long.data <- cbind(long.data, kronecker(z.t, rep(1, m)))

  # I HAD TROUBLE LETTING THIS LINE VARY WITH
  # DIFFERENT NUMBER OF z_t - right now it's
  # hard-coded for exactly 2 covariates z_t
  colnames(long.data)[8:9] = paste0("z_t", 1:p.logHR)

  ## STEP 4: generate c, the censoring mechanism
  c <- runif(n = n, min = min.c, max = max.c)
  delta <- as.numeric(t <= c)
  w <- ifelse(delta, t, c)

  # add censoring information to data
  long.data$c <- rep(c, each = m)
  long.data$delta <- rep(delta, each = m)
  long.data$w <- rep(w, each = m)

  return(long.data)
}

# # inspect data frame
# head(long.data)
# # censoring rate
# 1 - mean(long.data$delta)
# # another sanity check
# # truth = (0.5, 2, -1, 1.5)
# lm(formula = y ~ z_y1 + z_y2 + t, data = long.data) %>% coef()
