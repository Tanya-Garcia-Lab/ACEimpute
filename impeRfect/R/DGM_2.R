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
#' @param beta vector of regression coefficients on \code{z_y}, used for generation of \code{y}
#' @param alpha regression coefficient on (\code{age} - \code{t}) used for generation of \code{y}
#' @param sigma standard deviation of the random error, epsilon; default is 1
#' @param logHR vector of log hazard ratios on \code{z_t}, used for generation of \code{t}
#' @param rate.t rate parameter for baseline hazard of \code{t}; default is 1
#' @param rate.c rate parameter for exponential distribution of \code{c}; default is 1
#'
#' @export
DGM_2 = function(n = 1000, m = 3, b = NULL,
                 beta, alpha, sigma = 1,
                 logHR, rate.t = 1, rate.c = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)

  ## STEP 1: generate z.t, time-independent covariates
  # number of covariates in z.t
  p.logHR <- length(logHR)

  # since z.t are time-independent, the number of z.t
  #  for each subject is p.logHR for a total of n*p.logHR
  # simulate z.t as iid from Normal(0, 1)
  z.t <- matrix(data = rnorm(n = n*p.logHR, mean = 5, sd = 1),
                ncol = p.logHR)
  colnames(z.t) = paste0("z_t", 1:p.logHR)

  ## STEP 2: generate t from Cox simulation and age s
  t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = z.t, dist = "Exponential", lambda = rate.t)

  ## STEP 3: generate y, the longitudinal response, and z.y, the time-dependent covariates
  # intercept and coefficients on z.y
  p.beta = length(beta)
  long.data = impeRfect::DGM_1(n = n, m = m, b = b, beta = beta, sigma = sigma)

  long.data$t <- rep(t, each = m)
  # time_to_event variable
  long.data$time_to_event <- with(long.data, visit - 1 - t)

  # add alpha * (age - t) to Y
  long.data$y <- long.data$y + alpha * (long.data$time_to_event)

  # add time-independent covariates to data
  long.data <- cbind(long.data, kronecker(z.t, rep(1, m)))

  # I HAD TROUBLE LETTING THIS LINE VARY WITH
  # DIFFERENT NUMBER OF z_t - right now it's
  # hard-coded for exactly 2 covariates z_t
  colnames(long.data)[8:9] = paste0("z_t", 1:p.logHR)

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
