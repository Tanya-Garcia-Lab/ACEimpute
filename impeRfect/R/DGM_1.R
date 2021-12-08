#' Generate longitudinal dataset
#'
#' Generate data set for \code{n} clusters each of size \code{m}
#' with a response \code{y}, subject index \code{id}, as well as 
#' continuous and time-dependent covariates \code{z_y}
#'
#' @param n sample size; default is 1000
#' @param m number of visits per subject; default is 3
#' @param b \code{n}-dimensional vector of random intercepts; if NULL, random intercept is generated from N(0, 1)
#' @param beta vector of regression coefficients on \code{z_y}, used for generation of \code{y}
#' @param sigma standard deviation of the random error, epsilon; default is 1
#' 
#' @export
DGM_1 = function(n = 1000, m = 3, b = NULL, beta, sigma = 1) {
  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  
  # random error sd; random errors
  epsilon = rnorm(N, 0, sigma)
  
  # number of random effects
  # q = 1
  # Z = matrix(1, N)
  # if (q > 1) { Z = cbind(Z, X[, 1:(q - 1)]) }
  
  # random effects
  if (is.null(b)) { b = matrix(rnorm(n, 0, 1), n) }

  # response Y
  sample.data = data.frame(id = rep(1:n, each = m), y = beta[1] + rep(b, each = m) + epsilon)
  
  # fixed effects design matrix Z_y
  Z_y = rnorm(n = N*(p - 1), mean = 0, sd = 1) %>%
    matrix(nrow = N)
  colnames(Z_y) = paste0("z_y", 1:(p - 1))
  
  # add X %*% beta[2:p] to response and add X to data
  sample.data = sample.data %>%
    mutate(y = y + Z_y %*% beta[2:p]) %>%
    cbind(Z_y)
  
  return(sample.data)
}