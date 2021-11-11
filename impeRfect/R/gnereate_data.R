#' Generate longitudinal dataset
#'
#' Generate data set for \code{n} clusters each of size \code{m}
#' with a response Y, 1 binary, time-invariant covariate X1,
#' 1 continuous, time-variant covariate X2, and subject index id
#'
# model formula: Y = X1 + X2 + (1 | id)
# b is random, subject-specific intercept
# epsilon is random error
#' @param n sample size; default is 1000
#' @param m number of visits per subject; default is 3
#' @param b \code{n}-dimensional vector of random intercepts; if NULL, random intercept is generated from N(0, 1)
#' @param sigma standard deviation of the random error, epsilon; default is 1
#' @param p.X number of time-variant, continuous covariates
#' @param p.A number of time-invariant, binary covariates
#' @param beta p-dimensional vector of regression coefficients
# beta[1] is intercept
# if \code{p.X} > 0, beta[2:(\code{p.X} + 1)] are coefficients on X variables
# if \code{p.A} > 0, beta[(\code{p.X} + 2):p] are coefficients on A variables
#'
#' @export
generate_data = function(n = 1000, m = 3, b = NULL, sigma = 1,
                         p.X = 1, p.A = 1, beta) {
  if (length(beta) != 1 + p.X + p.A) { stop("length of beta not equal to 1 + p.X + p.A") }

  # total obs
  N = n*m
  # number of fixed effects
  p = length(beta)
  # random error sd; random errors
  sigma = 1
  epsilon = rnorm(N, 0, sigma)
  # number of random effects
  # q = 1
  # Z = matrix(1, N)
  # if (q > 1) { Z = cbind(Z, X[, 1:(q - 1)]) }
  # random effects
  if (is.null(b)) { b = matrix(rnorm(n, 0, 1), n) }
  # b = rt(n = n, df = 10)
  # response Y
  sample.data = data.frame(id = rep(1:n, each = m), Y = beta[1] + rep(b, each = m) + epsilon)
  if (p.X > 0) {
    X = rnorm(n = N*p.X, mean = 0, sd = 1) %>%
      matrix(nrow = N)
    colnames(X) = paste0("X", 1:p.X)
    sample.data$Y = sample.data$Y + X %*% beta[2:(p.X + 1)]
    sample.data = sample.data %>%
      cbind(X)
  }
  if (p.A > 0) {
    A = rbinom(n = n*p.A, size = 1, prob = 0.5) %>%
      matrix(nrow = n) %>%
      kronecker(Y = rep(1, m))
    colnames(A) = paste0("A", 1:p.A)
    sample.data$Y = sample.data$Y + A %*% beta[(p.X + 2):p]
    sample.data = sample.data %>%
      cbind(A)
  }

  return(sample.data)
}
